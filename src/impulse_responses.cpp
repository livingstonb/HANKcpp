#include <impulse_responses.h>
#include <parameters.h>
#include <model.h>
#include <steady_state.h>
#include <math.h>
#include <iostream>
#include <functional>
#include <utilities.h>

using namespace std::placeholders;

namespace {
	VectorXr get_AR1_path_logs(int T, double x0, double eps, double rho, const VectorXr& deltas, int nback);

	VectorXr get_AR1_path_levels(int T, double x0, double eps, double rho, const VectorXr& deltas, int nback);
}

void TransShock::setup() {
	if ( type == ShockType::none ) {
		std::cerr << "Shock type has not been set\n";
		throw 0;
	}

	if ( size == value_not_set ) {
		if ( type == ShockType::tfp_Y )
			size = -0.02;
		else if ( type == ShockType::monetary )
			size = 0.0025 / 4.0;
		else if ( type == ShockType::riskaver )
			size = 0.001;
	}

	if ( pers == value_not_set )
		pers = exp(-0.5);
}

IRF::IRF(const Parameters& p_, const Model& model_, const SteadyState& iss_) : p(p_), model(model_), iss(iss_) {}

void IRF::setup() {
	npricetrans = (p.nocc + 2) * Ttrans - 1;

	if ( (p.capadjcost > 0.0) | (p.invadjcost > 0.0) )
		npricetrans += Ttrans;

	if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
		++npricetrans;

	shock.setup();

	if ( permanentShock & (shock.type == ShockType::riskaver) )
		trans_equm.ss_riskaver = p.riskaver + shock.size;
	else
		trans_equm.ss_riskaver = p.riskaver;

	construct_delta_trans_vectors();
	set_shock_paths();
}

void IRF::construct_delta_trans_vectors() {
	double lb = (Ttrans + 1) * (deltatranstot - Ttrans * deltatransmin) / (Ttrans * (deltatranstot - deltatransmin));
	double la = (Ttrans + 1) * deltatransmin - deltatransmin * lb;

	double lit;
	deltatransvec = VectorXr(Ttrans);
	cumdeltatrans = VectorXr(Ttrans);
	for (int it=0; it<Ttrans; ++it) {
		lit = static_cast<double>(it) / (Ttrans + 1);
		cumdeltatrans(it) = la * lit / (1.0 - lb * lit);
	}

	deltatransvec(seq(1, Ttrans-1)) = cumdeltatrans(seq(1, Ttrans-1)) - cumdeltatrans(seq(0, Ttrans-2));
	deltatransvec(0) = cumdeltatrans(0);
}

void IRF::compute() {
	// Final steady state
	SteadyState finalss = iss;
	if ( !permanentShock ) {
		finalss.mode = SteadyState::SSType::final;
	}

	// Guess log deviations from steady state
	VectorXr xguess(npricetrans);
	if ( !permanentShock ) {
		xguess = VectorXr::Zero(npricetrans);

		if ( shock.type == ShockType::tfp_Y ) {
			VectorXr leveldev = fmax(p.capadjcost, 1.0) * trans_equm.tfp_Y / iss.tfp_Y;
			xguess(seq(0, Ttrans-1)) = leveldev.array().log();
		}
	}

	if ( solver == SolverType::broyden ) {
		std::function<void(int, const hank_float_type*, hank_float_type*)>
			obj_fn = std::bind(&IRF::transition_fcn, *this, _1, _2, _3);

		hank_float_type* z = new hank_float_type[npricetrans];
		for (int i=0; i<npricetrans; ++i)
			z[i] = 0.0;

		obj_fn(npricetrans, xguess.data(), z);

		delete[] z;
	}
	else {
		std::cerr << "Must select Broyden solver\n";
		throw 0;
	}
}

void IRF::transition_fcn(int n, const hank_float_type *x, hank_float_type *z) {
	// Guesses
	trans_equm.output = VectorXr(Ttrans);

	for (int i=0; i<Ttrans; ++i)
		trans_equm.output(i) = iss.output * x[i];

}

void IRF::set_shock_paths() {
	if ( shock.type == ShockType::tfp_Y )
		trans_equm.tfp_Y = get_AR1_path_logs(Ttrans, iss.tfp_Y, shock.size, shock.pers, deltatransvec, nendtrans);
	else if ( shock.type == ShockType::monetary ) {
		trans_equm.mpshock = get_AR1_path_levels(Ttrans, iss.mpshock, shock.size, shock.pers, deltatransvec, nendtrans);
		trans_equm.mpshock(seq(1, Ttrans-1)) = VectorXr::Constant(Ttrans-1, 0.0);
	}
	else if ( shock.type == ShockType::riskaver ) {
		if ( permanentShock )
			trans_equm.riskaver = VectorXr::Constant(Ttrans, trans_equm.ss_riskaver);
		else
			trans_equm.riskaver = get_AR1_path_logs(Ttrans, p.riskaver, shock.size, shock.pers, deltatransvec, nendtrans);
	}
}

namespace {
	VectorXr get_AR1_path_logs(int T, double x0, double eps, double rho, const VectorXr& deltas, int nback) {
		VectorXr y(T);
		double rho_it;

		y(0) = x0 * exp(eps);
		for (int it=1; it<T-nback; ++it) {
			rho_it = pow(rho, deltas(it-1));
			y(it) = pow(x0, 1.0 - rho_it) * pow(y(it-1), rho_it);
		}

		for (int it=T-nback; it<T; ++it)
			y(it) = x0;

		return y;
	}

	VectorXr get_AR1_path_levels(int T, double x0, double eps, double rho, const VectorXr& deltas, int nback) {
		VectorXr y(T);
		double rho_it;

		y(0) = x0 + eps;
		for (int it=1; it<T-nback; ++it) {
			rho_it = pow(rho, deltas(it-1));
			y(it) = x0 * (1.0 - rho_it) + y(it-1) * rho_it;
		}
		
		for (int it=T-nback; it<T; ++it)
			y(it) = x0;

		return y;
	}
}