#include <impulse_responses.h>
#include <parameters.h>
#include <model.h>
#include <steady_state.h>
#include <math.h>
#include <iostream>
#include <functional>
#include <utilities.h>
#include <hank_numerics.h>

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

void TransEquilibrium::set_array_sizes(const Parameters& p, int T) {
	tfp_Y = VectorXr::Zero(T);
	mpshock = VectorXr::Zero(T);
	riskaver = VectorXr::Zero(T);
	output = VectorXr::Zero(T);
	pi = VectorXr::Zero(T);
	rb = VectorXr::Zero(T);
	rnom = VectorXr::Zero(T);
	qcapital = VectorXr::Zero(T);
	pricelev = VectorXr::Zero(T);
	labor_occ = MatrixXr::Zero(T, p.nocc);
	priceadjust = VectorXr::Zero(T);
	pidot = VectorXr::Zero(T);
	logydot = VectorXr::Zero(T);
	elast = VectorXr::Zero(T);
	price_W = VectorXr::Zero(T);
}

Equilibrium::Equilibrium(const Parameters& p, const SteadyState& ss) {
	rb = p.rb;
	pi = p.pi;
	rnom = rb - pi;
	output = ss.output;
	labor_occ = ss.labor_occ;
	qcapital = 1.0;
}

IRF::IRF(const Parameters& p_, const Model& model_, const SteadyState& iss_) : p(p_), model(model_), iss(iss_) {}

void IRF::setup() {
	npricetrans = (p.nocc + 2) * Ttrans - 1;

	if ( (p.capadjcost > 0.0) | (p.invadjcost > 0.0) )
		npricetrans += Ttrans;

	if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
		++npricetrans;

	shock.setup();

	trans_equm.set_array_sizes(p, Ttrans);

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
	initial_equm = Equilibrium(p, iss);
	
	if ( !permanentShock ) {
		finalss.mode = SteadyState::SSType::final;
		final_equm = Equilibrium(p, iss);
	}
	else {
		// Compute final steady state
	}

	// Guess log deviations from steady state
	VectorXr xguess = VectorXr::Zero(npricetrans);
	if ( !permanentShock ) {
		if ( shock.type == ShockType::tfp_Y ) {
			// Guess output changes so capital constant
			VectorXr leveldev = fmax(p.capadjcost, 1.0) * trans_equm.tfp_Y / iss.tfp_Y;
			xguess(seq(0, Ttrans-1)) = leveldev.array().log();
		}
	}

	if ( solver == SolverType::broyden ) {
		std::function<void(int, const hank_float_type*, hank_float_type*)>
			obj_fn = std::bind(&IRF::transition_fcn, *this, _1, _2, _3);

		hank_float_type* z = new hank_float_type[npricetrans];
		HankUtilities::fillarr(z, 0.0, npricetrans);

		obj_fn(npricetrans, xguess.data(), z);

		double jacstepprice = 1.0e-6;
		hank_float_type* fjac = new hank_float_type[npricetrans * npricetrans];
		HankNumerics::jacobian_square(obj_fn, npricetrans, xguess.data(), z, fjac, jacstepprice);

		// HankUtilities::printvec(fjac, npricetrans * npricetrans);

		delete[] z;
		delete[] fjac;
	}
	else {
		std::cerr << "Must select Broyden solver\n";
		throw 0;
	}
}

void IRF::transition_fcn(int n, const hank_float_type *x, hank_float_type *z) {
	// Guesses
	bool set_rb = (shock.type != ShockType::monetary) | solveFlexPriceTransitions;
	bool set_pi = (shock.type == ShockType::monetary) & (!solveFlexPriceTransitions);

	assert( !(set_pi & set_rb) );
	assert( set_rb | set_pi );

	for (int it=0; it<Ttrans; ++it) {
		trans_equm.output(it) = initial_equm.output * exp(x[it]);

		// Guess for rb
		if ( set_rb & (it == Ttrans - 1) ) {
			trans_equm.rb(Ttrans-1) = final_equm.rb;
		}
		else if ( set_rb ) {
			trans_equm.rb(it) = initial_equm.rb + x[Ttrans + it];
		}

		// Guess for pi
		if ( set_pi & (it == Ttrans - 1) ) {
			trans_equm.pi(Ttrans-1) = final_equm.pi;
		}
		else if ( set_pi ) {
			trans_equm.pi(it) = initial_equm.pi + x[Ttrans + it];
		}

		// Guess for labor
		for (int io=0; io<p.nocc; ++io) {
			trans_equm.labor_occ(it, io) = initial_equm.labor_occ[io];
		}

		// Guess for capital
		for (int it=0; it<Ttrans; ++it) {
			if ( (p.capadjcost > 0) | (p.invadjcost > 0) ) {
				trans_equm.qcapital(it) = initial_equm.qcapital + exp(x[Ttrans * (p.nocc + 2) + it - 1]);
			}
			else
				trans_equm.qcapital(it) = 1.0;
		}

		// Guess for tax increase
		hank_float_type initlumpincr = 0;
		if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
			initlumpincr = x[Ttrans * (p.nocc + 2) - 1];

		// Inflation and nominal interest rate
		double ygap =  log(trans_equm.output(it) / initial_equm.output);
		if ( p.taylor.use_feedback_rule ) {
			if ( set_rb )
				trans_equm.rnom(it) = (initial_equm.rnom - p.taylor.coeff_pi * trans_equm.pi(it) + ygap + trans_equm.mpshock(it)) / (1.0 - p.taylor.coeff_pi);
			else if ( set_pi )
				trans_equm.rnom(it) = initial_equm.rnom + p.taylor.coeff_pi * trans_equm.pi(it) + p.taylor.coeff_y * ygap + trans_equm.mpshock(it);
		}
		else {
			// Partial adjustment rule
			if ( it == 0 )
				trans_equm.rnom(it) = initial_equm.rnom + trans_equm.mpshock(it);
			else
				trans_equm.rnom(it) = (trans_equm.rnom(it-1) + deltatransvec(it-1) * p.taylor.pers
					* (initial_equm.rnom + p.taylor.coeff_pi * trans_equm.pi(it) + p.taylor.coeff_y * ygap + trans_equm.mpshock(it))) / (1.0 + deltatransvec(it-1) * p.taylor.pers); 
		}

		if ( set_rb )
			trans_equm.pi(it) = trans_equm.rnom(it) - trans_equm.rb(it);
		else
			trans_equm.rb(it) = trans_equm.rnom(it) - trans_equm.pi(it);

		// Guess price level
		if ( it == 0 )
			trans_equm.pricelev(it) = 1.0;
		else
			trans_equm.pricelev(it) = trans_equm.pricelev(it-1) / (1.0 - deltatransvec(it-1) * trans_equm.pi(it-1));

		// Guess price adj cost
		if ( solveFlexPriceTransitions )
			trans_equm.priceadjust(it) = 0.0;
		else
			trans_equm.priceadjust(it) = (p.priceadjcost / 2.0) * pow(trans_equm.pi(it), 2) * trans_equm.output(it);

		// Guess wholesale price
		if ( solveFlexPriceTransitions )
			trans_equm.price_W(it) = 1.0 - 1.0 / trans_equm.elast(it);
	}

	// Guess inflation and output growth
	trans_equm.pidot(Ttrans-1) = 0;
	trans_equm.logydot(Ttrans-1) = 0;

	for (int it=Ttrans-2; it>0; --it) {
		trans_equm.pidot(it) = (trans_equm.pi(it+1) - trans_equm.pi(it)) / deltatransvec(it);
		trans_equm.logydot(it) = (log(trans_equm.output(it+1)) - log(trans_equm.output(it))) / deltatransvec(it);
	}



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