#include <impulse_responses.h>
#include <parameters.h>
#include <model.h>
#include <equilibrium.h>
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

IRF::IRF(const Parameters& p_, const Model& model_, const EquilibriumElement& iss_) : p(p_), model(model_), initial_equm(iss_) {}

void IRF::setup() {
	npricetrans = (p.nocc + 2) * Ttrans - 1;

	if ( (p.capadjcost > 0.0) | (p.invadjcost > 0.0) )
		npricetrans += Ttrans;

	if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
		++npricetrans;

	shock.setup();

	trans_equm.reset(Ttrans);

	VectorXr riskaver;
	if ( permanentShock & (shock.type == ShockType::riskaver) )
		riskaver = VectorXr::Constant(Ttrans, p.riskaver + shock.size);
	else
		riskaver = VectorXr::Constant(Ttrans, p.riskaver);

	for (int it=0; it<Ttrans; ++it) {
		trans_equm[it].riskaver = riskaver[it];
	}

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
	if ( !permanentShock ) {
		final_equm = initial_equm;

	}
	else {
		// Compute final steady state
		// final_equm.create_final_steady_state();
	}

	// Guess log deviations from steady state
	VectorXr xguess = VectorXr::Zero(npricetrans);
	if ( !permanentShock ) {
		if ( shock.type == ShockType::tfp_Y ) {
			// Guess output changes so capital constant
			for (int it=0; it<Ttrans; ++it) {
				xguess(it) = log(fmax(p.capadjcost, 1.0) * trans_equm[it].tfp_Y / initial_equm.tfp_Y);
			}
		}
	}

	if ( solver == SolverType::broyden ) {
		std::function<void(int, const hank_float_type*, hank_float_type*)>
			obj_fn = std::bind(&IRF::transition_fcn, this, _1, _2, _3);

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
	make_transition_guesses(n, x, z);
}

void IRF::make_transition_guesses(int n, const hank_float_type *x, hank_float_type *z) {
	// Guesses
	bool set_rb = (shock.type != ShockType::monetary) | solveFlexPriceTransitions;
	bool set_pi = (shock.type == ShockType::monetary) & (!solveFlexPriceTransitions);

	assert( !(set_pi & set_rb) );
	assert( set_rb | set_pi );

	for (int it=0; it<Ttrans; ++it) {
		trans_equm[it].output = initial_equm.output * exp(x[it]);

		// Guess for rb
		if ( set_rb & (it == Ttrans - 1) ) {
			trans_equm[Ttrans-1].rb = final_equm.rb;
		}
		else if ( set_rb ) {
			trans_equm[it].rb = initial_equm.rb + x[Ttrans + it];
		}

		// Guess for pi
		if ( set_pi & (it == Ttrans - 1) ) {
			trans_equm[Ttrans-1].rb = final_equm.pi;
		}
		else if ( set_pi ) {
			trans_equm[it].pi = initial_equm.pi + x[Ttrans + it];
		}

		// Guess for labor
		trans_equm[it].labor_occ = initial_equm.labor_occ;

		// Guess for capital
		for (int it=0; it<Ttrans; ++it) {
			if ( (p.capadjcost > 0) | (p.invadjcost > 0) ) {
				trans_equm[it].qcapital = initial_equm.qcapital + exp(x[Ttrans * (p.nocc + 2) + it - 1]);
			}
			else
				trans_equm[it].qcapital = 1.0;
		}

		// Guess for tax increase
		hank_float_type initlumpincr = 0;
		if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
			initlumpincr = x[Ttrans * (p.nocc + 2) - 1];

		// Inflation and nominal interest rate
		double ygap =  log(trans_equm[it].output / initial_equm.output);
		if ( p.taylor.use_feedback_rule ) {
			if ( set_rb )
				trans_equm[it].rnom = (initial_equm.rnom - p.taylor.coeff_pi * trans_equm[it].pi + ygap + trans_equm[it].mpshock) / (1.0 - p.taylor.coeff_pi);
			else if ( set_pi )
				trans_equm[it].rnom = initial_equm.rnom + p.taylor.coeff_pi * trans_equm[it].pi + p.taylor.coeff_y * ygap + trans_equm[it].mpshock;
		}
		else {
			// Partial adjustment rule
			if ( it == 0 )
				trans_equm[it].rnom = initial_equm.rnom + trans_equm[it].mpshock;
			else
				trans_equm[it].rnom = (trans_equm[it-1].rnom + deltatransvec(it-1) * p.taylor.pers
					* (initial_equm.rnom + p.taylor.coeff_pi * trans_equm[it].pi + p.taylor.coeff_y * ygap + trans_equm[it].mpshock)) / (1.0 + deltatransvec(it-1) * p.taylor.pers); 
		}

		if ( set_rb )
			trans_equm[it].pi = trans_equm[it].rnom - trans_equm[it].rb;
		else
			trans_equm[it].rb = trans_equm[it].rnom - trans_equm[it].pi;

		// Guess price level
		if ( it == 0 )
			trans_equm[it].pricelev = 1.0;
		else
			trans_equm[it].pricelev = trans_equm[it-1].pricelev / (1.0 - deltatransvec(it-1) * trans_equm[it-1].pi);

		// Guess price adj cost
		if ( solveFlexPriceTransitions )
			trans_equm[it].priceadjust = 0.0;
		else
			trans_equm[it].priceadjust = (p.priceadjcost / 2.0) * pow(trans_equm[it].pi, 2) * trans_equm[it].output;
	}

	// Guess inflation and output growth
	trans_equm[Ttrans-1].pidot = 0;
	trans_equm[Ttrans-1].logydot = 0;

	for (int it=Ttrans-2; it>0; --it) {
		trans_equm[it].pidot = (trans_equm[it].pi - trans_equm[it].pi) / deltatransvec(it);
		trans_equm[it].logydot = (log(trans_equm[it+1].output) - log(trans_equm[it].output)) / deltatransvec(it);
	}

	for (int it=0; it<Ttrans; ++it) {
		// Guess wholesale price
		if ( solveFlexPriceTransitions ) {
			trans_equm[it].price_W = 1.0 - 1.0 / trans_equm[it].elast;
			trans_equm[it].firmdiscount = -1e6;
		}
		else {
			// Phillips curve for marginal costs
			switch ( p.firm_discount_rate_type ) {
				case FirmDiscountRateType::rho:
					trans_equm[it].firmdiscount = p.rho;
					break;
				case FirmDiscountRateType::rb_iss:
					trans_equm[it].firmdiscount = initial_equm.rb;
					break;
				case FirmDiscountRateType::ra_iss:
					trans_equm[it].firmdiscount = initial_equm.ra;
					break;
				case FirmDiscountRateType::rb_trans:
					trans_equm[it].firmdiscount = trans_equm[it].rb;
					break;
				case FirmDiscountRateType::ra_trans:
					trans_equm[it].firmdiscount = trans_equm[it].ra;
					break;
			}

			// Guess wholesale price
			hank_float_type mkup = (trans_equm[it].elast - 1.0) / trans_equm[it].elast;
			trans_equm[it].price_W = (trans_equm[it].firmdiscount - trans_equm[it].logydot) * trans_equm[it].pi * p.priceadjcost / trans_equm[it].elast
				+ mkup - trans_equm[it].pidot * p.priceadjcost / trans_equm[it].elast;
			trans_equm[it].price_W = fmin(max_price_W, fmax(min_price_W, trans_equm[it].price_W));
		}
	}

}

void IRF::set_shock_paths() {
	VectorXr tfp_Y, mpshock, riskaver;
	hank_float_type initial_mpshock = 0;
	if ( shock.type == ShockType::tfp_Y )
		tfp_Y = get_AR1_path_logs(Ttrans, initial_equm.tfp_Y, shock.size, shock.pers, deltatransvec, nendtrans);
	else if ( shock.type == ShockType::monetary ) {
		mpshock = get_AR1_path_levels(Ttrans, initial_mpshock, shock.size, shock.pers, deltatransvec, nendtrans);
		mpshock(seq(1, Ttrans-1)) = VectorXr::Constant(Ttrans-1, 0.0);
	}
	else if ( shock.type == ShockType::riskaver ) {
		if ( permanentShock )
			riskaver = VectorXr::Constant(Ttrans, p.riskaver + shock.size);
		else
			riskaver = get_AR1_path_levels(Ttrans, p.riskaver, shock.size, shock.pers, deltatransvec, nendtrans);
	}

	for (int it=0; it<Ttrans; ++it) {
		trans_equm[it].tfp_Y = 0;
		trans_equm[it].mpshock = 0;
		trans_equm[it].riskaver = 0;

		if ( shock.type == ShockType::tfp_Y )
			trans_equm[it].tfp_Y = tfp_Y[it];
		else if ( shock.type == ShockType::monetary )
			trans_equm[it].mpshock = mpshock[it];
		else if ( shock.type == ShockType::riskaver )
			trans_equm[it].riskaver = riskaver[it];
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