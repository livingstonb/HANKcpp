#include <impulse_responses.h>

#include <parameters.h>
#include <model.h>
#include <equilibrium.h>
#include <iostream>
#include <functional>
#include <hank_numerics.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <distribution_statistics.h>
#include <upwinding.h>
#include <assert.h>
#include <model_functions.h>

#include <math.h>

#include <cminpack_wrapper.h>

using namespace std::placeholders;

using SolverArgsIRF = HANK::UniquePtrContainer<const Parameters, const Model, const EquilibriumInitial, EquilibriumFinal, const IRF>;

namespace {
	std::vector<hank_float_type> construct_cum_delta_trans_vector(double deltatranstot, double deltatransmin, int Ttrans);

	std::vector<hank_float_type> get_AR1_path_logs(int T, double x0, double eps, double rho, const std::vector<hank_float_type>& deltas, int nback);

	std::vector<hank_float_type> get_AR1_path_levels(int T, double x0, double eps, double rho, const std::vector<hank_float_type>& deltas, int nback);

	hank_float_type get_firmdiscount(FirmDiscountRateType discount_type, const EquilibriumInitial& initial_equm, const EquilibriumTrans& trans_equm);

	void make_transition_guesses(const Parameters& p, const hank_float_type* x, IRF* irf);

	int final_steady_state_obj_fn(void* solver_args_voidptr, int /* n */, const real *x, real *fvec, int /* iflag */ );
}

void TransShock::setup() {
	if ( type == ShockType::none ) {
		std::cerr << "Shock type has not been set\n";
		throw 0;
	}

	if ( size == HANK::ValueNotSet ) {
		if ( type == ShockType::tfp_Y )
			size = -0.02;
		else if ( type == ShockType::monetary )
			size = 0.0025 / 4.0;
		else if ( type == ShockType::riskaver )
			size = 0.001;
	}

	if ( pers == HANK::ValueNotSet )
		pers = exp(-0.5);
}

IRF::IRF(const Parameters& p_, const Model& model_, const EquilibriumInitial& iss_) : p(p_), model(model_), initial_equm(iss_) {}

void IRF::setup()
{
	npricetrans = (p.nocc + 2) * Ttrans - 1;

	if ( (p.capadjcost > 0.0) | (p.invadjcost > 0.0) )
		npricetrans += Ttrans;

	if ( p.adjGovBudgetConstraint == GovBCAdjType::debt )
		npricetrans += 1;

	shock.setup();

	flextransition = p.solveFlexPriceTransition;
	stickytransition = p.solveStickyPriceTransition;

	cumdeltatrans = construct_cum_delta_trans_vector(deltatranstot, deltatransmin, Ttrans);
	deltatransvec.push_back(cumdeltatrans[0]);
	for (int it=1; it<Ttrans; ++it)
		deltatransvec.push_back(cumdeltatrans[it] - cumdeltatrans[it-1]);
}

void IRF::compute()
{
	// Final steady state
	if ( !permanentShock )
		final_equm_ptr.reset(new EquilibriumFinal(initial_equm));
	else {
		find_final_steady_state();
	}

	for (int it=0; it<Ttrans; ++it)
		trans_equm.push_back(EquilibriumBase(initial_equm));
	set_shock_paths();

	// Guess log deviations from steady state
	std::vector<hank_float_type> xguess(npricetrans);
	std::fill(xguess.begin(), xguess.end(), 0);
	if ( (shock.type == ShockType::tfp_Y)  & !permanentShock ) {
		// Guess output changes so capital constant
		for (int it=0; it<Ttrans; ++it) {
			xguess[it] = log(fmax(p.capadjcost, 1.0) * trans_equm[it].tfp_Y / initial_equm.tfp_Y);
		}
	}
	trans_equm.clear();

	std::cout << '\n';
	HANK::horzline();
	std::cout << "FINDING TRANSITION EQUILIBRIUM\n\n";

	if ( solver == SolverType::broyden ) {
		std::function<void(int, const hank_float_type*, hank_float_type*)>
			obj_fn = std::bind(&IRF::transition_fcn, this, _1, _2, _3);

		std::vector<hank_float_type> z(npricetrans);
		std::fill(z.begin(), z.end(), 0);

		std::cout << " - Computing objective function at initial guess\n";
		obj_fn(npricetrans, xguess.data(), z.data());

		double jacstepprice = 1.0e-6;
		hank_float_type fjac[npricetrans * npricetrans];

		std::cout << " - Computing jacobian\n";
		computingJacobian = true;
		HankNumerics::jacobian_square(obj_fn, npricetrans, xguess.data(), z.data(), fjac, jacstepprice);
		computingJacobian = false;

		std::cout << " - Beginning optimization with Broyden Backstep routine\n";
		HankNumerics::broyden_backstep(obj_fn, npricetrans, xguess.data(),
			z.data(), fjac, maxIterTrans, tolTransition);
	}
	else {
		std::cerr << "Must select Broyden solver\n";
		throw 0;
	}
}

void IRF::transition_fcn(int /* n */, const hank_float_type *x, hank_float_type *fvec)
{
	trans_equm.clear();
	for (int it=0; it<Ttrans; ++it) {
		trans_equm.push_back(EquilibriumBase(initial_equm));
		trans_equm[it].tdelta = deltatransvec[it];
	}

	set_shock_paths();

	make_transition_guesses(p, x, this);
	compute_remaining_variables();

	solve_trans_equilibrium(trans_equm, p, initial_equm, *final_equm_ptr, deltatransvec.data());

	// Solve policies backward
	HJB hjb(p, model, final_equm_ptr->V);
	for (int it=Ttrans-2; it>=0; --it) {
		hjb.update(trans_equm[it]);
		trans_equm[it].V = hjb.V;
		trans_equm[it].policies = *hjb.optimal_decisions;
	}

	// Solve distribution forward
	std::vector<DistributionStatistics> trans_stats;
	trans_stats.push_back(DistributionStatistics(p, model, initial_equm.policies, initial_equm.density));
	StationaryDist sdist(initial_equm.density);
	for (int it=0; it<Ttrans-1; ++it) {
		sdist.gtol = 1.0e-9;
		sdist.compute(p, model, trans_equm[it], trans_equm[it].policies);

		trans_stats.push_back(DistributionStatistics(p, model, trans_equm[it].policies, sdist.density));
	}

	// Set residuals
	int ix = 0;
	std::vector<std::string> equation_names, variable_names;
	std::vector<hank_float_type> variable_values;

	fvec[ix] = trans_equm[0].capital / initial_equm.capital - 1.0;
	equation_names.push_back("Capital market clearing");
	variable_names.push_back("capital");
	variable_values.push_back(trans_equm[0].capital);
	++ix;

	for (int it=1; it<Ttrans; ++it) {
		fvec[ix] = trans_stats[it].Ea / (trans_equm[it].valcapital + trans_equm[it].equity_A) - 1.0;
		equation_names.push_back("Illiquid asset market clearing (" + std::to_string(it) + ")");
		variable_names.push_back("E[a]");
		variable_values.push_back(trans_stats[it].Ea);
		++ix;

		fvec[ix] = trans_stats[it].Eb / trans_equm[it].bond - 1.0;
		equation_names.push_back("Liquid asset market clearing (" + std::to_string(it) + ")");
		variable_names.push_back("E[b]");
		variable_values.push_back(trans_stats[it].Eb);
		++ix;
	}

	for (int it=0; it<Ttrans; ++it) {
		for (int io=0; io<p.nocc; ++io) {
			fvec[ix] = trans_stats[it].Elabor_occ[io] * model.occdist[io] / trans_equm[it].labor_occ[io] - 1.0;
			equation_names.push_back("Labor market clearing (" + std::to_string(it) + "," + std::to_string(io) + ")");
			variable_names.push_back("E[labor_io]");
			variable_values.push_back(trans_stats[it].Elabor_occ[io]);
			++ix;
		}

		if ( (p.capadjcost > 0) | (p.invadjcost > 0) ) {
			fvec[ix] = trans_equm[it].inv_cap_ratio / (trans_equm[it].investment / trans_equm[it].capital) - 1.0;
			equation_names.push_back("Inv-capital market clearing (" + std::to_string(it) + ")");
			variable_names.push_back("I/K");
			variable_values.push_back(trans_equm[it].inv_cap_ratio);
			++ix;
		}
	}

	// HANK::OptimStatus optim_status(equation_names, variable_names, fvec, variable_values);
	// HANK::print(optim_status);
	if ( !computingJacobian )
		HANK::print(HANK::OptimNorm(HANK::norm(fvec, npricetrans)));
}

void IRF::compute_remaining_variables()
{
	// Guesses
	bool set_rb = (shock.type != ShockType::monetary) | flextransition;
	bool set_pi = (shock.type == ShockType::monetary) & (!flextransition);

	assert( !(set_pi & set_rb) );
	assert( set_rb | set_pi );

	for (int it=0; it<Ttrans; ++it) {
		// Inflation and nominal interest rate
		double ygap = log(trans_equm[it].output / initial_equm.output);
		if ( p.taylor.use_feedback_rule ) {
			if ( set_rb ) {
				trans_equm[it].rnom = (initial_equm.rnom - p.taylor.coeff_pi * trans_equm[it].rb
					+ p.taylor.coeff_y * ygap + trans_equm[it].mpshock) / (1.0 - p.taylor.coeff_pi);
			}
			else if ( set_pi ) {
				trans_equm[it].rnom += p.taylor.coeff_pi * trans_equm[it].rb
					+ p.taylor.coeff_y * ygap + trans_equm[it].mpshock;
			}
		}
		else {
			// Partial adjustment rule
			if ( it == 0 )
				trans_equm[it].rnom += trans_equm[it].mpshock;
			else {
				trans_equm[it].rnom = (trans_equm[it-1].rnom + deltatransvec[it-1] * p.taylor.pers
					* (initial_equm.rnom + p.taylor.coeff_pi * trans_equm[it].pi + p.taylor.coeff_y * ygap + trans_equm[it].mpshock))
					/ (1.0 + deltatransvec[it-1] * p.taylor.pers); 
			}
		}

		if ( set_rb )
			trans_equm[it].pi = trans_equm[it].rnom - trans_equm[it].rb;
		else
			trans_equm[it].rb = trans_equm[it].rnom - trans_equm[it].pi;

		// Price level
		if ( it == 0 )
			trans_equm[it].pricelev = 1.0;
		else
			trans_equm[it].pricelev = trans_equm[it-1].pricelev / (1.0 - deltatransvec[it-1] * trans_equm[it-1].pi);

		// Guess price adj cost
		if ( flextransition )
			trans_equm[it].priceadjust = 0.0;
		else
			trans_equm[it].priceadjust = ModelFunctions::priceadjcost(trans_equm[it].pi, trans_equm[it].output, p.priceadjcost);
	}

	// Inflation and output growth
	trans_equm[Ttrans-1].pidot = 0;
	trans_equm[Ttrans-1].logydot = 0;

	for (int it=Ttrans-2; it>=0; --it) {
		trans_equm[it].pidot = (trans_equm[it+1].pi - trans_equm[it].pi) / deltatransvec[it];
		trans_equm[it].logydot = (log(trans_equm[it+1].output) - log(trans_equm[it].output)) / deltatransvec[it];
	}

	for (int it=0; it<Ttrans; ++it) {
		// Guess wholesale price
		if ( flextransition ) {
			trans_equm[it].price_W = 1.0 - 1.0 / p.elast;
			trans_equm[it].firmdiscount = -1e6;
		}
		else {
			// Phillips curve for marginal costs: set initial_equm[it].firmdiscount
			trans_equm[it].firmdiscount = get_firmdiscount(p.firm_discount_rate_type, initial_equm, trans_equm[it]);

			// Guess wholesale price
			hank_float_type mkup = (p.elast - 1.0) / p.elast;
			trans_equm[it].price_W = (trans_equm[it].firmdiscount - trans_equm[it].logydot)
				* trans_equm[it].pi * p.priceadjcost / p.elast
				+ mkup - trans_equm[it].pidot * p.priceadjcost / p.elast;
			trans_equm[it].price_W = fmin(max_price_W, fmax(min_price_W, trans_equm[it].price_W));
		}
	}

}

void IRF::set_shock_paths()
{
	std::vector<hank_float_type> mpshock, tfp_Y, riskaver;
	for (int it=0; it<Ttrans; ++it) {
		mpshock.push_back(0);
		tfp_Y.push_back(initial_equm.tfp_Y);
		riskaver.push_back(p.riskaver);
	}

	hank_float_type initial_mpshock = 0;
	if ( shock.type == ShockType::tfp_Y )
		tfp_Y = get_AR1_path_logs(Ttrans, initial_equm.tfp_Y, shock.size, shock.pers, deltatransvec, nendtrans);
	else if ( shock.type == ShockType::monetary ) {
		mpshock = get_AR1_path_levels(Ttrans, initial_mpshock, shock.size, shock.pers, deltatransvec, nendtrans);
	}
	else if ( shock.type == ShockType::riskaver ) {
		if ( permanentShock )
			std::fill(riskaver.begin(), riskaver.end(), p.riskaver + shock.size);
		else
			riskaver = get_AR1_path_levels(Ttrans, p.riskaver, shock.size, shock.pers, deltatransvec, nendtrans);
	}

	for (int it=0; it<Ttrans; ++it) {
		trans_equm[it].tfp_Y = tfp_Y[it];
		trans_equm[it].mpshock = mpshock[it];
		trans_equm[it].riskaver = riskaver[it];
	}
}

void IRF::find_final_steady_state()
{
	int n = p.nocc + 2;
	std::vector<hank_float_type> x(n);

	x[0] = initial_equm.capital;

	for (int io=0; io<p.nocc; ++io)
		x[io + 1] = initial_equm.labor_occ[io];

	x[p.nocc + 1] = log(initial_equm.rb);

	SolverArgsIRF args;
	args.ptr1.reset(&p);
	args.ptr2.reset(&model);
	args.ptr3.reset(&initial_equm);
	args.ptr5.reset(this);

	cminpack_hybrd1_wrapper(final_steady_state_obj_fn, &args, n, x.data());

	final_equm_ptr = std::move(args.ptr4);
}

namespace {
	std::vector<hank_float_type> construct_cum_delta_trans_vector(double deltatranstot, double deltatransmin, int Ttrans)
	{
		double lb = (Ttrans + 1) * (deltatranstot - Ttrans * deltatransmin) / (Ttrans * (deltatranstot - deltatransmin));
		double la = (Ttrans + 1) * deltatransmin - deltatransmin * lb;
		std::vector<hank_float_type> cumdeltatrans(Ttrans);

		for (int it=0; it<Ttrans; ++it) {
			double lit = static_cast<double>(it + 1) / (Ttrans + 1);
			cumdeltatrans[it] = la * lit / (1.0 - lb * lit);
		}

		return cumdeltatrans;
	}

	std::vector<hank_float_type> get_AR1_path_logs(int T, double x0, double eps, double rho, const std::vector<hank_float_type>& deltas, int nback) {
		std::vector<hank_float_type> y(T);
		double rho_it;

		y[0] = x0 * exp(eps);
		for (int it=1; it<T-nback; ++it) {
			rho_it = pow(rho, deltas[it-1]);
			y[it] = pow(x0, 1.0 - rho_it) * pow(y[it-1], rho_it);
		}

		for (int it=T-nback; it<T; ++it)
			y[it] = x0;

		return y;
	}

	std::vector<hank_float_type> get_AR1_path_levels(int T, double x0, double eps, double rho, const std::vector<hank_float_type>& deltas, int nback) {
		std::vector<hank_float_type> y(T);
		double rho_it;

		y[0] = x0 + eps;
		for (int it=1; it<T-nback; ++it) {
			rho_it = pow(rho, deltas[it-1]);
			y[it] = x0 * (1.0 - rho_it) + y[it-1] * rho_it;
		}
		
		for (int it=T-nback; it<T; ++it)
			y[it] = x0;

		return y;
	}

	hank_float_type get_firmdiscount(FirmDiscountRateType discount_type, const EquilibriumInitial& initial_equm, const EquilibriumTrans& trans_equm_el)
	{
		switch ( discount_type ) {
			case FirmDiscountRateType::rho:
				return initial_equm.rho;
			case FirmDiscountRateType::rb_iss:
				return initial_equm.rb;
			case FirmDiscountRateType::ra_iss:
				return initial_equm.ra;
			case FirmDiscountRateType::rb_trans:
				return trans_equm_el.rb;
			case FirmDiscountRateType::ra_trans:
				return trans_equm_el.ra;
		}

		std::cerr << "Logic error\n";
		throw 0;
	}

	void make_transition_guesses(const Parameters& p, const hank_float_type* x, IRF* irf)
	{
		bool guess_rb = (irf->shock.type != ShockType::monetary) | irf->flextransition;
		bool guess_pi = (irf->shock.type == ShockType::monetary) & (!irf->flextransition);

		assert( !(guess_pi & guess_rb) );
		assert( guess_rb | guess_pi );

		const EquilibriumFinal& final_equm = *(irf->final_equm_ptr);
		auto& trans_equm = irf->trans_equm;
		int Ttrans = irf->Ttrans;
		int ix = 0;

		// Guess for output
		for (int it=0; it<Ttrans; ++it) {
			trans_equm[it].output *= exp(x[it]);
			++ix;
		}

		// Guess for rb or pi
		for (int it=0; it<Ttrans-1; ++it) {
			if ( guess_rb )
				trans_equm[it].rb += x[ix];
			else if ( guess_pi )
				trans_equm[it].pi += x[ix];

			++ix;
		}

		if ( guess_rb )
			trans_equm[Ttrans-1].rb = final_equm.rb;
		else if ( guess_pi )
			trans_equm[Ttrans-1].pi = final_equm.pi;

		// Guess for occupation-specific labor
		for (int it=0; it<Ttrans; ++it) {
			for (int io=0; io<final_equm.nocc; ++io) {
				trans_equm[it].labor_occ[io] *= exp(x[ix]);
				++ix;
			}
		}

		// Guess for capital
		for (int it=0; it<Ttrans; ++it) {
			if ( (p.capadjcost > 0) | (p.invadjcost > 0) ) {
				trans_equm[it].qcapital += exp(x[ix]);
				++ix;
			}
		}
	}

	int final_steady_state_obj_fn(void* solver_args_voidptr, int /* n */, const real *x, real *fvec, int /* iflag */ )
	{
		SolverArgsIRF& solver_args = *(SolverArgsIRF *) solver_args_voidptr;
		const Parameters& p = *(solver_args.ptr1);
		const Model& model = *(solver_args.ptr2);
		const EquilibriumInitial& iss = *(solver_args.ptr3);
		const IRF& irf = *(solver_args.ptr5);
		
		solver_args.ptr4.reset(new EquilibriumFinal(iss));
		EquilibriumFinal& final_ss = *solver_args.ptr4;

		if ( irf.permanentShock & (irf.shock.type == ShockType::riskaver) )
			final_ss.riskaver += irf.shock.size;

		final_ss.solve(p, iss, x);	

		HJB hjb(p, model, final_ss);
		hjb.iterate(final_ss);

		StationaryDist sdist;
		sdist.gtol = 1.0e-9;
		sdist.compute(p, model, final_ss, *hjb.optimal_decisions);

		DistributionStatistics stats(p, model, *hjb.optimal_decisions, sdist.density);

		std::vector<std::string> variable_names, equation_names;
		std::vector<hank_float_type> variable_values;

		// Function value error
		fvec[0] = stats.Ea / (final_ss.capital + final_ss.equity_A) - 1.0;
		equation_names.push_back("Illiquid asset market clearing");
		variable_values.push_back(stats.Ea);
		variable_names.push_back("E[a]");

		for (int io=0; io<p.nocc; ++io) {
			fvec[io+1] = stats.Elabor_occ[io] * model.occdist[io] / final_ss.labor_occ[io] - 1.0;
			equation_names.push_back("Labor market clearing (" + std::to_string(io) + ")");
			variable_values.push_back(stats.Elabor_occ[io] * model.occdist[io]);
			variable_names.push_back("labor_occ_" + std::to_string(io));
		}

		fvec[p.nocc+1] = stats.Eb / final_ss.equity_B - 1.0;
		equation_names.push_back("Liquid asset marking clearing");
		variable_values.push_back(stats.Eb);
		variable_names.push_back("E[b]");

		HANK::OptimStatus optim_status(equation_names, variable_names, fvec, variable_values);
		HANK::print(optim_status);

		return 0;
	}
}