#include <bellman.h>

namespace {
	constexpr bool is_stationary_pt_or_limit(double Vb) {
		return (Vb <= ValueFnDerivatives::StationaryPtOrLimit);
	}

	void check_progress(double vdiff, int freq, int ii, double vtol) {
	if ( ii == 0 )
		return;
	else if ( (ii == 1) | (ii % freq == 0) ) {
		std::cout << "Iteration " << ii << ", diff = " << vdiff << '\n';
	}

	if ( vdiff <= vtol )
		std::cout << "Converged after " << ii << " iterations." << '\n';
	}

	boost3d make_value_guess(const Model& model, const SteadyState& ss) {
		const Parameters& p = model.p;

		boost3d V(model.dims);
		double lc, u, llabdisutil = 0.0;
		const double sep_constant = 1.0 / 3.0;
		double_array wageexpr, bdriftnn;

		bdriftnn = model.get_rb_effective().array() * model.bgrid.array();
		bdriftnn = bdriftnn.max(0.0);

		switch (p.laborsupply) {
			case LaborType::sep:
				wageexpr = ss.netwagegrid.array() * sep_constant;
				llabdisutil = model.labdisutil(sep_constant, ss.chi);
				break;
			case LaborType::none:
				wageexpr = ss.netwagegrid.array();
				break;
			case LaborType::ghh:
				wageexpr = ss.netwagegrid.array().pow(1 + p.frisch);
				break;
		}

		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				for (int iy=0; iy<p.ny; ++iy) {
					lc = wageexpr(iy) + p.lumptransfer + bdriftnn(ib);
					u = model.util(lc);
					V[ia][ib][iy] = (u - llabdisutil) / (p.rho + p.deathrate);
				}
			}
		}

		return V;
	}

	class Drifts {
		public:
			Drifts() {}
			Drifts(double s, double d, double areturn, double acost, bool kfe) {
				if ( kfe ) {
					aB = fmin(d + areturn, 0.0);
					aF = fmax(d + areturn, 0.0);
					bB = fmin(s - d - acost, 0.0);
					bF = fmax(s - d - acost, 0.0);
				}
				else {
					aB = fmin(d, 0.0) + fmin(areturn, 0.0);
					aF = fmax(d, 0.0) + fmax(areturn, 0.0);
					bB = fmin(-d - acost, 0) + fmin(s, 0.0);
					bF = fmax(-d - acost, 0) + fmax(s, 0.0);
				}
			}
			double aB, aF, bB, bF;
	};

	Upwinding::DepositUpwind optimal_deposits(const Model& model, double Va, double Vb, double a) {
		Upwinding::DepositUpwind dupwind;

		dupwind.d = model.adjcosts.cost1inv(Va / Vb - 1.0, a);
		dupwind.d = fmin(dupwind.d, model.p.dmax);

		dupwind.Hd = Va * dupwind.d - Vb * (dupwind.d + model.adjcosts.cost(dupwind.d, a));
		return dupwind;
	}
}

HJB::HJB(const Model& model_, const SteadyState& ss) : model(model_), p(model_.p), V(model.dims) {
	V = make_value_guess(model, ss);
}

void HJB::iterate(const SteadyState& ss) {
	int ii = 0;
	double lVdiff = 1.0;

	boost_array_shape<double, 3> flat_dims = {{model.ntot, 1, 1}};

	boost_array_type<double, 3> lastV(flat_dims);
	boost_array_type<double, 3> newV(flat_dims);
	lastV = reshape_array(V, {model.ntot, 1, 1});

	while ( (ii < maxiter) & (lVdiff > vtol) ) {
		auto policies = update_policies(ss);
		update_value_fn(ss, policies);

		newV = reshape_array(V, {model.ntot, 1, 1});
		lVdiff = (boost2eigen(lastV) - boost2eigen(newV)).lpNorm<Eigen::Infinity>();
		lastV = reshape_array(V, {model.ntot, 1, 1});

		check_progress(lVdiff, dispfreq, ii, vtol);

		++ii;
	}
}

Upwinding::Policies HJB::update_policies(const SteadyState& ss) {
	ValueFnDerivatives derivs;
	double chi = ss.chi;

	Upwinding::ConUpwind upwindB, upwindF, upwind0, upwindBad;
	upwindBad.s = 0.0;
	upwindBad.Hc = -1.0e12;

	Upwinding::DepositUpwind depositFB, depositBF, depositBB, depositBad;
	depositBad.Hd = -1.0e12;

	Upwinding::Policies policies(model.dims);

	double gbdrift, gnetwage, idioscale, illiq, labdisutil;
	bool worth_adjusting;
	const double prof_common = p.lumptransfer + p.profdistfracL * (1.0 - p.corptax) * ss.profit;

	bool labor_is_separable = (p.laborsupply == LaborType::sep);
	const bool scale_wrt_ss = (!p.scaleDisutilityIdio) & p.prodispshock
		& p.prodDispScaleDisutility & labor_is_separable;

	double prof_keep;
	if ( p.taxHHProfitIncome )
		prof_keep = 1.0 - p.labtax;
	else
		prof_keep = 1.0;

	double_vector profW = prof_keep * p.profdistfracW
			* (1.0 - p.corptax) * ss.profit * model.profsharegrid.array();

	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();

	std::function<Upwinding::ConUpwind(double, double, double, double)> opt_c;
	if ( p.laborsupply == LaborType::none ) {
		opt_c = [this] (double Vb, double bdrift, double netwage, double) {
				return optimal_consumption_no_laborsupply(Vb, bdrift, netwage);
			};
		}
	else if ( p.laborsupply == LaborType::sep ) {
		opt_c = [this, chi] (double Vb, double bdrift, double netwage, double idioscale) {
				return optimal_consumption_sep_labor(Vb, bdrift, netwage, chi, idioscale);
			};
		}
	else if ( p.laborsupply == LaborType::ghh ) {
		opt_c = [this, chi] (double Vb, double bdrift, double netwage, double idioscale) {
				return optimal_consumption_ghh_labor(Vb, bdrift, netwage, chi, idioscale);
			};
		}
	else
		throw "Logic error";

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			for (int iy=0; iy<p.ny; ++iy) {
				derivs = compute_derivatives(ia, ib, iy);

				if ( p.scaleDisutilityIdio )
					idioscale = model.yprodgrid(iy);
				else if ( scale_wrt_ss )
					idioscale = model.yprodgrid(iy) / ss.yprodgrid(iy);
				else
					idioscale = 1.0;

				gbdrift = bdrift(ib) + prof_common + profW(iy);
				gnetwage = ss.netwagegrid(iy);

				// CONSUMPTION UPWINDING

				// Upwinding forward
				if ( ib < p.nb - 1 ) {
					upwindF = opt_c(derivs.VbF, gbdrift, gnetwage, idioscale);
					upwindF.valid = ( upwindF.s > 0.0 );
				}
				else
					upwindF = upwindBad;

				// Upwinding backward
				if ( ib > 0 )
					upwindB = opt_c(derivs.VbB, gbdrift, gnetwage, idioscale);
				else
					upwindB = opt_c(derivs.StationaryPtOrLimit, gbdrift, gnetwage, idioscale);
				upwindB.valid = ( upwindB.s < 0.0 );

				// No drift
				upwind0 = opt_c(derivs.StationaryPtOrLimit, gbdrift, gnetwage, idioscale);

				// Update c, s, and h
				policies.update_c(ia, ib, iy, upwindF, upwindB, upwind0);

				// DEPOSIT UPWINDING
				illiq = model.agrid(ia);

				// Deposit decision: a forward, b backward
				if ( (ia < p.na - 1) & (ib > 0) ) {
					depositFB = optimal_deposits(model, derivs.VaF, derivs.VbB, illiq);
					depositFB.valid = ( (depositFB.d > 0) & (depositFB.Hd > 0) );
				}
				else
					depositFB = depositBad;

				// Deposit decision: a backward, b forward
				if ( (ia > 0) & (ib < p.nb - 1) ) {
					depositBF = optimal_deposits(model, derivs.VaB, derivs.VbF, illiq);
					worth_adjusting = ( depositBF.d <= -model.adjcosts.cost(depositBF.d, illiq) );
					depositBF.valid = ( worth_adjusting & (depositBF.Hd > 0) );
				}
				else
					depositBF = depositBad;

				// Deposit decision: a backward, b backward
				if ( ia > 0 ) {
					labdisutil = idioscale * model.labdisutil(upwindB.h, chi);

					if ( p.laborsupply == LaborType::ghh )
						derivs.VbB = model.util1(upwindB.c - labdisutil / p.labwedge);
					else
						derivs.VbB = model.util1(upwindB.c);

					depositBB = optimal_deposits(model, derivs.VaB, derivs.VbB, illiq);
					worth_adjusting = ( depositBB.d > -model.adjcosts.cost(depositBB.d, illiq) );
					depositBB.valid = ( worth_adjusting & (depositBB.d <= 0) & (depositBB.Hd > 0));
				}
				else
					depositBB = depositBad;

				// Update d
				policies.update_d(ia, ib, iy, depositFB, depositBF, depositBB);

				// Update u
				labdisutil = idioscale * model.labdisutil(policies.h[ia][ib][iy], chi);
				if ( p.laborsupply == LaborType::none )
					policies.u[ia][ib][iy] = model.util(policies.c[ia][ib][iy]);
				if ( p.laborsupply == LaborType::sep )
					policies.u[ia][ib][iy] = model.util(policies.c[ia][ib][iy]) - labdisutil / p.labwedge;
				else
					policies.u[ia][ib][iy] = model.util(policies.c[ia][ib][iy] - labdisutil / p.labwedge);

				V[ia][ib][iy] = 0.9 * V[ia][ib][iy];
			}
		}
	}
	return policies;
}

ValueFnDerivatives HJB::compute_derivatives(int ia, int ib, int iy) const {
	ValueFnDerivatives d;

	// Forward derivatives
	if (ia < p.na - 1) {
		d.VaF = (V[ia+1][ib][iy] - V[ia][ib][iy]) / model.dagrid(ia);
		d.VaF = fmax(d.VaF, dVamin);
	}

	if (ib < p.nb - 1) {
		d.VbF = (V[ia][ib+1][iy] - V[ia][ib][iy]) / model.dbgrid(ib);
		d.VbF = fmax(d.VbF, dVbmin);
	}

	// Backward derivatives
	if (ia > 0) {
		d.VaB = (V[ia][ib][iy] - V[ia-1][ib][iy]) / model.dagrid(ia-1);
		d.VaB = fmax(d.VaB, dVamin);
	}

	if (ib > 0) {
		d.VbB = (V[ia][ib][iy] - V[ia][ib-1][iy]) / model.dbgrid(ib-1);
		d.VbB = fmax(d.VbB, dVbmin);
	}

	return d;
}

Upwinding::ConUpwind HJB::optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const {
	Upwinding::ConUpwind upwind;
	upwind.h = 1.0;

	if ( is_stationary_pt_or_limit(Vb)) {
		upwind.c = model.util1inv(Vb);
		upwind.s = bdrift + upwind.h * netwage - upwind.c;
	}
	else {
		upwind.s = 0.0;
		upwind.c = bdrift + upwind.h * netwage;
	}

	if ( upwind.c > 0.0 )
		upwind.Hc = model.util(upwind.c) + Vb * upwind.s;
	else
		upwind.Hc = -1.0e12;

	return upwind;
}

Upwinding::ConUpwind HJB::optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	Upwinding::ConUpwind upwind;

	if ( is_stationary_pt_or_limit(Vb) ) {
		upwind.h = model.labdisutil1inv(p.labwedge * netwage * Vb / idioscale, chi);
	}
	else {
		upwind.s = 0.0;
		const double hmin = std::max(-bdrift / netwage + 1.0e-5, 0.0);
		const double wscale = 1.0;
		const double hmax = ( p.imposeMaxHours ) ? 1 : 100;
		const double v1_at_min = model.util1BC(hmin, chi, bdrift, netwage, wscale);
		const double v1_at_max = model.util1BC(hmax, chi, bdrift, netwage, wscale);

		if ( v1_at_max <= 0.0 )
			upwind.h = hmax;
		else if ( v1_at_min >= 0.0 )
			upwind.h = hmin;
		else {
			std::function<double(double)> objective = [=] (double h) {
				return model.util1BC(h, chi, bdrift, netwage, wscale);
			};
			const double facc = 1.0e-8;
			upwind.h = HankNumerics::rtsec(objective, hmin, hmax, facc);
		}
	}

	if ( p.imposeMaxHours )
		upwind.h = fmin(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( is_stationary_pt_or_limit(Vb) ) {
		upwind.c = model.util1inv(Vb);
		upwind.s = bdrift + upwind.h * netwage - upwind.c;
	}
	else {
		upwind.c = bdrift + upwind.h * netwage;
	}

	if ( upwind.c > 0.0 )
		upwind.Hc = model.util(upwind.c) - labdisutil / p.labwedge + Vb * upwind.s;
	else
		upwind.Hc = -1.0e12;

	return upwind;
}

Upwinding::ConUpwind HJB::optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	Upwinding::ConUpwind upwind;

	upwind.h = model.labdisutil1inv(p.labwedge * netwage / idioscale, chi);
	if ( p.imposeMaxHours )
		upwind.h = fmax(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( is_stationary_pt_or_limit(Vb) ) {
		upwind.c = model.util1inv(Vb) + labdisutil / p.labwedge;
		upwind.s = bdrift + upwind.h * netwage - upwind.c;
	}
	else {
		upwind.s = 0;
		upwind.c = bdrift + upwind.h * netwage;
	}

	if ( upwind.c - labdisutil / p.labwedge > 0.0 )
		upwind.Hc = model.util(upwind.c - labdisutil / p.labwedge) + Vb * upwind.s;
	else
		upwind.Hc = -1.0e12;

	return upwind;
}

void HJB::update_value_fn(const SteadyState& ss, const Upwinding::Policies& policies) {
	double_vector bvec(model.ntot);
	double_vector ycol, vcol;
	boost3d::index_gen indices;
	boost1d vcol_boost = new_array<double, 1>({p.ny});
	double d, s, acost, areturn, val, val1, val2;
	int iab;
	Drifts drifts;
	bool kfe = false;

	double_vector adriftvec = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdriftvec = model.get_rb_effective().array() * model.bgrid.array();

	int na = p.na;
	auto to_ab_index = [na](int iia, int iib) { return iia + na * iib; };

	for ( int iy=0; iy<p.ny; ++iy ) {
		triplet_list Aentries;
		Aentries.reserve(5 * p.na * p.nb);
		ycol = model.prodmarkovscale * model.ymarkovoff.row(iy);

		iab = 0;
		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				d = policies.d[ia][ib][iy];
				s = policies.d[ia][ib][iy];
				acost = model.adjcosts.cost(d, model.agrid[ia]);
				areturn = adriftvec(ia);

				// Vector of constants
				vcol_boost = V[indices[ia][ib][range()]];
				vcol = boost2eigen(vcol_boost);
				bvec(iab) = delta * policies.u[ia][ib][iy] + V[ia][ib][iy] + delta * ycol.dot(vcol);

				// Compute drifts
				drifts = Drifts(s, d, areturn, acost, kfe);

				// Matrix entries
				if ( ia > 0 )
					val = -drifts.aB / model.dagrid(ia-1);
				else
					val = 0.0;

				if ( val != 0.0 )
					Aentries.push_back(triplet_type(iab, to_ab_index(ia-1, ib), val));

				if ( ib > 0 )
					val = -drifts.bB / model.dbgrid(ib-1);
				else
					val = 0.0;

				if ( val != 0.0 )
					Aentries.push_back(triplet_type(iab, to_ab_index(ia, ib-1), val));

				// Matrix entries -- diagonal
				if ( ia == 0 )
					val1 = -drifts.aF / model.dagrid(ia);
				else if ( ia == p.na - 1 )
					val1 = drifts.aB / model.dagrid(ia-1);
				else
					val1 = drifts.aB / model.dagrid(ia-1) - drifts.aF / model.dagrid(ia);

				if ( ib == 0 )
					val2 = -drifts.bF / model.dbgrid(ib);
				else if ( ib == p.nb - 1 )
					val2 = drifts.bB / model.dbgrid(ib-1);
				else
					val2 = drifts.bB / model.dbgrid(ib-1) - drifts.bF / model.dbgrid(ib);

				Aentries.push_back(triplet_type(iab, iab, val1 + val2));

				++iab;
			}
		}
	}
}