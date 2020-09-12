#include <bellman.h>
#include <parameters.h>
#include <model.h>
#include <upwinding.h>
#include <steady_state.h>
#include <hank_numerics.h>
#include <utilities.h>

#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>
#include <hank_boost_eigen_routines.h>
#include <transition_matrix.h>

#include <algorithm>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

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

	StdVector3d<double> make_value_guess(const Model& model, const SteadyState& ss) {
		const Parameters& p = model.p;

		StdVector3d<double> V(p.na, p.nb, p.ny);
		double lc, u, llabdisutil = 0.0;
		const double sep_constant = 1.0 / 3.0;
		double_array wageexpr, bdriftnn;

		bdriftnn = model.get_rb_effective().array() * model.bgrid.array();

		switch ( p.laborsupply ) {
			case LaborType::sep:
				wageexpr = vector2eigenv(ss.netwagegrid).array() * sep_constant;
				llabdisutil = model.labdisutil(sep_constant, ss.chi);
				break;
			case LaborType::none:
				wageexpr = vector2eigenv(ss.netwagegrid).array();
				break;
			case LaborType::ghh:
				wageexpr = vector2eigenv(ss.netwagegrid).array().pow(1 + p.frisch);
				break;
		}

		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				for (int iy=0; iy<p.ny; ++iy) {
					lc = wageexpr(iy) + p.lumptransfer + bdriftnn(ib);
					u = model.util(lc);
					V(ia,ib,iy) = (u - llabdisutil) / (p.rho + p.deathrate);
				}
			}
		}

		return V;
	}

	Upwinding::DepositUpwind optimal_deposits(const Model& model, double Va, double Vb, double a) {
		Upwinding::DepositUpwind dupwind;

		dupwind.d = model.adjcosts.cost1inv(Va / Vb - 1.0, a);
		dupwind.d = fmin(dupwind.d, model.p.dmax);

		dupwind.Hd = Va * dupwind.d - Vb * (dupwind.d + model.adjcosts.cost(dupwind.d, a));
		return dupwind;
	}
}

HJB::HJB(const Model& model_, const SteadyState& ss) : model(model_), p(model_.p), V(p.na, p.nb, p.ny), optimal_decisions(model.dims) {
	V = make_value_guess(model, ss);
}

void HJB::iterate(const SteadyState& ss) {
	int ii = 0;
	double lVdiff = 1.0;

	double_vector lastV = to_eigenv(V);
	Upwinding::Policies policies(model.dims);
	while ( (ii < maxiter) & (lVdiff > vtol) ) {
		policies = update_policies(ss);
		update_value_fn(ss, policies);

		double_vector newV = to_eigenv(V);
		lVdiff = (lastV - newV).lpNorm<Eigen::Infinity>();
		lastV = to_eigenv(V);

		check_progress(lVdiff, dispfreq, ii, vtol);

		++ii;
	}

	if ( ii == maxiter )
		std::cout << "HJB did not converge" << '\n';

	optimal_decisions = policies;
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

	bool labor_is_ghh = (p.laborsupply == LaborType::ghh);
	const bool scale_wrt_ss = (!p.scaleDisutilityIdio) & p.prodispshock
		& p.prodDispScaleDisutility & labor_is_ghh;

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
					idioscale = model.yprodgrid(iy) / ss.yprodgrid[iy];
				else
					idioscale = 1.0;

				gbdrift = bdrift(ib) + prof_common + profW(iy);
				gnetwage = ss.netwagegrid[iy];

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
					if ( ib == 0 ) {
						labdisutil = idioscale * model.labdisutil(upwindB.h, chi);

						if ( p.laborsupply == LaborType::ghh )
							derivs.VbB = model.util1(upwindB.c - labdisutil / p.labwedge);
						else
							derivs.VbB = model.util1(upwindB.c);
					}
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
			}
		}
	}
	return policies;
}

ValueFnDerivatives HJB::compute_derivatives(int ia, int ib, int iy) const {
	ValueFnDerivatives d;

	// Forward derivatives
	if (ia < p.na - 1) {
		d.VaF = (V(ia+1,ib,iy) - V(ia,ib,iy)) / model.dagrid(ia);
		d.VaF = fmax(d.VaF, dVamin);
	}

	if (ib < p.nb - 1) {
		d.VbF = (V(ia,ib+1,iy) - V(ia,ib,iy)) / model.dbgrid(ib);
		d.VbF = fmax(d.VbF, dVbmin);
	}

	// Backward derivatives
	if (ia > 0) {
		d.VaB = (V(ia,ib,iy) - V(ia-1,ib,iy)) / model.dagrid(ia-1);
		d.VaB = fmax(d.VaB, dVamin);
	}

	if (ib > 0) {
		d.VbB = (V(ia,ib,iy) - V(ia,ib-1,iy)) / model.dbgrid(ib-1);
		d.VbB = fmax(d.VbB, dVbmin);
	}

	return d;
}

Upwinding::ConUpwind HJB::optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const {
	Upwinding::ConUpwind upwind;
	upwind.h = 1.0;

	if ( !is_stationary_pt_or_limit(Vb)) {
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

	if ( !is_stationary_pt_or_limit(Vb) ) {
		upwind.h = model.labdisutil1inv(p.labwedge * netwage * Vb / idioscale, chi);
	}
	else {
		upwind.s = 0.0;
		const double hmin = fmax(-bdrift / netwage + 1.0e-5, 0.0);
		const double hmax = ( p.imposeMaxHours ) ? 1 : 100;
		const double v1_at_min = model.util1BC(hmin, chi, bdrift, netwage, idioscale);
		const double v1_at_max = model.util1BC(hmax, chi, bdrift, netwage, idioscale);

		if ( (v1_at_max > 0) & (v1_at_min < 0) ) {
			std::function<double(double)> objective = [=] (double h) {
				return model.util1BC(h, chi, bdrift, netwage, idioscale);
			};
			const double facc = 1.0e-8;
			upwind.h = HankNumerics::rtsec(objective, hmin, hmax, facc);
		}
		else if ( v1_at_max <= 0.0 )
			upwind.h = hmax;
		else if ( v1_at_min >= 0.0 )
			upwind.h = hmin;
		else
			throw "Logic error";
	}

	if ( p.imposeMaxHours )
		upwind.h = fmin(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( !is_stationary_pt_or_limit(Vb) ) {
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

	if ( !is_stationary_pt_or_limit(Vb) ) {
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
	double_vector bvec(p.nb * p.na);
	double_vector ycol, vcol(p.ny);
	boost_index indices;
	int iab;
	Bellman::Drifts drifts;
	bool kfe = false;

	sparse_matrix ldiagmat, sparseI = speye(p.na * p.nb);

	double_vector adriftvec = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdriftvec = model.get_rb_effective().array() * model.bgrid.array();

	int na = p.na;
	for ( int iy=0; iy<p.ny; ++iy ) {
		triplet_list Aentries;
		ycol = model.prodmarkovscale * model.ymarkovoff.row(iy);

		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				iab = TO_INDEX_1D(ia, ib, p.na);

				// Vector of constants, bvec
				for (int iy2=0; iy2<p.ny; ++iy2)
					vcol[iy2] = V(ia,ib,iy2);

				bvec(iab) = delta * policies.u[ia][ib][iy] + V(ia,ib,iy) + delta * ycol.dot(vcol);
			}
		}

		// sparse_matrix A = construct_A_matrix(ss, policies, iy, kfe);
		sparse_matrix A = construct_transition_matrix(p, model, ss.ra, policies, iy, kfe);

		// Construct B matrix = I + delta * (rho * I - A)
		double ldiag = 1.0 + delta * (p.rho + p.deathrate) - delta * model.prodmarkovscale * model.ymarkovdiag(iy,iy);
		ldiagmat = sparseI * ldiag;
		sparse_matrix B = ldiagmat - delta * A;
		B.makeCompressed();

		sparse_solver solver;
		solver.compute(B);
		if ( solver.info() != Eigen::Success )
			throw "Sparse solver failure";

		double_vector v_vec = solver.solve(bvec);
		if ( solver.info() != Eigen::Success )
			throw "Sparse solver failure";

		for (int ia=0; ia<p.na; ++ia)
			for (int ib=0; ib<p.nb; ++ib)
				V(ia,ib,iy) = v_vec[TO_INDEX_1D(ia, ib, na)];
	}
}