#include <bellman.h>
#include <parameters.h>
#include <model.h>
#include <transition_matrix.h>
#include <upwinding.h>
#include <equilibrium.h>
#include <hank_numerics.h>
#include <utilities.h>
#include <adjustment_costs.h>

#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>

#include <functional>
#include <algorithm>
#include <math.h>
#include <string>

#include <hank_macros.h>

using namespace std::placeholders;

namespace {
	constexpr bool is_stationary_pt_or_limit(double Vb);

	void check_progress(double vdiff, int freq, int ii, double vtol);

	vector3dr make_value_guess(const Model& model, const Equilibrium& ss, double riskaver);

	Upwinding::DepositUpwind optimal_deposits(const Model& model, double Va, double Vb, double a, double illprice);
}

HJB::HJB(const Model& model_, const Equilibrium& ss) : model(model_), p(model_.p), V(p.na, p.nb, model.ny), optimal_decisions(model.dims) {
	riskaver = ss.riskaver;
	V = make_value_guess(model, ss, riskaver);
}

void HJB::iterate(const Equilibrium& ss) {
	int ii = 0;
	double lVdiff = 1.0;

	VectorXr lastV = as_eigen<VectorXr>(V);
	Upwinding::Policies policies(model.dims);
	while ( (ii < maxiter) & (lVdiff > vtol) ) {
		policies = update_policies(ss);
		update_value_fn(ss, policies);

		VectorXr newV = as_eigen<VectorXr>(V);
		lVdiff = (lastV - newV).lpNorm<Eigen::Infinity>();
		lastV = as_eigen<VectorXr>(V);

		check_progress(lVdiff, dispfreq, ii, vtol);

		++ii;
	}

	if ( ii == maxiter )
		std::cout << "HJB did not converge" << '\n';

	optimal_decisions = policies;

	if ( global_hank_options->print_diagnostics )
		print_variables();
}

Upwinding::Policies HJB::update_policies(const Equilibrium& ss) {
	ValueFnDerivatives derivs;
	double chi = p.chi;

	Upwinding::ConUpwind upwindB, upwindF, upwind0, upwindBad;
	upwindBad.s = 0.0;

	Upwinding::DepositUpwind depositFB, depositBF, depositBB, depositBad;

	Upwinding::Policies policies(model.dims);

	double gbdrift, gnetwage, idioscale, illiq, labdisutil;
	bool worth_adjusting;

	double prof_keep;
	if ( p.taxHHProfitIncome )
		prof_keep = 1.0 - p.labtax;
	else
		prof_keep = 1.0;

	ArrayXr proftot = prof_keep * p.profdistfracW
		* (1.0 - p.corptax) * ss.profit * as_eigen_map<const ArrayXr>(model.profsharegrid);
	proftot += ss.lumptransfer + p.profdistfracL * (1.0 - p.corptax) * ss.profit;

	Eigen::Map<const ArrayXr> bgridvec(model.bgrid.data(), model.bgrid.size());
	VectorXr bdrift = as_eigen<ArrayXr>(model.get_rb_effective(ss.rb, ss.rborr)) * bgridvec;

	std::function<Upwinding::ConUpwind(double, double, double, double)> opt_c;
	if ( p.endogLabor )
		opt_c = std::bind(&HJB::optimal_consumption_sep_labor, this, _1, _2, _3, chi, _4);
	else
		opt_c = std::bind(&HJB::optimal_consumption_no_laborsupply, this, _1, _2, _3);

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			for (int iy=0; iy<model.ny; ++iy) {
				derivs = compute_derivatives(ia, ib, iy);

				idioscale = 1; // model.yprodgrid(iy);

				gbdrift = bdrift(ib) + proftot(iy);
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
				illiq = model.agrid[ia] * ss.illprice;

				// Deposit decision: a forward, b backward
				if ( (ia < p.na - 1) & (ib > 0) ) {
					depositFB = optimal_deposits(model, derivs.VaF, derivs.VbB, illiq, ss.illprice);
					depositFB.valid = ( (depositFB.d > 0) & (depositFB.Hd > 0) );
				}
				else
					depositFB = depositBad;

				// Deposit decision: a backward, b forward
				if ( (ia > 0) & (ib < p.nb - 1) ) {
					depositBF = optimal_deposits(model, derivs.VaB, derivs.VbF, illiq, ss.illprice);
					worth_adjusting = ( depositBF.d <= -model.adjcosts->cost(depositBF.d, illiq) );
					depositBF.valid = ( worth_adjusting & (depositBF.Hd > 0) );
				}
				else
					depositBF = depositBad;

				// Deposit decision: a backward, b backward
				if ( ia > 0 ) {
					if ( ib == 0 )
						derivs.VbB = model.util1(upwindB.c, riskaver);

					depositBB = optimal_deposits(model, derivs.VaB, derivs.VbB, illiq, ss.illprice);
					worth_adjusting = ( depositBB.d > -model.adjcosts->cost(depositBB.d, illiq) );
					depositBB.valid = ( worth_adjusting & (depositBB.d <= 0) & (depositBB.Hd > 0));
				}
				else
					depositBB = depositBad;

				// Update d
				policies.update_d(ia, ib, iy, depositFB, depositBF, depositBB);

				// Update u
				labdisutil = idioscale * model.labdisutil(policies.h(ia,ib,iy), chi);
				if ( p.endogLabor )
					policies.u(ia,ib,iy) = model.util(policies.c(ia,ib,iy), riskaver) - labdisutil / p.labwedge;
				else
					policies.u(ia,ib,iy) = model.util(policies.c(ia,ib,iy), riskaver);
			}
		}
	}
	return policies;
}

ValueFnDerivatives HJB::compute_derivatives(int ia, int ib, int iy) const {
	ValueFnDerivatives d;

	// Forward derivatives
	if (ia < p.na - 1) {
		d.VaF = (V(ia+1,ib,iy) - V(ia,ib,iy)) / model.dagrid[ia];
		d.VaF = fmax(d.VaF, dVamin);
	}

	if (ib < p.nb - 1) {
		d.VbF = (V(ia,ib+1,iy) - V(ia,ib,iy)) / model.dbgrid[ib];
		d.VbF = fmax(d.VbF, dVbmin);
	}

	// Backward derivatives
	if (ia > 0) {
		d.VaB = (V(ia,ib,iy) - V(ia-1,ib,iy)) / model.dagrid[ia-1];
		d.VaB = fmax(d.VaB, dVamin);
	}

	if (ib > 0) {
		d.VbB = (V(ia,ib,iy) - V(ia,ib-1,iy)) / model.dbgrid[ib-1];
		d.VbB = fmax(d.VbB, dVbmin);
	}

	return d;
}

Upwinding::ConUpwind HJB::optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const {
	Upwinding::ConUpwind upwind;
	upwind.h = 1.0;

	if ( !is_stationary_pt_or_limit(Vb)) {
		upwind.c = model.util1inv(Vb, riskaver);
		upwind.s = bdrift + upwind.h * netwage - upwind.c;
	}
	else {
		upwind.s = 0.0;
		upwind.c = bdrift + upwind.h * netwage;
	}

	if ( upwind.c > 0.0 )
		upwind.Hc = model.util(upwind.c, riskaver) + Vb * upwind.s;
	else
		upwind.Hc = -1.0e12;

	return upwind;
}

Upwinding::ConUpwind HJB::optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	Upwinding::ConUpwind upwind;

	if ( !is_stationary_pt_or_limit(Vb) ) {
		upwind.h = model.labdisutil1inv(p.labwedge * netwage * Vb, chi);
	}
	else {
		double hmin = fmax(-bdrift / netwage + 1.0e-5, 0.0);
		double hmax = ( p.imposeMaxHours ) ? 1 : 100;
		double v1_at_min = model.util1BC(hmin, riskaver, chi, bdrift, netwage, idioscale);
		double v1_at_max = model.util1BC(hmax, riskaver, chi, bdrift, netwage, idioscale);

		if ( (v1_at_max > 0) & (v1_at_min < 0) ) {
			std::function<double(double)> objective = std::bind(
				&Model::util1BC, &model, _1, riskaver, chi, bdrift, netwage, idioscale);
			double facc = 1.0e-8;
			upwind.h = HankNumerics::rtsec(objective, hmin, hmax, facc);
		}
		else if ( v1_at_max <= 0 )
			upwind.h = hmax;
		else if ( v1_at_min >= 0 )
			upwind.h = hmin;
		else {
			std::cerr << "v1_at_min = " << v1_at_min << '\n';
			std::cerr << "v1_at_max = " << v1_at_max << '\n';
			std::cerr << "Logic error\n";
			throw 0;
		}
	}

	if ( p.imposeMaxHours )
		upwind.h = fmin(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( !is_stationary_pt_or_limit(Vb) ) {
		upwind.c = model.util1inv(Vb, riskaver);
		upwind.s = bdrift + upwind.h * netwage - upwind.c;
	}
	else {
		upwind.c = bdrift + upwind.h * netwage;
		upwind.s = 0.0;
	}

	if ( (upwind.c > 0.0) & (!is_stationary_pt_or_limit(Vb)) )
		upwind.Hc = model.util(upwind.c, riskaver) - labdisutil / p.labwedge + Vb * upwind.s;
	else if ( upwind.c > 0.0 )
		upwind.Hc = model.util(upwind.c, riskaver) - labdisutil / p.labwedge;
	else
		upwind.Hc = -1.0e12;

	return upwind;
}

void HJB::update_value_fn(const Equilibrium& ss, const Upwinding::Policies& policies) {
	VectorXr bvec(p.nb * p.na);
	VectorXr ycol, vcol(model.ny);
	int iab;
	bool kfe = false;

	SparseXd ldiagmat, sparseI = speye(p.na * p.nb);

	for ( int iy=0; iy<model.ny; ++iy ) {
		ycol = model.prodmarkovscale * model.matrices->ymarkovoff.row(iy);

		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				iab = TO_INDEX_1D(ia, ib, p.na, p.nb);

				// Vector of constants, bvec
				for (int iy2=0; iy2<model.ny; ++iy2)
					vcol[iy2] = V(ia,ib,iy2);

				bvec(iab) = delta * policies.u(ia,ib,iy) + V(ia,ib,iy) + delta * ycol.dot(vcol);
			}
		}

		SparseMatContainer Acont = construct_transition_matrix(p, model, ss.ra, ss.illprice, ss.illpricedot, policies, iy, kfe);
		SparseXd& A = Acont.get();

		// Construct B matrix = I + delta * (rho * I - A)
		double ldiag = 1.0 + delta * (p.rho + p.deathrate) - delta * model.prodmarkovscale * model.matrices->ymarkovdiag(iy,iy);
		ldiagmat = sparseI * ldiag;
		SparseXd B = ldiagmat - delta * A;
		B.makeCompressed();

		sparse_solver solver;
		solver.compute(B);
		if ( solver.info() != Eigen::Success ) {
			std::cerr << "Sparse solver failure" << '\n';
			throw 0;
		}

		Eigen::VectorXd bvecd = bvec.cast<double>();
		VectorXr v_vec = solver.solve(bvecd).cast<hank_float_type>();
		if ( solver.info() != Eigen::Success ) {
			std::cerr << "Sparse solver failure" << '\n';
			throw 0;
		}

		for (int ia=0; ia<p.na; ++ia)
			for (int ib=0; ib<p.nb; ++ib)
				V(ia,ib,iy) = v_vec[TO_INDEX_1D(ia, ib, na, p.nb)];
	}
}

void HJB::print_variables() const {
	std::cout << '\n';
	HankUtilities::horzline();
	std::cout << "SELECTED OUTPUT FROM BELLMAN:\n";

	std::vector<std::string> names;
	std::vector<double> values;

	names.push_back("c(1,1,1)");
	values.push_back(optimal_decisions.c(1,1,1));

	names.push_back("c(1,10,1)");
	values.push_back(optimal_decisions.c(1,10,1));

	names.push_back("c(10,1,1)");
	values.push_back(optimal_decisions.c(10,1,1));

	names.push_back("c(20,20,5)");
	values.push_back(optimal_decisions.c(20,20,5));

	names.push_back("d(1,1,1)");
	values.push_back(optimal_decisions.d(1,1,1));

	names.push_back("d(1,10,1)");
	values.push_back(optimal_decisions.d(1,10,1));

	names.push_back("d(10,1,1)");
	values.push_back(optimal_decisions.d(10,1,1));

	names.push_back("d(20,20,5)");
	values.push_back(optimal_decisions.d(20,20,5));

	HankUtilities::print_values(names, values);

	HankUtilities::horzline();
}

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

	vector3dr make_value_guess(const Model& model, const Equilibrium& ss, double riskaver) {
		const Parameters& p = model.p;

		vector3dr V(p.na, p.nb, model.ny);
		double lc, u;

		auto bgridvec = as_eigen_map<const ArrayXr>(model.bgrid);
		VectorXr bdriftnn = as_eigen<ArrayXr>(model.get_rb_effective(ss.rb, ss.rborr)) * bgridvec;

		auto agridvec = as_eigen_map<const ArrayXr>(model.agrid);
		ArrayXr adriftnn = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * agridvec;

		for (int ia=0; ia<p.na; ++ia) {
			for (int ib=0; ib<p.nb; ++ib) {
				for (int iy=0; iy<model.ny; ++iy) {
					lc = ss.netwagegrid[iy] + ss.lumptransfer + bdriftnn(ib) + adriftnn(ia);
					u = model.util(lc, riskaver);
					V(ia,ib,iy) = u / (p.rho + p.deathrate);
				}
			}
		}

		return V;
	}

	Upwinding::DepositUpwind optimal_deposits(const Model& model, double Va, double Vb, double a, double illprice) {
		Upwinding::DepositUpwind dupwind;
		dupwind.d = model.adjcosts->cost1inv(Va / (Vb * illprice) - 1.0, a);

		dupwind.Hd = Va * dupwind.d / illprice - Vb * (dupwind.d + model.adjcosts->cost(dupwind.d, a));
		return dupwind;
	}
}