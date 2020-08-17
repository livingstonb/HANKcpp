#include <bellman.h>

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
		update(ss);

		newV = reshape_array(V, {model.ntot, 1, 1});
		lVdiff = (boost2eigen(lastV) - boost2eigen(newV)).lpNorm<Eigen::Infinity>();
		lastV = reshape_array(V, {model.ntot, 1, 1});

		check_progress(lVdiff, dispfreq, ii, vtol);

		++ii;
	}
}

void HJB::update(const SteadyState& ss) {
	OptConsumptionArgs optc_args;
	double VaF, VaB, VbF, VbB, prof_keep;
	double prof_common = p.lumptransfer + p.profdistfracL * (1.0 - p.corptax) * ss.profit;

	if ( p.taxHHProfitIncome )
		prof_keep = 1.0 - p.labtax;
	else
		prof_keep = 1.0;

	double_vector profW = prof_keep * p.profdistfracW
			* (1.0 - p.corptax) * ss.profit * model.profsharegrid.array();

	double_vector adrift = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			for (int iy=0; iy<p.ny; ++iy) {
				compute_derivatives(VaF, VbF, VaB, VbB, ia, ib, iy);

				// if (p.prodispshock)

				optc_args.gbdrift = bdrift(ib) + prof_common + profW(iy);
				optc_args.gnetwage = netwagegrid(iy);


				optimalConsumption(VbB, cB, hB, sB, HcB)


				V[ia][ib][iy] = 0.9 * V[ia][ib][iy];

			}
		}
	}

}

void HJB::compute_derivatives(
	double& VaF, double& VbF, double& VaB, double& VbB,
	int ia, int ib, int iy) const {
	// Forward derivatives
	if (ia < p.na - 1) {
		VaF = (V[ia+1][ib][iy] - V[ia][ib][iy]) / model.dagrid(ia);
		VaF = std::max(VaF, dVamin);
	}

	if (ib < p.nb - 1) {
		VbF = (V[ia][ib+1][iy] - V[ia][ib][iy]) / model.dbgrid(ib);
		VbF = std::max(VbF, dVbmin);
	}

	// Backward derivatives
	if (ia > 0) {
		VaB = (V[ia][ib][iy] - V[ia-1][ib][iy]) / model.dagrid(ia-1);
		VaB = std::max(VaB, dVamin);
	}

	if (ib > 0) {
		VbB = (V[ia][ib][iy] - V[ia][ib-1][iy]) / model.dbgrid(ib-1);
		VbB = std::max(VbB, dVbmin);
	}
}

ConUpwind HJB::optimal_consumption(double gbdrift, double gnetwage) {
	ConUpwind output;
	output.s = 0.0;

	if (p.laborsupply == LaborType::none) {
		output.h = 1.0;
		output.c = gbdrift * lh * gnetwage;

		if (lc > 0.0)
			output.Hc = model.util(output.c);
		else
			output.Hc = 1.0e12;
	}
	else if (p.laborsupply == LaborType::sep) {
		double hmin = std::max(-gbdrift / gnetwage + 1.0e-5, 0.0);

		if (p.imposeMaxHours)
			double hmax = 1.0;
		else
			double hmax = 100.0;

	}
	else if (p.laborsupply == LaborType::ghh) {

	}

	return output;
}

boost_array_type<double, 3> make_value_guess(const Model& model, const SteadyState& ss) {
	const Parameters& p = model.p;

	boost_array_type<double, 3> V(model.dims);
	double lc, u, llabdisutil = 0.0;
	const double sep_constant = 1.0 / 3.0;
	double_array wageexpr;
	double_array bdriftnn;

	bdriftnn = model.get_rb_effective().array() * model.bgrid.array();
	bdriftnn = bdriftnn.max(0.0);

	switch (p.laborsupply) {
		case LaborType::sep:
			wageexpr = ss.netwagegrid.array() * sep_constant;
			llabdisutil = ss.chi * pow(sep_constant, 1.0 + 1.0 / p.frisch) / (1.0 + 1.0 / p.frisch);
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

void check_progress(double vdiff, int freq, int ii, double vtol) {
	if ( ii == 0 )
		return;
	else if ( (ii == 1) | (ii % freq == 0) ) {
		std::cout << "Iteration " << ii << ", diff = " << vdiff << '\n';
	}

	if ( vdiff <= vtol )
		std::cout << "Converged after " << ii << " iterations." << '\n';
}
