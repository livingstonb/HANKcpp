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
	ConUpwind upwindB, upwindF;
	double VaF, VaB, VbF, VbB, prof_keep;
	double gbdrift, gnetwage, idioscale, chi = ss.chi;
	const double prof_common = p.lumptransfer + p.profdistfracL * (1.0 - p.corptax) * ss.profit;

	bool labor_is_separable = (p.laborsupply == LaborType::sep);
	const bool scale_wrt_ss = (!p.scaleDisutilityIdio) & p.prodispshock
		& p.prodDispScaleDisutility & labor_is_separable;

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

				if ( p.scaleDisutilityIdio )
					idioscale = model.yprodgrid(iy);
				else if ( scale_wrt_ss )
					idioscale = model.yprodgrid(iy) / ss.yprodgrid(iy);
				else
					idioscale = 1.0;

				gbdrift = bdrift(ib) + prof_common + profW(iy);
				gnetwage = ss.netwagegrid(iy);

				upwindB = optimal_consumption(VbB, gbdrift, gnetwage, chi, idioscale);
				upwindF = optimal_consumption(VbF, gbdrift, gnetwage, chi, idioscale);


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

ConUpwind HJB::optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	switch ( p.laborsupply ) {
		case LaborType::none:
			return optimal_consumption_no_laborsupply(Vb, bdrift, netwage);
			break;
		case LaborType::sep:
			return optimal_consumption_sep_labor(Vb, bdrift, netwage, chi, idioscale);
			break;
		case LaborType::ghh:
			return optimal_consumption_ghh_labor(Vb, bdrift, netwage, chi, idioscale);
			break;
	}

	throw "Logic error";
}

ConUpwind HJB::optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const {
	ConUpwind upwind;
	upwind.h = 1.0;

	if ( Vb > -999.0 ) {
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

ConUpwind HJB::optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	ConUpwind upwind;

	if ( Vb > -999.0 ) {
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
		upwind.h = std::min(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( Vb > -999.0 ) {
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

ConUpwind HJB::optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	ConUpwind upwind;

	upwind.h = model.labdisutil1inv(p.labwedge * netwage / idioscale, chi);
	if ( p.imposeMaxHours )
		upwind.h = std::max(upwind.h, 1.0);

	double labdisutil = idioscale * model.labdisutil(upwind.h, chi);

	if ( Vb > -999.0 ) {
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

void check_progress(double vdiff, int freq, int ii, double vtol) {
	if ( ii == 0 )
		return;
	else if ( (ii == 1) | (ii % freq == 0) ) {
		std::cout << "Iteration " << ii << ", diff = " << vdiff << '\n';
	}

	if ( vdiff <= vtol )
		std::cout << "Converged after " << ii << " iterations." << '\n';
}
