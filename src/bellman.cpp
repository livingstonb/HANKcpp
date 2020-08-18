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

	class Policies {
		public:
			Policies(const boost_array_shape<double, 3>& dims) : c(dims), h(dims), s(dims) {};

			boost_array_type<double, 3> c, h, s;

			void update(int ia, int ib, int iy, const ConUpwind& uwF, const ConUpwind& uwB, const ConUpwind& uw0) {
				bool not_backward = (!uwB.valid) | (uwF.Hc >= uwB.Hc);
				bool not_forward = (!uwF.valid) | (uwB.Hc >= uwF.Hc);
				bool forward_better_than_nothing = uwF.Hc >= uw0.Hc;
				bool backward_better_than_nothing = uwB.Hc >= uw0.Hc;

				ConUpwind uw_selected;
				if ( uwF.valid & not_backward & forward_better_than_nothing )
					uw_selected = uwF;
				else if ( uwB.valid & not_forward & backward_better_than_nothing )
					uw_selected = uwB;
				else
					uw_selected = uw0;

				c[ia][ib][iy] = uw_selected.c;
				h[ia][ib][iy] = uw_selected.h;
				s[ia][ib][iy] = uw_selected.s;
			}
	};
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
		update(ss);

		newV = reshape_array(V, {model.ntot, 1, 1});
		lVdiff = (boost2eigen(lastV) - boost2eigen(newV)).lpNorm<Eigen::Infinity>();
		lastV = reshape_array(V, {model.ntot, 1, 1});

		check_progress(lVdiff, dispfreq, ii, vtol);

		++ii;
	}
}



void HJB::update(const SteadyState& ss) {
	ValueFnDerivatives derivs;

	ConUpwind upwindB, upwindF, upwind0, upwindBad;
	upwindBad.s = 0.0;
	upwindBad.Hc = -1.0e12;

	Policies policies(model.dims);

	double gbdrift, gnetwage, idioscale;
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

	double_vector adrift = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();

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

				// Upwinding forward
				if ( ib < p.nb - 1 )
					upwindF = optimal_consumption(derivs.VbF, gbdrift, gnetwage, ss.chi, idioscale);
				else
					upwindF = upwindBad;
				upwindF.valid = ( upwindF.s > 0.0 );

				// Upwinding backward
				if ( ib > 0 )
					upwindB = optimal_consumption(derivs.VbB, gbdrift, gnetwage, ss.chi, idioscale);
				else
					upwindB = optimal_consumption(derivs.StationaryPtOrLimit, gbdrift, gnetwage, ss.chi, idioscale);
				upwindB.valid = ( upwindB.s < 0.0 );

				// No drift
				upwind0 = optimal_consumption(derivs.StationaryPtOrLimit, gbdrift, gnetwage, ss.chi, idioscale);

				policies.update(ia, ib, iy, upwindF, upwindB, upwind0);

				V[ia][ib][iy] = 0.9 * V[ia][ib][iy];
			}
		}
	}

}

ValueFnDerivatives HJB::compute_derivatives(int ia, int ib, int iy) const {
	ValueFnDerivatives d;

	// Forward derivatives
	if (ia < p.na - 1) {
		d.VaF = (V[ia+1][ib][iy] - V[ia][ib][iy]) / model.dagrid(ia);
		d.VaF = std::max(d.VaF, dVamin);
	}

	if (ib < p.nb - 1) {
		d.VbF = (V[ia][ib+1][iy] - V[ia][ib][iy]) / model.dbgrid(ib);
		d.VbF = std::max(d.VbF, dVbmin);
	}

	// Backward derivatives
	if (ia > 0) {
		d.VaB = (V[ia][ib][iy] - V[ia-1][ib][iy]) / model.dagrid(ia-1);
		d.VaB = std::max(d.VaB, dVamin);
	}

	if (ib > 0) {
		d.VbB = (V[ia][ib][iy] - V[ia][ib-1][iy]) / model.dbgrid(ib-1);
		d.VbB = std::max(d.VbB, dVbmin);
	}

	return d;
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

ConUpwind HJB::optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	ConUpwind upwind;

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
		upwind.h = std::min(upwind.h, 1.0);

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

ConUpwind HJB::optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const {
	ConUpwind upwind;

	upwind.h = model.labdisutil1inv(p.labwedge * netwage / idioscale, chi);
	if ( p.imposeMaxHours )
		upwind.h = std::max(upwind.h, 1.0);

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


