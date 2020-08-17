#include <bellman.h>

HJB::HJB(const Model& model, const SteadyState& ss) : V(model.dims) {
	V = make_value_guess(model, ss);
}

void HJB::iterate(const Model& model, const SteadyState& ss) {
	const Parameters& p = model.p;

	double_vector adrift = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();


}

boost_array_type<double, 3> make_value_guess(const Model& model, const SteadyState& ss) {
	const Parameters& p = model.p;

	boost_array_type<double, 3> V(model.dims);
	double lc, u, llabdisutil = 0.0;
	double sep_constant = 1.0 / 3.0;
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