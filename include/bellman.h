#ifndef _BELLMAN_H
#define _BELLMAN_H

class HJB {
	public:
		HJB(const Model& model, const SteadyState& ss);

		void iterate();

		array_type<double, 3> V;
};

array_type<double, 3> make_value_guess();

void iterate(const Model& model, const SteadyState& ss) {
	const Parameters& p = model.p;

	double_vector adrift = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();


}


#endif