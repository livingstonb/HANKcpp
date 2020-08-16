#ifndef _BELLMAN_H
#define _BELLMAN_H

class HJB {
	public:
		HJB(const Model& model, const SteadyState& ss)

		void iterate();
};

void iterate(const Model& model, const SteadyState& ss) {
	const Parameters& p = model.p;

	double_vector adrift = (ss.ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();
	double_vector bdrift = model.get_rb_effective().array() * model.bgrid.array();


}


#endif