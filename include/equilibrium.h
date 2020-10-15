#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank_types.h>
#include <memory>

class Parameters;

class Model;

template<typename Eltype>
class HankArray {
	public:
		std::unique_ptr<Eltype[]> data = nullptr;

		HankArray() {
			reset(1);
		}

		HankArray(int n) {
			reset(n);
		}

		void reset(int n) {
			data.reset(new Eltype[n]);
			T = n;
		}

		Eltype& operator[](int i) {return data[i];}

		Eltype operator[](int i) const {return data[i];}

		int T;
};

class FactorQuantities {
	public:
		FactorQuantities() {}

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN;

		hank_float_type capshareY, capshareN, capfracY, capfracN, labor_Y, labor_N, labor;
};

class EquilibriumOnePeriod {
	public:
		hank_float_type alpha_Y, alpha_N, price_W, drs_Y, drs_N;
};

class Equilibrium : public HankArray<EquilibriumOnePeriod> {
	public:
		void create_initial_steady_state(const Parameters& p);

		HankArray<FactorQuantities> compute_factors(const Model& model, const std::vector<hank_float_type>& labor_occ);

		int nocc;
};



#endif