#include <equilibrium.h>
#include <parameters.h>
#include <model.h>

void Equilibrium::create_initial_steady_state(const Parameters& p) {
	reset(1);
	data[0].alpha_Y = p.alpha_Y;
	data[0].alpha_N = p.alpha_N;
	data[0].price_W = 1.0 - 1.0 / p.elast;
	data[0].drs_Y = p.drs_Y;
	data[0].drs_N = p.drs_N;
	nocc = p.nocc;
}

HankArray<FactorQuantities> Equilibrium::compute_factors(const Model& model, const std::vector<hank_float_type>& labor_occ) {
	HankArray<FactorQuantities> factors_array;

	for (int it=0; it<T; ++it) {
		FactorQuantities& factors = factors_array[it];
		factors.capshareY = data[it].alpha_Y * data[it].price_W * data[it].drs_Y;
		factors.capshareN = data[it].alpha_N * (1.0 - data[it].price_W) * data[it].drs_N;
		factors.capfracY = factors.capshareY / (factors.capshareY + factors.capshareN);
		factors.capfracN = factors.capshareN / (factors.capshareY + factors.capshareN);

		factors.labshareY.resize(nocc);
		factors.labshareN.resize(nocc);
		factors.labfracY.resize(nocc);
		factors.labfracN.resize(nocc);
		factors.labor_Y = 1.0;
		factors.labor_N = 1.0;
		factors.labor = 0;
		for (int io=0; io<nocc; ++io) {
			factors.labshareY[io] = (1.0 - data[it].alpha_Y) * data[it].price_W * data[it].drs_Y * model.occYsharegrid[io];
			factors.labshareN[io] = (1.0 - data[it].alpha_N) * (1.0 - data[it].price_W) * data[it].drs_N * model.occNsharegrid[io];
			factors.labfracY[io] = factors.labshareY[io] / (factors.labshareY[io] + factors.labshareN[io]);
			factors.labfracN[io] = factors.labshareN[io] / (factors.labshareY[io] + factors.labshareN[io]);
			factors.labor_Y *= pow(factors.labfracY[io] * labor_occ[io], model.occYsharegrid[io]);
			factors.labor_N *= pow(factors.labfracN[io] * labor_occ[io], model.occNsharegrid[io]);
			factors.labor += labor_occ[io] * model.occdist[io];
		}
	}

	return factors_array;
}

