#include <model.h>

Model::Model(Parameters params) {
	p = params;
	make_grids();
}

void Model::make_grids() {
	bgrid = powerSpacedGrid<vector>(p.nb, p.bmin, p.bmax, p.bcurv);
	agrid = powerSpacedGrid<vector>(p.na, p.amin, p.amax, p.acurv);

	auto occgrids = occupationGrid<vector>(p);
	occgrid = occgrids.first;
	occdist = occgrids.second;
}