#include <model.h>

Model::Model(Parameters params) {
	p = params;
	make_grids();
}

void Model::make_grids() {
	bgrid = powerSpacedGrid(p.nb, p.bmin, p.bmax, p.bcurv);
	agrid = powerSpacedGrid(p.na, p.amin, p.amax, p.acurv);

	auto occgrids = occupationGrid(p);
	occgrid = occgrids.first;
	occdist = occgrids.second;
}