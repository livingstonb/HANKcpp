#include <model.h>

Model::Model(Parameters params) {
	p = params;
	make_grids();
}

void Model::make_grids() {
	bgrid = PowerSpacedGrid(p.nb, p.bmin, p.bmax, p.bcurv);
	agrid = PowerSpacedGrid(p.na, p.amin, p.amax, p.acurv);
}