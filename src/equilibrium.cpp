#include <equilibrium.h>
#include <parameters.h>
#include <model.h>

Equilibrium::Equilibrium(const Parameters& p_, const Model& model_)
	: p(p_), model(model_) {}

void Equilibrium::initial_ss() {
	profit = 0.0;
	taxrev = p.labtax * model.wage_occ.dot(model.labor_occ)
}