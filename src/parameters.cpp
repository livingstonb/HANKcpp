#include <parameters.h>

// namespace {
// 	void copy_array(const double[]& from, double[]& to, int n) {
// 		for (int i=0; i<j; ++i)
// 			to[i] = from[i];
// 	}
// }

Parameters::setup(const Options& opts) {
	nb = nb_pos + nb_neg;

	ny = nprod * nocc;
	naby = na * nb * ny;
	nab = na * nb;

	profdistfracA = alpha_Y;

	rborr = rb + borrwedge;

	kappa_d[3] = kappa_w[3];
	switch ( opts.depositCostMode ) {
		case depositCostMode::symmetric:
			// Set kappa_d equal to kappa_w
			kappafc_d = kappafc_w;
			std::copy(std::begin(kappa_w), std::end(kappa_w), std::begin(kappa_d));
			break;
		case depositCostMode::no_deposit_cost:
			kappafc_d = 0;
			kappa_d[0] = 0;
			kappa_d[1] = 100; // This is the optimal depost rate, i.e. higher for lower cost
			kappa_d[2] = 1;
			kappa_d[4] = 0;
			break;
		case depositCostMode::custom:
			// Leave deposit adj cost parameters at current values
			break;
	}
};