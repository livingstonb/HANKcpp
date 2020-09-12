#include <parameters.h>
#include <algorithm>
#include <vector>
#include <math.h>

void Parameters::setup(const Options& opts) {

	if ( opts.fast ) {
		nb_pos = 10;
		nb_neg = 5;
		na = 10;
	}

	if ( borrowing )
		nb = nb_pos + nb_neg;
	else
		nb = nb_pos;

	ny = nprod * nocc;
	naby = na * nb * ny;
	nab = na * nb;

	rborr = rb + borrwedge;

	kappa_w[1] = pow((1.0-kappa_w[0]) * (1.0+kappa_w[2]), -1.0 / kappa_w[2]);
	kappa_d[3] = kappa_w[3];
	switch ( opts.depositCostMode ) {
		case DepositCostMode::symmetric:
			// Set kappa_d equal to kappa_w
			kappa_d_fc = kappa_w_fc;
			kappa_d = kappa_w;
			break;
		case DepositCostMode::no_deposit_cost:
			kappa_d_fc = 0;
			kappa_d[0] = 0;
			kappa_d[1] = 100; // This is the optimal depost rate, i.e. higher for lower cost
			kappa_d[2] = 1;
			kappa_d[4] = 0;
			break;
		case DepositCostMode::custom:
			// Leave deposit adj cost parameters at current values
			break;
	}
};