#include <parameters.h>

Parameters::Parameters() {
	nb = nb_pos + nb_neg;

	ny = nprod * nocc;
	naby = na * nb * ny;
	nab = na * nb;

	profdistfracA = alpha_Y;

	rborr = rb + borrwedge;
};