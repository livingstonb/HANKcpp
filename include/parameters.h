#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <iostream>

class Parameters {
	public:
		Parameters();

		bool Borrowing = true;

		// --- GRIDS ---
		// agrid
		int na = 40;
		double amin = 0.0;
		double acurv = 0.2;
		double amax = 100.0;

		// bgrid
		int nb_pos = 40;
		int nb_neg = 10;
		int nb;
		double bmin = 0.0;
		double bmax = 50.0;
		double bcurv = 0.2;

		// Other grid parameters
		int nprod = 15;
		int nocc = 6;
		int ny;
		int naby;
		int nab;

		// --- PREFERENCES ---
		// Elasticity of substitution
		double elast = 1.1;

		// --- PRODUCTION PARAMETERS ---
		double drs_Y = 1.0;
		double alpha_Y = 0.333;
		double tfp_Y = 1.0;
		double drs_N = 0.0;
		double alpha_N = 0.333;
		double tfp_N = 1.0;
		double depreciation = 0.07 / 4.0;

		// --- OTHER FIRM-SIDE PARAMETERS --
		// Fraction of profits to illiquid equity (set to alpha)
		double profdistfracA;

		// --- OTHER MODEL PARAMETERS ---
		// Steady state output gap
		double ssgap = 0.0;
		double priceadjcost = 100.0;

		// --- OTHER SOLUTION PARAMETERS ---
		double cmin = 1.0e-5;
		double dmax = 1.0e10;
		double facc = 1.0e-10;

		// --- CALIBRATION TARGETS ---
		double targetMeanIll = 11.68;
		double targetMeanLiq = 0.8;

		// number of time periods
		int Ttransition = 200;
};

#endif