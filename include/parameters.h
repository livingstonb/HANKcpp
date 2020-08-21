#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <iostream>
#include <options.h>
#include <hank.h>

class Parameters {
	public:
		Parameters() {}

		void setup(const Options& opts);

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
		double blim = -1.0;

		// Other grid parameters
		int nprod = 15;
		int nocc = 6;
		int ny;
		int naby;
		int nab;

		// --- PREFERENCES ---
		double elast = 1.1;
		double frisch = 1.0;
		double riskaver = 1.0;
		double rho = 0.01444;
		double prefshock = 1.0;
		LaborType laborsupply = LaborType::sep;

		// --- PRODUCTION PARAMETERS ---
		double drs_Y = 1.0;
		double alpha_Y = 0.333;
		double tfp_Y = 1.0;
		double drs_N = 0.0;
		double alpha_N = 0.333;
		double tfp_N = 1.0;
		double depreciation = 0.07 / 4.0;

		// --- OTHER FIRM-SIDE PARAMETERS ---
		double taxHHProfitIncome = 1.0;
		double profdistfracA = 0.333; // Fraction of profits to illiquid equity (set to alpha)
		double profdistfracB = 0.0;
		double profdistfracW = 0.667;
		double profdistfracL = 0.0;
		double priceadjcost = 100.0;

		// --- LABOR MARKET PARAMETERS ---
		bool imposeMaxHours = true;
		double labwedge = 1.0;
		bool scaleDisutilityIdio = false;

		// --- WITHDRAWAL COSTS ---
		bool exponential_adjcosts = false;
		double kappa_w_fc = 0.0;
		double kappa_w[5] = {
			0.0, // kappa0_w
			0.04336, // kappa1_w
			0.40176, // kappa2_w
			0.03 * 2.92 / 4.0, // kappa3_w
			0.0 // kappa4_w
		};

		// --- DEPOSIT COSTS --
		double kappa_d_fc = 0.2;
		double kappa_d[5] = {
			0.0, // kappa0_d
			0.6, // kappa1_d
			1.0, // kappa2_d
			-1.0e10, // unused but will be set equal to kappa3_w
			0.0 // kappa4_d
		};

		// --- OTHER MODEL PARAMETERS ---
		bool borrowing = true;
		// Steady state output gap
		double ssgap = 0.0;
		double meanlabeff = 1.0;
		double hourtarget = 1.0 / 3.0;
		bool perfectAnnuityMarkets = true;
		double deathrate = 1.0 / (4.0 * 45.0);
		double rb = 0.02 / 4.0;
		double borrwedge = 0.025;
		double rborr = -1.0e5; // Will be set to rb + borrwedge
		double corptax = 0.0;
		double labtax = 0.30;
		double lumptransfer = 0.05;
		bool prodispshock = false;
		bool prodDispScaleDisutility = true;
		AdjustCostFnRatioMode adjCostRatioMode = AdjustCostFnRatioMode::max;

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