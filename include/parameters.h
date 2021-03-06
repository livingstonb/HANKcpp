#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <hank_config.h>
#include <array>
#include <hank.h>
#include <string>

class Options;

class TaylorRule {
	public:
		hank_float_type coeff_pi = 1.25;

		hank_float_type coeff_y = 0;

		hank_float_type pers = 0.5; // Ignored if using feedback rule

		bool use_feedback_rule = true;
};

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

enum class GovBCAdjType { spending, lumptax, debt, proptax, ftpl, none, fiscal };

enum class FirmDiscountRateType { rho, rb_iss, ra_iss, rb_trans, ra_trans };

class WealthTarget
{
	public:
		enum class Type { mean, median, none };

		WealthTarget(WealthTarget::Type type_, double value_) : type(type_), value(value_) {}

		Type type;

		double value;

		bool is_mean() const {return (type == Type::mean);}

		bool is_median() const {return (type == Type::median);}
};

class Parameters : public HankBase {
	public:
		Parameters() {}

		void setup(const Options& opts);

		void update();

		std::map<std::string, hank_float_type> variables_map() const override;

		std::string title() const override {return "PARAMETER VALUES";}

		bool Borrowing = true;

		// --- GRIDS ---
		// agrid
		int na = 40;
		double amin = 0.0;
		double acurv = 0.15;
		double amax = 100;

		// bgrid
		int nb_pos = 40;
		int nb_neg = 10;
		int nb;
		double bmin = 0.0;
		double bmax = 50.0;
		double bcurv = 0.35;
		double bcurv_neg = 0.4;
		double blim = -1.0;

		// Other grid parameters
		int nocc = 4;
		int nab;

		std::string income_dir = "";

		// --- PREFERENCES ---
		double elast = 5;
		double frisch = 1.0;
		bool adjustProdGridFrisch = true;
		double adjFrischGridFrac = 0.85;
		double riskaver = 1.5;
		double rho = 0.014444;
		double prefshock = 1.0;
		bool endogLabor = true;

		// --- PRODUCTION PARAMETERS ---
		double drs_Y = 0.9;
		double alpha_Y = 0.333;
		double tfp_Y = 1.0;
		double drs_N = 0.5;
		double alpha_N = 0.333;
		double tfp_N = 1.0;
		double depreciation = 0.05 / 4.0;

		// --- OTHER FIRM-SIDE PARAMETERS ---
		bool taxHHProfitIncome = true;
		double profdistfracA = 0.333; // Fraction of profits to illiquid equity (set to alpha)
		double profdistfracB = 0.0;
		double profdistfracW = 0.667;
		double profdistfracL = 0.0;
		double priceadjcost = 100.0;
		bool make_profit_correction = true;
		double capadjcost = 0.0;
		double invadjcost = 0.0;
		FirmDiscountRateType firm_discount_rate_type = FirmDiscountRateType::ra_iss;

		// --- LABOR MARKET PARAMETERS ---
		bool imposeMaxHours = true;
		double labwedge = 1.0;
		bool scaleDisutilityIdio = false;

		// --- WITHDRAWAL COSTS ---
		bool exponential_adjcosts = false;
		double kappa_w_fc = 0.0;
		std::array<double, 5> kappa_w = {
			0.0, // kappa0_w
			0.0, // kappa1_w
			0.15, // kappa2_w
			0.01, // kappa3_w
			0.0 // kappa4_w
		};

		// --- DEPOSIT COSTS --
		DepositCostMode depositCostMode = DepositCostMode::symmetric;
		double kappa_d_fc = 0.2;
		std::array<double, 5> kappa_d = {
			0.0, // kappa0_d
			0.6, // kappa1_d
			1.0, // kappa2_d
			-1.0e10, // unused but will be set equal to kappa3_w
			0.0 // kappa4_d
		};

		// --- OTHER MODEL PARAMETERS ---
		bool oneAssetNoCapital = false;
		bool borrowing = true;
		double ssgap = 0.0;
		double meanlabeff = 1.0;
		double hourtarget = 1.0 / 3.0;
		bool perfectAnnuityMarkets = true;
		double deathrate = 0.02 / 4;
		double rb = 0.02 / 4.0;
		double borrwedge = 0.025;
		double rborr = -1.0e5; // Will be set to rb + borrwedge
		bool prodispshock = false;
		AdjustCostFnRatioMode adjCostRatioMode = AdjustCostFnRatioMode::max;
		double target_KY_ratio;
		double chi;
		double USGDPperHH = 146435 / 4.0;
		double pi = 0; // Steady state inflation target

		// --- GOVERNMENT PARAMETERS ---
		double corptax = 0.0;
		double labtax = 0.25;
		double lumptransfer = 0.05;
		GovBCAdjType adjGovBudgetConstraint = GovBCAdjType::fiscal;
		TaylorRule taylor;
		double govdebtadjwedge = 0.10;

		// --- OTHER SOLUTION PARAMETERS ---
		double cmin = 1.0e-5;
		double dmax = 100;
		double facc = 1.0e-10;
		bool solveFlexPriceTransition = false;
		bool solveStickyPriceTransition = true;

		// --- CALIBRATION TARGETS ---
		WealthTarget illiqWealthTarget = WealthTarget(WealthTarget::Type::mean, 512183);
		WealthTarget liqWealthTarget = WealthTarget(WealthTarget::Type::median, 3500);
		double targetMeanIllGuess = 512183;

		// number of time periods
		// int Ttransition = 200;
};

#endif