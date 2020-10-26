#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <hank_config.h>
#include <hank.h>
#include <upwinding.h>
#include <math.h>

// Forward declarations
class Model;

class Equilibrium;

class Parameters;

// Class for solving the HJB
class HJB {
	private:
		Upwinding::Policies update_policies(const Equilibrium& ss);

		void update_value_fn(const Equilibrium& ss, const Upwinding::Policies& policies);

		Upwinding::ConUpwind optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		Upwinding::ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const;

		Upwinding::ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

	public:
		HJB(const Model& model_, const Equilibrium& ss);

		void iterate(const Equilibrium& ss);

		void print_variables() const;

		const Model& model;

		const Parameters& p;

		vector3dr V;

		int maxiter = 500;

		int dispfreq = 50;

		double vtol = 1.0e-8;

		double delta = 1.0e6;

		double dVamin = 1.0e-8;

		double dVbmin = 1.0e-8;

		double riskaver;

		Upwinding::Policies optimal_decisions;
};

#endif