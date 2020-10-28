#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <hank_config.h>
#include <hank.h>
#include <math.h>
#include <memory>

// Forward declarations
class Model;

class Equilibrium;

class EquilibriumTrans;

class Parameters;

namespace Upwinding
{
	class Policies;

	class ConUpwind;
}

// Class for solving the HJB
class HJB {
	private:
		Upwinding::Policies update_policies(const Equilibrium& ss);

		void update_value_fn(const Equilibrium& ss, const Upwinding::Policies& policies);

		Upwinding::ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage, double /* idioscale */) const;

		Upwinding::ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		const Model& model;

		const Parameters& p;

	public:
		HJB(const Parameters& p_, const Model& model_, const vector3dr& V);
	
		HJB(const Parameters& p_, const Model& model_, const Equilibrium& ss);

		void iterate(const Equilibrium& ss);

		void update(const EquilibriumTrans& equm);

		vector3dr V;

		int maxiter = 500;

		int dispfreq = 50;

		double vtol = 1.0e-8;

		double delta = 1.0e6;

		double dVamin = 1.0e-8;

		double dVbmin = 1.0e-8;

		double riskaver;

		std::unique_ptr<Upwinding::Policies> optimal_decisions;
};

#endif