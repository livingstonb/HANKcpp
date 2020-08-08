
#include <iostream>
#include <procedures.h>
#include <parameters.h>
#include <model.h>
#include <initial_steady_state.h>
#include <mkl.h>
// #include <cblas.h>

class Grids {
	public:
		Grids() {};

		double acurv = 0.15;
		double bcurv_pos = 0.35;
		double bcurv_neg = 0.4;
		double amax = 2000.0;
		double bmax = 40.0;
};



int main () {
	Parameters params = Parameters();
	Model model = Model(params);
	
	solve_initial_steady_state(model);

	std::vector<double> vec1 {1.0, 2.0, 3.0};
	std::vector<double> vec2 {1.0, -10.0, 1.0};

	// double dot = ddot(3, vec1.data(), 1, vec2.data(), 1);
	double dot = vdot(vec1, vec2);
	std::cout << dot << '\n';
}





// !PARAMETERS FOR GRID CONSTRUCTION
// real(8), parameter :: agridparam = 0.15 !0.9 !0.15	!for a: approaches linear as goes to 1, approaches L shaped as goes to 0
// real(8), parameter :: bgridparam = 0.35 !0.25 !0.35	!for b pos: approaches linear as goes to 1, approaches L shaped as goes to 0
// real(8), parameter :: bgridparamNEG = 0.4 !for b neg: approaches linear as goes to 1, approaches L shaped as goes to 0
// real(8), parameter :: amax  = 2000.0 !0.1 !2000.0 !multiple of quarterly output
// real(8), parameter :: bmax  = 40.0 !500.0 !40.0

// !OTHER PARAMETERS
// real(8), parameter :: cmin = 1.0e-5 !minimum consumption for natural borrowing limit
// real(8), parameter :: dmax = 1.0e10 !maximum deposit rate, for numerical stability while converging
// real(8), parameter :: facc = 1.0e-10 !1.0e-6

// integer, parameter :: Ttransition = 78 !200 !no. time steps for the transition (each step is can be a different number of time units)

// END MODULE Parameters