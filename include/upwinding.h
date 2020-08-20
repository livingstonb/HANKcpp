#ifndef _UPWINDING_H
#define _UPWINDING_H

#include <hank.h>

namespace Upwinding {

struct ConUpwind {
	double c, h, s, Hc = 1.0e-10;
	bool valid = false;
};

struct DepositUpwind {
	double d, Hd = 1.0e-10;
	bool valid = false;

	bool at_least_as_good_as(const DepositUpwind& du) const {
		return ( (!du.valid) | (Hd >= du.Hd) );
	}
};

class Policies {
	public:
		Policies(const boost3dshape& dims) : c(dims), h(dims), s(dims), d(dims), u(dims) {};

		boost3d c, h, s, d, u;

		void update_c(int ia, int ib, int iy, const ConUpwind& uwF, const ConUpwind& uwB, const ConUpwind& uw0) {
			bool not_backward = (!uwB.valid) | (uwF.Hc >= uwB.Hc);
			bool not_forward = (!uwF.valid) | (uwB.Hc >= uwF.Hc);
			bool forward_better_than_nothing = uwF.Hc >= uw0.Hc;
			bool backward_better_than_nothing = uwB.Hc >= uw0.Hc;

			ConUpwind uw_selected;
			if ( uwF.valid & not_backward & forward_better_than_nothing )
				uw_selected = uwF;
			else if ( uwB.valid & not_forward & backward_better_than_nothing )
				uw_selected = uwB;
			else
				uw_selected = uw0;

			c[ia][ib][iy] = uw_selected.c;
			h[ia][ib][iy] = uw_selected.h;
			s[ia][ib][iy] = uw_selected.s;
		}

		void update_d(int ia, int ib, int iy, const DepositUpwind& uFB,
			const DepositUpwind& uBF, const DepositUpwind& uBB) {
			bool chooseFB = uFB.valid & uFB.at_least_as_good_as(uBF) & uFB.at_least_as_good_as(uBB);
			bool chooseBF = uBF.valid & uBF.at_least_as_good_as(uFB) & uBF.at_least_as_good_as(uBB);
			bool chooseBB = uBB.valid & uFB.at_least_as_good_as(uBF) & uBB.at_least_as_good_as(uFB);

			if ( chooseFB )
				d[ia][ib][iy] = uFB.d;
			else if ( chooseBF )
				d[ia][ib][iy] = uBF.d;
			else if ( chooseBB )
				d[ia][ib][iy] = uBB.d;
			else if ( (!uFB.valid) & (!uBF.valid) & (!uBB.valid) )
				d[ia][ib][iy] = 0;
			else
				throw "Error while upwinding deposits";
		}
};

}

#endif