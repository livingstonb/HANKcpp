#ifndef _UPWINDING_H
#define _UPWINDING_H

#include <hank_config.h>
#include <hank_types.h>

namespace Upwinding {

struct ConUpwind {
	double c, h, s, Hc = -1.0e12;
	bool valid = false;
};

struct DepositUpwind {
	double d, Hd = -1.0e12;
	bool valid = false;

	bool at_least_as_good_as(const DepositUpwind& du) const {
		return ( (!du.valid) | (Hd >= du.Hd) );
	}
};

class Policies {
	public:
		Policies(const std::vector<int> dims) : c(dims), h(dims), s(dims), d(dims), u(dims) {};

		vector3dr c, h, s, d, u;

		void update_c(int ia, int ib, int iy, const ConUpwind& uwF, const ConUpwind& uwB, const ConUpwind& uw0);

		void update_d(int ia, int ib, int iy, const DepositUpwind& uFB,
			const DepositUpwind& uBF, const DepositUpwind& uBB);
};

}

#endif