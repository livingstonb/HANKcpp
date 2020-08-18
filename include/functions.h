#ifndef _HANK_FUNCTIONS_H
#define _HANK_FUNCTIONS_H

#include <cmath>

namespace HankFunctions {

inline double utility(double c, double prefshock, double riskaver) {
	if (riskaver == 1.0)
		return prefshock * log(c);
	else
		return prefshock * pow(c, 1.0 - riskaver) / (1.0 - riskaver);
}

inline double utility1(double c, double prefshock, double riskaver) {
	return prefshock * pow(c, -riskaver);
}

inline double utility1inv(double u, double prefshock, double riskaver) {
	return pow(u / prefshock, -1.0 / riskaver);
}


inline double labor_disutility(double h, double frisch, double chi) {
	return chi * pow(h, 1 + 1.0 / frisch) / (1 + 1.0 / frisch);
}

inline double labor_disutility1(double h, double frisch, double chi) {
	return chi * pow(h, 1.0 / frisch);
}

inline double labor_disutility1inv(double du, double frisch, double chi) {
	return pow(du / chi, frisch);
}

}

#endif