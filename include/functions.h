#ifndef _HANK_FUNCTIONS_H
#define _HANK_FUNCTIONS_H

#include <hank_config.h>
#include <math.h>

enum class FnType { function, deriv, deriv_inv };

namespace HankFunctions {

// template<FnType ftype=FnType::function>
// hank_float_type utility(hank_float_type c, hank_float_type prefshock, hank_float_type riskaver) {
// 	if (riskaver == 1.0)
// 		return prefshock * log(c);
// 	else
// 		return prefshock * pow(c, 1.0 - riskaver) / (1.0 - riskaver);
// }

// template<>
// hank_float_type utility<FnType::deriv>(hank_float_type c, hank_float_type prefshock, hank_float_type riskaver) {
// 	if (riskaver == 1.0)
// 		return prefshock * log(c);
// 	else
// 		return prefshock * pow(c, 1.0 - riskaver) / (1.0 - riskaver);
// }

inline hank_float_type utility(hank_float_type c, hank_float_type prefshock, hank_float_type riskaver) {
	if (riskaver == 1.0)
		return prefshock * log(c);
	else
		return prefshock * pow(c, 1.0 - riskaver) / (1.0 - riskaver);
}

inline hank_float_type utility1(hank_float_type c, hank_float_type prefshock, hank_float_type riskaver) {
	return prefshock * pow(c, -riskaver);
}

inline hank_float_type utility1inv(hank_float_type u, hank_float_type prefshock, hank_float_type riskaver) {
	return pow(u / prefshock, -1.0 / riskaver);
}

inline hank_float_type labor_disutility(hank_float_type h, hank_float_type frisch, hank_float_type chi) {
	return chi * pow(h, 1 + 1.0 / frisch) / (1 + 1.0 / frisch);
}

inline hank_float_type labor_disutility1(hank_float_type h, hank_float_type frisch, hank_float_type chi) {
	return chi * pow(h, 1.0 / frisch);
}

inline hank_float_type labor_disutility1inv(hank_float_type du, hank_float_type frisch, hank_float_type chi) {
	return pow(du / chi, frisch);
}

inline hank_float_type capadjcost(hank_float_type x, hank_float_type adjcost, hank_float_type deprec) {
	return 0.5 * adjcost * pow(x - deprec, 2);
}

inline hank_float_type capadjcost1(hank_float_type x, hank_float_type adjcost, hank_float_type deprec) {
	return adjcost * (x - deprec);
}

inline hank_float_type capadjcost1inv(hank_float_type x, hank_float_type adjcost, hank_float_type deprec) {
	return x / adjcost + deprec;
}

}

#endif