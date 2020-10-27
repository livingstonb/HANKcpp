#ifndef _HANK_FUNCTIONS_H
#define _HANK_FUNCTIONS_H

#include <math.h>

enum class FnType { function, deriv, deriv_inv };

namespace ModelFunctions {

template<typename T1, typename T2, typename T3>
auto  utility(T1 c, T2 prefshock, T3 riskaver) -> decltype(prefshock * log(c))
{
	if (riskaver == 1.0)
		return prefshock * log(c);
	else
		return prefshock * pow(c, 1.0 - riskaver) / (1.0 - riskaver);
}

template<typename T1, typename T2, typename T3>
auto utility1(T1 c, T2 prefshock, T3 riskaver) -> decltype(prefshock * pow(c, -riskaver))
{
	return prefshock * pow(c, -riskaver);
}

template<typename T1, typename T2, typename T3>
auto utility1inv(T1 u, T2 prefshock, T3 riskaver) -> decltype(pow(u / prefshock, -1.0 / riskaver))
{
	return pow(u / prefshock, -1.0 / riskaver);
}

template<typename T1, typename T2, typename T3>
auto labor_disutility(T1 h, T2 frisch, T3 chi) -> decltype(chi * pow(h, 1 + 1.0 / frisch) / (1 + 1.0 / frisch))
{
	return chi * pow(h, 1 + 1.0 / frisch) / (1 + 1.0 / frisch);
}

template<typename T1, typename T2, typename T3>
auto labor_disutility1(T1 h, T2 frisch, T3 chi) -> decltype(chi * pow(h, 1.0 / frisch))
{
	return chi * pow(h, 1.0 / frisch);
}

template<typename T1, typename T2, typename T3>
auto labor_disutility1inv(T1 du, T2 frisch, T3 chi) -> decltype(pow(du / chi, frisch))
{
	return pow(du / chi, frisch);
}

template<typename T1, typename T2, typename T3>
auto capadjcost(T1 x, T2 adjcost, T3 deprec) -> decltype(adjcost * pow(x - deprec, 2))
{
	return 0.5 * adjcost * pow(x - deprec, 2);
}

template<typename T1, typename T2, typename T3>
auto capadjcost1(T1 x, T2 adjcost, T3 deprec) -> decltype(adjcost * (x - deprec))
{
	return adjcost * (x - deprec);
}

template<typename T1, typename T2, typename T3>
auto capadjcost1inv(T1 x, T2 adjcost, T3 deprec) -> decltype(x / adjcost + deprec)
{
	return x / adjcost + deprec;
}

template<typename T1, typename T2, typename T3>
auto priceadjcost(T1 pi, T2 output, T3 adjcostparam) -> decltype((adjcostparam / 2.0) * pow(pi, 2) * output)
{
	return (adjcostparam / 2.0) * pow(pi, 2) * output;
}

}

#endif