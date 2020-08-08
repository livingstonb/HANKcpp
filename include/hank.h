#ifndef _HANK_H
#define _HANK_H

#include <boost/multi_array.hpp>
#include <Eigen/Core>

typedef std::vector<double> vector;

typedef Eigen::Map<Eigen::MatrixXd> map_type;

// Boost
template <size_t N>
using array_type = boost::multi_array<double, N>;

template <size_t N>
using array_shape = boost::array<typename array_type<N>::index, N>;

typedef boost::multi_array_types::index_range range;

#endif