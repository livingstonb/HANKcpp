#ifndef _HANK_H
#define _HANK_H

#include <boost/multi_array.hpp>
#include <Eigen/Core>

typedef std::vector<double> vector;

// template <typename T>
// using smart_ptr = std::shared_ptr<T>;

typedef Eigen::Map<Eigen::MatrixXd> map_type;
typedef Eigen::VectorXd double_vector;
typedef Eigen::MatrixXd double_matrix;

// Boost
template <typename T, size_t N>
using boost_array_type = boost::multi_array<T, N>;

template <typename T, size_t N>
using boost_array_shape = boost::array<typename boost_array_type<T, N>::index, N>;

typedef boost::multi_array_types::index_range range;

#endif