#ifndef _HANK_H
#define _HANK_H

#include <boost/multi_array.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>

typedef std::vector<double> vector;

// Eigen
using Eigen::all;
using Eigen::last;
using Eigen::seq;
typedef Eigen::Map<Eigen::MatrixXd> map_type;
typedef Eigen::VectorXd double_vector;
typedef Eigen::ArrayXd double_array;
typedef Eigen::MatrixXd double_matrix;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> bool_vector;
typedef Eigen::SparseMatrix<double> sparse_matrix;
typedef Eigen::Triplet<double> triplet_type;
typedef std::vector<triplet_type> triplet_list;

// Boost
template <typename T, size_t N>
using boost_array_type = boost::multi_array<T, N>;

template <typename T, size_t N>
using boost_array_shape = boost::array<typename boost_array_type<T, N>::index, N>;

using boost3d = boost_array_type<double, 3>;
using boost1d = boost_array_type<double, 1>;
using boost3dshape = boost_array_shape<double, 3>;
using boost1dshape = boost_array_shape<double, 1>;

typedef boost::multi_array_types::index_range range;

// Enumerations
enum class LaborType { none, sep, ghh };

enum class AdjustCostFnRatioMode { none, linear, max };

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

// Declare which type to use for grids
typedef double_vector grid_type;

#endif