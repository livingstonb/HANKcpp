#ifndef HANK_EIGEN_DENSE
#define HANK_EIGEN_DENSE

#include <Eigen/Core>

using Eigen::seq;

typedef Eigen::Map<Eigen::MatrixXd> map_type;

typedef Eigen::VectorXd double_vector;

typedef Eigen::ArrayXd double_array;

typedef Eigen::MatrixXd double_matrix;

typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> bool_vector;

typedef double_vector grid_type;

#endif