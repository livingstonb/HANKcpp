
#include <hank_eigen_sparse.h>

// Forward declarations
class Model;

class HJB;

class StationaryDist {
	public:
		StationaryDist() {}

		void compute(const Model& model, const HJB& hjb, const std::vector<sparse_matrix>& A);

		double delta = 1.0e6;
};