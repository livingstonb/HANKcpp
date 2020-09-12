#ifndef _TRANSITION_MATRIX_H
#define _TRANSITION_MATRIX_H

#include <hank_eigen_sparse.h>

class Parameters;

class Model;

namespace Upwinding { class Policies; }

sparse_matrix construct_transition_matrix(const Parameters& p, const Model& model, double ra,
	const Upwinding::Policies& policies, int iy, bool kfe);

sparse_matrix get_kfe_transition_matrix(const Parameters& p, const Model& model, double ra,
	const Upwinding::Policies& policies, int iy);

#endif