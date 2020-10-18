#ifndef _TRANSITION_MATRIX_H
#define _TRANSITION_MATRIX_H

class Parameters;

class Model;

namespace Upwinding { class Policies; }

class SparseMatContainer;

SparseMatContainer construct_transition_matrix(const Parameters& p, const Model& model, double ra,
	double illprice, const Upwinding::Policies& policies, int iy, bool kfe);

SparseMatContainer get_kfe_transition_matrix(const Parameters& p, const Model& model, double ra,
	double illprice, const Upwinding::Policies& policies, int iy);

#endif