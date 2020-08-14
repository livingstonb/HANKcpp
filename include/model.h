#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include <parameters.h>
#include <procedures.h>
#include <hank.h>

std::vector<double> read_matrix(const std::string& file_loc);

std::size_t find_multiple(const std::string& line, int pos);

void fix_rounding(double_matrix& mat);

double_vector compute_grid_deltas(const double_vector& grid);

class ModelBase
{
	public:
		ModelBase(const Parameters p, const std::string& income_dir) {
			make_grids(p);
			create_income_process(income_dir, p);
			create_combined_variables();
		}

		double_vector bgrid_;
		double_vector bdelta_;
		double_vector agrid_;
		double_vector adelta_;
		double_vector occgrid_;
		double_vector occdist_;

		double_vector logprodgrid_;
		double_vector prodgrid_;
		double_vector proddist_;
		double_matrix prodmarkov_;

		double_vector yprodgrid_;
		double_vector yoccgrid_;
		double_vector ydist_;
		double_matrix ymarkov_;

		int nocc_;
		int nprod_;

		void make_grids(const Parameters& p);
		void create_income_process(const std::string& income_dir, const Parameters& p);
		void create_combined_variables();
};

class Model : private ModelBase {
	public:
		Model(const Parameters p_, const std::string& income_dir)
			: ModelBase(p_, income_dir), p(p_) {};

		const Parameters p;

		const double_vector& bgrid = bgrid_;
		const double_vector& bdelta = bdelta_;
		const double_vector& agrid = agrid_;
		const double_vector& adelta = adelta_;
		const double_vector& occgrid = occgrid_;
		const double_vector& occdist = occdist_;

		const double_vector& logprodgrid = logprodgrid_;
		const double_vector& prodgrid = prodgrid_;
		const double_vector& proddist = proddist_;
		const double_matrix& prodmarkov = prodmarkov_;

		const double_vector& yprodgrid = yprodgrid_;
		const double_vector& yoccgrid = yoccgrid_;
		const double_vector& ydist = ydist_;
		const double_matrix& ymarkov = ymarkov_;

		const int nb = p.nb;
		const int na = p.na;
		const int nocc = nocc_;
		const int nprod = nprod_;

		double_vector get_rb_effective() const;
};

#endif