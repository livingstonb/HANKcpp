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
		const double_vector& logprodgrid = logprodgrid_;
		const double_vector& prodgrid = prodgrid_;
		const double_vector& proddist = proddist_;
		const double_matrix& prodmarkov = prodmarkov_;
		const int nb = p.nb;
		const int na = p.na;

		const double drs_Y = p.drs_Y;
		const double drs_N = p.drs_N;
		const double alpha_Y = p.alpha_Y;
		const double alpha_N = p.alpha_N;
		const double meanlabeff = p.meanlabeff;
		const double hourtarget = p.hourtarget;

		const double rb = p.rb;
		const double rborr = p.rborr;
		const double deathrate = p.deathrate;
		const bool perfectAnnuityMarkets = p.perfectAnnuityMarkets;

		double_vector get_rb_effective() const;
};

#endif