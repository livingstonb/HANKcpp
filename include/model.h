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

class ModelBase
{
	public:
		ModelBase(const Parameters params, const std::string& income_dir) : p(params) {
			make_grids(p);
			create_income_process(income_dir);
		}

		const Parameters p;

		double_vector bgrid_;
		double_vector agrid_;
		double_vector occgrid_;
		double_vector occdist_;

		double_vector logprodgrid_;
		double_vector prodgrid_;
		double_vector proddist_;
		double_matrix prodmarkov_;

		void make_grids(const Parameters& params);
		void create_income_process(const std::string& income_dir);

		int nb_ = p.nb;
		int na_ = p.na;
};

class Model : private ModelBase {
	public:
		Model(const Parameters p, const std::string& income_dir)
			: ModelBase(p, income_dir) {};
		const double_vector& bgrid = bgrid_;
		const double_vector& agrid = agrid_;
		const double_vector& occgrid = occgrid_;
		const double_vector& logprodgrid = logprodgrid_;
		const double_vector& prodgrid = prodgrid_;
		const double_vector& proddist = proddist_;
		const double_matrix& prodmarkov = prodmarkov_;
		const int nb = p.nb;
		const int na = p.na;
		const double drs_Y = p.drs_Y;
};

#endif