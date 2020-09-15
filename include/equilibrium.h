#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

class Parameters;

class Model;

class Equilibrium {
	public:
		Equilibrium(const Parameters& p, const Model& model);

		const Parameters& p;

		const Model& model;
};

#endif