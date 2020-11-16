//
// Created by xidad on 26.10.2020.
//

#ifndef DENDRIVEV3_SOLVER_H
#define DENDRIVEV3_SOLVER_H

#include <random>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>

#include "Domain.h"
#include "BoundingRect.h"
#include "MyTypedefs.h"

using namespace myTypedefs;

namespace TVDLimitersFunctions
{
	double sgn(const double &x);

	double minmodFunction(const double &f_im1, const double &f_i, const double &f_ip1);
}

enum TVDLimitersTypes
{
	MINMOD, MC, SUPERBEE
};

class Solver
{
private:
	double h, c, tMax, g;
	int NX, NY;

	bool **domainMesh;
	double **uPrev;
	double **uNext;

	Domain *domain;
	BoundingRect *boundingRect;


private:
	double iToX(const int &i);

	double jToY(const int &j);

public:
	Solver(BoundingRect &boundingRect, Domain &domain, const double &h, const double &c);

	~Solver();

	void solve(double(&ICF)(const double &, const double &), const TVDLimitersTypes &type, const double &dt,
	           const double &tMax);

	std::vector<std::vector<double>> getNeighbours(const int &j, const int &i);

	double uWave(double(&TVDLimiterFunction)(const double &, const double &, const double &), const double &f_im1,
	             const double &f_i, const double &f_ip1);

	void exportDendrite(std::fstream &file);

	void exportComputationRegion(std::fstream &file);

	void exportData(std::fstream &file);

	void save(const std::string &name);
};


#endif //DENDRIVEV3_SOLVER_H
