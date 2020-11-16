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

enum TVDLimitersTypes
{
	MINMOD, MC, SUPERBEE
};

class Solver
{
private:
	double h, c, CN, dt;
	int NX, NT, N, saveTStep;

	double **u;
	double *uTempNext, *uTempPrevious;

	BoundingRect *boundingRect;


private:
	double iToX(const int &i);

	double jToY(const int &j);

	double uf(const int &i);

private:
	double sgn(const double &x);

	double minmodFunction(const int &i);

	double uWavePlusHalf(double (Solver::*TVDLimiterFunction)(const int &), const int &i);

public:
	Solver(BoundingRect &boundingRect, const double &h, const double &c, const double &saveTRate,
	       const double &CN = 0.1);

	~Solver();

	void solve(double(&ICF)(const double &, const double &), const TVDLimitersTypes &type);


	void exportData(std::fstream &file);

	void save(const std::string &name);
};


#endif //DENDRIVEV3_SOLVER_H
