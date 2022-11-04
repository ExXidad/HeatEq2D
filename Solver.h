#ifndef DENDRIVEV3_SOLVER_H
#define DENDRIVEV3_SOLVER_H

#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <chrono>

#include "Domain.h"
#include "BoundingRect.h"
#include "MyTypedefs.h"

using namespace myTypedefs;

class Solver
{
private:
	float h;
	int NX, NY;

	bool **domainMesh;
	float **T;
	float **tmpT;

	bool (*constTempFCond)(const double &, const double &);

	double (*tempF)(const double &, const double &);

	double (*lambdaF)(const double &, const double &);

	double (*sourceF)(const double &, const double &, const double &);

	double T0;
	float **r, **d, **q;
	float alpha, beta, delOld, delNew, del0;

	bool firstEPUpdRun = true;
	bool saveProgressFlag = true;
	bool reachedTopEdgeFlag = false;

	Domain *domain;
	BoundingRect *boundingRect;

	std::vector<vec2i> seedParticles{};
	vec2i shift, secondaryShift;

private:
	std::mt19937 gen;
	std::uniform_int_distribution<> randI;
	std::uniform_int_distribution<> randJ;
	std::uniform_int_distribution<> uid;
	std::uniform_real_distribution<> urd;

private:
	double iToX(const int &i);

	double jToY(const int &j);

	void updateElectricPotential(const float &absError);

	float scalarProduct(float **x, float **y);

	__attribute__((always_inline)) bool computationAreaContains(const int &j, const int &i);

	void applyOperatorB(float **result, float **x);

	void applyOperatorNoB(float **result, float **x);

public:
	float getH() const;

	int getNx() const;

	int getNy() const;

	void setSaveProgressFlag(bool saveProgressFlag);

public:
	Solver(BoundingRect &boundingRect, Domain &domain, const float &h,
		   double(&tempF)(const double &, const double &),
		   bool(&constTempF)(const double &, const double &),
		   double(&lambdaF)(const double &, const double &),
		   double(&sourceF)(const double &, const double &, const double &),
		   const double T0 = 1.
	);

	~Solver();

	void solve(const double tol = 1e-5);

	void exportComputationRegion(std::fstream &file);

	void exportTemp(std::fstream &file, const bool includeCoord = false);

	void exportR(std::fstream &file);

	void printArray(float **arr);

	void printArray(bool **arr);
};


#endif //DENDRIVEV3_SOLVER_H
