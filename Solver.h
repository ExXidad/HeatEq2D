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
	double h;
	int NX, NY;

	double U = 4;
	double mu = 6.4 * pow(10,-8), D = 1.648 * pow(10,-9), Z = 1;

	bool **domainMesh;
	bool **dendrite;
	double **electricPotential;
	double **electricPotentialTemporary;
	double **r, **d, **q;
	double alpha, beta, delOld, delNew, del0;
	bool firstEPUpdRun = true;

	double **electricFieldI;
	double **electricFieldJ;

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

	bool collidesAnything(const int &j, const int &i);

	__attribute__((always_inline)) bool dendriteOrDomainContains(const int &j, const int &i);

	int neighbours4(const int &j, const int &i);

	void randomShift(vec2i &shiftVar, const int &j, const int &i);

	void updateElectricPotential(const double &absError);

	void updateElectricField();

	double scalarProduct(double **x, double **y);

	__attribute__((always_inline)) bool computationAreaContains(const int &j, const int &i);

	void applyOperatorB(double **result, double **x);

	void applyOperatorNoB(double **result, double **x);

public:
	Solver(BoundingRect &boundingRect, Domain &domain, const double &h);

	~Solver();

	void addNucleus(const int &j, const int &i);

	void randomSeed(const double &fraction);

	void solve(const double &fraction, const double &reactionProbability);

	void exportDendrite(std::fstream &file);

	void exportComputationRegion(std::fstream &file);

	void exportData(std::fstream &file);

	void exportPotential(std::fstream &file);

	void exportField(std::fstream &file);

	void exportR(std::fstream &file);

	void printArray(double **arr);

	void printArray(bool **arr);
};


#endif //DENDRIVEV3_SOLVER_H
