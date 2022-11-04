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

	float U;
	float mu = 6.4 * pow(10,-8), D = 1.648 * pow(10,-9), Z = 1;

	bool **domainMesh;
	bool **dendrite;
	float **electricPotential;
	float **electricPotentialTemporary;
	float **r, **d, **q;
	float alpha, beta, delOld, delNew, del0;

	bool firstEPUpdRun = true;
	bool saveProgressFlag = true;
	bool reachedTopEdgeFlag = false;

	float **electricFieldI;
	float **electricFieldJ;

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
	float iToX(const int &i);

	float jToY(const int &j);

	bool collidesAnything(const int &j, const int &i);

	__attribute__((always_inline)) bool dendriteOrDomainContains(const int &j, const int &i);

	int neighbours4(const int &j, const int &i);

	void randomShift(vec2i &shiftVar, const int &j, const int &i);

	void updateElectricPotential(const float &absError);

	void updateElectricField();

	float scalarProduct(float **x, float **y);

	__attribute__((always_inline)) bool computationAreaContains(const int &j, const int &i);

	void applyOperatorB(float **result, float **x);

	void applyOperatorNoB(float **result, float **x);

public:
	float getH() const;

	int getNx() const;

	int getNy() const;

	float getU() const;

	float getMu() const;

	float getD() const;

	void setSaveProgressFlag(bool saveProgressFlag);

public:
	Solver(BoundingRect &boundingRect, Domain &domain, const float &h, const float &U);

	~Solver();

	void addNucleus(const int &j, const int &i);

	void randomSeed(const float &fraction);

	void solve(const float &fraction, const float &reactionProbability);

	void exportDendrite(std::fstream &file);

	void exportComputationRegion(std::fstream &file);

	void exportData(std::fstream &file);

	void exportPotential(std::fstream &file);

	void exportField(std::fstream &file);

	void exportR(std::fstream &file);

	void printArray(float **arr);

	void printArray(bool **arr);
};


#endif //DENDRIVEV3_SOLVER_H
