//
// Created by xidad on 26.10.2020.
//

#ifndef DENDRIVEV3_SOLVER_H
#define DENDRIVEV3_SOLVER_H

#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>

#include "Domain.h"
#include "BoundingRect.h"
#include "MyTypedefs.h"

using namespace myTypedefs;

class Solver
{
private:
	double h;
	int NX, NY;

	double U = 10;
	double mu = 0.1, D = 1, Z = 1;

	bool **domainMesh;
	bool **dendrite;
	double **electricPotential;
	double **electricPotentialTemporary;
	double **r, **d, **q;
	double alpha, beta, rNormSqPrev, rNormSq;

	vec2d **electricField;

	Domain *domain;
	BoundingRect *boundingRect;

	std::vector<vec2i> seedParticles{};

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

	bool dendriteOrDomainContains(const int &j, const int &i);

	int neighbours4(const int &j, const int &i);

	vec2i rectifyJI(const int &j, const int &i);

	vec2i randomShift(const int &j, const int &i);

	void updateElectricPotential(const double &absError);

	void updateElectricField();

	double scalarProduct(double **x, double **y);

	bool computationAreaContains(const int &j, const int &i);

	void applyOperator(double **result, double **x, double (Solver::*accessFunction)(const int &, const int &, double **x));

	bool bcf(const int &j, const int &i);

	double accessFuncIncludingBC(const int &j, const int &i, double **x);

	double accessFuncNotIncludingBC(const int &j, const int &i, double **x);

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

	void printArray(double **arr);

	void printArray(bool **arr);
};


#endif //DENDRIVEV3_SOLVER_H
