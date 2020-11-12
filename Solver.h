//
// Created by xidad on 26.10.2020.
//

#ifndef DENDRIVEV3_SOLVER_H
#define DENDRIVEV3_SOLVER_H

#include <random>
#include <ctime>
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

	bool **domainMesh;
	bool **dendrite;

	Domain *domain;
	BoundingRect *boundingRect;

	std::vector<double> transitionProbabilities{0.246667,0.246667,0.246667,0.26};
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

	vec2i randomShift();

public:
	Solver(BoundingRect &boundingRect, Domain &domain, const double &h);

	~Solver();

	void addNucleus(const int &j, const int &i);

	void randomSeed(const int &N);

	void solve(const int &N, const double &reactionProbability);

	void exportDendrite(std::fstream &file);

	void exportComputationRegion(std::fstream &file);

	void exportData(std::fstream &file);
};


#endif //DENDRIVEV3_SOLVER_H
