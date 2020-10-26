//
// Created by xidad on 26.10.2020.
//

#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const double &h)
{
	this->h = h;
	this->boundingRect = &boundingRect;
	this->domain = &domain;

	NY = static_cast<int>(boundingRect.getYSize() / h) - 1;
	NX = static_cast<int>(boundingRect.getXSize() / h) - 1;

	domainMesh = new bool *[NY];
	for (int i = 0; i < NX; ++i)
		domainMesh[i] = new bool[NX];


	dendrite = new bool *[NY];
	for (int i = 0; i < NX; ++i)
		dendrite[i] = new bool[NX];

	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			if (domain.unionContains(iToX(i), jToY(j)))
				domainMesh[j][i] = true;
		}
	}


	gen = std::mt19937(time(nullptr));
	randI = std::uniform_int_distribution<>(0, NX - 1);
	urd = std::uniform_real_distribution<>(0, 1);
}

Solver::~Solver()
{
	for (int j = 0; j < NY; ++j)
		delete[] domainMesh[j];
	delete[] domainMesh;

	for (int j = 0; j < NY; ++j)
		delete[] dendrite[j];
	delete[] dendrite;
}

double Solver::iToX(const int &i)
{
	return boundingRect->getSize()[0][0] + h / 2 + i * h;
}

double Solver::jToY(const int &j)
{
	return boundingRect->getSize()[1][1] - (h / 2 + j * h);
}

void Solver::addNucleus(const int &j, const int &i)
{
	dendrite[j][i] = true;
}

void Solver::solve(const int &N, const double &reactionProbability)
{
	for (int i = 0; i < N; ++i) {
		vec2i particle{0, randI(gen)};
		while (true) {
			vec2i shift = randomShift();

			if (!collides(particle[0], particle[1])) {
				particle += shift;
				if (!boundingRect->contains(iToX(particle[1]), jToY(particle[0]))) particle = {0, randI(gen)};
			}

			if (collides(particle[0], particle[1]))
				if (urd(gen) <= reactionProbability) {
					dendrite[particle[0]][particle[1]] = true;
					if (i % std::max(1, static_cast<int>(1. * N / 100)) == 0)
						std::cout << "Progress: " << 1. * i / N * 100 << "%" << std::endl;
					break;
				} else {
					vec2i secondaryShift = randomShift();
					while (contains((particle + secondaryShift)[0],
									(particle + secondaryShift)[1])) secondaryShift = randomShift();
				}
		}
	}
}

bool Solver::collides(const int &j, const int &i)
{
	for (int k = -1; k <= 1; ++k)
		for (int l = -1; l <= 1; ++l)
			if (abs(k) xor abs(l)) {
				int newJ = j + k, newI = i + l;

				if (NY <= newJ) newJ = NY - 1;
				if (newJ < 0) newJ = 0;

				if (NX <= newI) newI = NY - 1;
				if (newI < 0) newI = 0;

				if (dendrite[newJ][newI] || domainMesh[newJ][newI]) return true;
			}
	return false;
}

bool Solver::contains(const int &j, const int &i)
{
	return dendrite[j][i] || domainMesh[j][i];
}

void Solver::exportDendrite(std::fstream &file)
{
	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			int val = dendrite[j][i];
			file << val << "\t";
		}
		file << std::endl;
	}
}

void Solver::exportComputationRegion(std::fstream &file)
{
	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			int val = domainMesh[j][i];
			file << val << "\t";
		}
		file << std::endl;
	}
}

void Solver::exportData(std::fstream &file)
{
	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			if (domainMesh[j][i]) file << 1 << "\t";
			else if (dendrite[j][i]) file << 2 << "\t";
			else file << 0 << "\t";
		}
		file << std::endl;
	}
}

vec2i Solver::randomShift()
{
	int dI = 0, dJ = 0;
	double randNumber = urd(gen);
	for (int k = 0; k < transitionProbabilities.size(); ++k) {
		if (randNumber <= transitionProbabilities[k]) {
			if (k == 0) --dJ;
			if (k == 1) --dI;
			if (k == 2) ++dJ;
			if (k == 3) ++dI;
			break;
		} else {
			randNumber -= transitionProbabilities[k];
		}
	}
	return vec2i(dJ, dI);
}
