//
// Created by xidad on 26.10.2020.
//

#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const double &h, const double &c)
{
	this->h = h;
	this->c = c;
	this->boundingRect = &boundingRect;
	this->domain = &domain;

	NY = static_cast<int>(boundingRect.getYSize() / h) - 1;
	NX = static_cast<int>(boundingRect.getXSize() / h) - 1;

	domainMesh = new bool *[NY];
	for (int i = 0; i < NX; ++i)
		domainMesh[i] = new bool[NX];


	uPrev = new double *[NY];
	for (int i = 0; i < NX; ++i)
		uPrev[i] = new double[NX];

	uNext = new double *[NY];
	for (int i = 0; i < NX; ++i)
		uNext[i] = new double[NX];

	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			if (domain.contains(iToX(i), jToY(j)))
				domainMesh[j][i] = true;
		}
	}
}

Solver::~Solver()
{
	for (int j = 0; j < NY; ++j)
		delete[] domainMesh[j];
	delete[] domainMesh;

	for (int j = 0; j < NY; ++j)
		delete[] uPrev[j];
	delete[] uPrev;

	for (int j = 0; j < NY; ++j)
		delete[] uNext[j];
	delete[] uNext;
}

double Solver::iToX(const int &i)
{
	return boundingRect->getSize()[0][0] + h / 2 + i * h;
}

double Solver::jToY(const int &j)
{
	return boundingRect->getSize()[1][1] - (h / 2 + j * h);
}

void Solver::solve(double(&ICF)(const double &, const double &), const TVDLimitersTypes &type, const double &dt,
                   const double &tMax)
{
	int N = static_cast<int>(tMax / dt) + 1;
	this->tMax = tMax;
	this->g = c * dt / h;

	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			uPrev[j][i] = ICF(iToX(i), jToY(j));
		}
	}

	double
	(*TVDLimiterFunction)(const double &, const double &, const double &) = &TVDLimitersFunctions::minmodFunction;

	switch (type) {
		case MINMOD:
			TVDLimiterFunction = &TVDLimitersFunctions::minmodFunction;
			break;

		case MC:
			break;

		case SUPERBEE:
			break;
	}

	save(std::to_string(0));

	for (int tIteration = 0; tIteration < N; ++tIteration) {
		for (int j = 0; j < NY; ++j) {
			for (int i = 0; i < NX; ++i) {
				if (!domainMesh[j][i]) {
					std::vector<std::vector<double>> neighboursP = getNeighbours(j, i);
					double uWaveXp = uWave(*TVDLimiterFunction, neighboursP[0][0], uPrev[j][i], neighboursP[1][0]);
					double uWaveYp = uWave(*TVDLimiterFunction, neighboursP[1][0], uPrev[j][i], neighboursP[1][1]);

					std::vector<std::vector<double>> neighboursM = getNeighbours(j - 1, i - 1);
					double uWaveXm = uWave(*TVDLimiterFunction, neighboursM[0][0], uPrev[j][i], neighboursM[1][0]);
					double uWaveYm = uWave(*TVDLimiterFunction, neighboursM[1][0], uPrev[j][i], neighboursM[1][1]);

					uNext[j][i] = uPrev[j][i] - g * ((uWaveXp - uWaveXm) / 2 + (uWaveYp - uWaveYm) / 2);
				}
			}
		}
		double **tmpPointer = uPrev;
		uPrev = uNext;
		uNext = tmpPointer;

		if (tIteration % static_cast<int>(N / 10) == 0) {
			std::cout << "Progress: " << 1. * tIteration / N * 100 << "%" << std::endl;
			save(std::to_string(tIteration));
		}
	}
}


void Solver::exportDendrite(std::fstream &file)
{
	for (int j = 0; j < NY; ++j) {
		for (int i = 0; i < NX; ++i) {
			int val = uPrev[j][i];
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
			file << uPrev[j][i] << "\t";
		}
		file << std::endl;
	}
}

std::vector<std::vector<double>> Solver::getNeighbours(const int &j, const int &i)
{
	double uNextXNeighbour, uPrevXNeighbour;
	double uNextYNeighbour, uPrevYNeighbour;

	if (i + 1 > NX - 1 || j < 0 || domainMesh[j][i + 1])
		uNextXNeighbour = 0;
	else
		uNextXNeighbour = uPrev[j][i + 1];

	if (i - 1 < 0 || j < 0 || domainMesh[j][i - 1])
		uPrevXNeighbour = 0;
	else
		uPrevXNeighbour = uPrev[j][i - 1];

	if (j + 1 > NY - 1 || i < 0 || domainMesh[j + 1][i])
		uNextYNeighbour = 0;
	else
		uNextYNeighbour = uPrev[j + 1][i];

	if (j - 1 < 0 || i < 0 || domainMesh[j - 1][i])
		uPrevYNeighbour = 0;
	else
		uPrevYNeighbour = uPrev[j - 1][i];

	return {{uPrevXNeighbour, uNextXNeighbour},
	        {uPrevYNeighbour, uNextYNeighbour}};
}

double Solver::uWave(double(&TVDLimiterFunction)(const double &, const double &, const double &), const double &f_im1,
                     const double &f_i, const double &f_ip1)
{
	return f_i + (1 - g) / 2 * TVDLimiterFunction(f_im1, f_i, f_ip1);
}

void Solver::save(const std::string &name)
{
	std::fstream file;
	file.open(name, std::ios::out);
	exportData(file);
	file.close();
}


double TVDLimitersFunctions::sgn(const double &x)
{
	if (x > 0) return 1;
	else if (x < 0) return -1;
	return 0;
}

double TVDLimitersFunctions::minmodFunction(const double &f_im1, const double &f_i, const double &f_ip1)
{
	return std::min(std::abs(f_ip1 - f_i), std::abs(f_i - f_im1)) * sgn(f_ip1 - f_i);
}