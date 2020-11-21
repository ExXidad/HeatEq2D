//
// Created by xidad on 26.10.2020.
//

#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const double &h)
{
	this->h = h;
	this->boundingRect = &boundingRect;
	this->domain = &domain;

	NY = static_cast<int>(boundingRect.getYSize() / h);
	NX = static_cast<int>(boundingRect.getXSize() / h);

	domainMesh = new bool *[NY];
	for (int i = 0; i < NY; ++i)
		domainMesh[i] = new bool[NX];


	dendrite = new bool *[NY];
	for (int i = 0; i < NY; ++i)
		dendrite[i] = new bool[NX];

	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			if (boundingRect.contains(iToX(i), jToY(j))
				&&
				domain.unionContains(iToX(i), jToY(j)))
			{
				domainMesh[j][i] = true;
			}
		}
	}

	electricPotential = new double *[NY];
	for (int i = 0; i < NY; ++i)
		electricPotential[i] = new double[NX];

	electricPotentialTemporary = new double *[NY];
	for (int j = 0; j < NY; ++j)
	{
		electricPotentialTemporary[j] = new double[NX];
		for (int i = 0; i < NX; ++i)
		{
			electricPotentialTemporary[j][i] = U * (1 - static_cast<double>(j) / NY);
			electricPotentialTemporary[j][i] = accessFuncIncludingBC(j, i, electricPotentialTemporary);
		}
	}

	r = new double *[NY];
	for (int i = 0; i < NY; ++i)
		r[i] = new double[NX];

	d = new double *[NY];
	for (int i = 0; i < NY; ++i)
		d[i] = new double[NX];

	q = new double *[NY];
	for (int i = 0; i < NY; ++i)
		q[i] = new double[NX];

	electricField = new vec2d *[NY];
	for (int i = 0; i < NY; ++i)
		electricField[i] = new vec2d[NX];


	gen = std::mt19937(time(nullptr));
	randI = std::uniform_int_distribution<>(0, NX - 1);
	randJ = std::uniform_int_distribution<>(0, NY - 1);
	uid = std::uniform_int_distribution<>(0, 1);
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

	for (int j = 0; j < NY; ++j)
		delete[] electricField[j];
	delete[] electricField;

	for (int j = 0; j < NY; ++j)
		delete[] r[j];
	delete[] r;

	for (int j = 0; j < NY; ++j)
		delete[] d[j];
	delete[] d;

	for (int j = 0; j < NY; ++j)
		delete[] electricPotential[j];
	delete[] electricPotential;

	for (int j = 0; j < NY; ++j)
		delete[] electricPotentialTemporary[j];
	delete[] electricPotentialTemporary;
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

void Solver::solve(const double &fraction, const double &reactionProbability)
{
	double **x1 = new double *[NY];
	for (int j = 0; j < NY; ++j)
	{
		x1[j] = new double[NX];
		for (int i = 0; i < NX; ++i)
		{
			x1[j][i] = urd(gen);
		}
	}
	double **x2 = new double *[NY];
	for (int j = 0; j < NY; ++j)
	{
		x2[j] = new double[NX];
		for (int i = 0; i < NX; ++i)
		{
			x2[j][i] = urd(gen);
		}
	}

	applyOperatorNoB(q, x1);
	std::cout << "x2^T*A*x1\t" << scalarProduct(x2, q) << std::endl;


	applyOperatorNoB(q, x2);
	std::cout << "x1^T*A*x2\t" << scalarProduct(x1, q) << std::endl;


	for (int j = 0; j < NY; ++j)
		delete[] x1[j];
	delete[] x1;

	for (int j = 0; j < NY; ++j)
		delete[] x2[j];
	delete[] x2;

	int N = static_cast
			<int>(fraction * NX * NY);

	// Delete dendrites treated as seed particles
	for (const auto &seedParticle : seedParticles)
	{
		dendrite[seedParticle[0]][seedParticle[1]] = false;
	}

	// Find initial EF
	updateElectricPotential(pow(10, -10));
	updateElectricField();

	for (int i = 0; i < N; ++i)
	{
		vec2i particle;
		if (seedParticles.empty())
		{
			// If no seed particles remained - generate random one
			particle = {0, randI(gen)};
		}
		else
		{
			// Otherwise pick one
			particle = seedParticles[seedParticles.size() - 1];
			seedParticles.erase(seedParticles.end());
		}
		while (true)
		{
			vec2i shift = randomShift(particle[0], particle[1]);

			// If doesn't collide anything - step
			if (!collidesAnything(particle[0], particle[1]))
			{
				particle += shift;
				if (!computationAreaContains(particle[0], particle[1]))
				{
					particle = {0, randI(gen)};
				}

			}
			else

				// If reaction succeeded
			if (urd(gen) <= reactionProbability)
			{
				// Relaxation process
				if (dendrite[particle[0] + 1][particle[1]])
				{
					int n1 = -1, n2 = -1;

					if (computationAreaContains(particle[0], particle[1] + 1) &&
						!dendrite[particle[0]][particle[1] + 1])
					{
						if (!dendrite[particle[0] + 1][particle[1] + 1] &&
							dendrite[particle[0] + 2][particle[1] + 1])
						{
							n1 = neighbours4(particle[0] + 1, particle[1] + 1);
						}
					}

					if (computationAreaContains(particle[0], particle[1] - 1) &&
						!dendrite[particle[0]][particle[1] - 1])
					{
						if (!dendrite[particle[0] + 1][particle[1] - 1] &&
							dendrite[particle[0] + 2][particle[1] - 1])
						{
							n2 = neighbours4(particle[0] + 1, particle[1] - 1);
						}
					}

					if (n1 == -1 && n2 == -1)
					{ dendrite[particle[0]][particle[1]] = true; }
					else if (n2 > n1)
					{ dendrite[particle[0] + 1][particle[1] - 1] = true; }
					else if (n1 > n2)
					{ dendrite[particle[0] + 1][particle[1] + 1] = true; }
					else if (n1 == n2)
					{ dendrite[particle[0] + 1][particle[1] - 1 + uid(gen) * 2] = true; }
				}
				else
				{
					dendrite[particle[0]][particle[1]] = true;
				}

				// Progress output
				if (i % std::max(1, static_cast<int>(1. * N / 100)) == 0)
				{
					std::cout << "Progress: " << 1. * i / N * 100 << "%" << std::endl;

					// Update EF
					updateElectricPotential(0.01);
					updateElectricField();
				}
				break;
			}
			else
			{
				vec2i secondaryShift = randomShift(particle[0], particle[1]);
				vec2i tmpParticle = particle + secondaryShift;
				if (!computationAreaContains(tmpParticle[0], tmpParticle[1]))
				{
					particle = {0, randI(gen)};
				}
				else if (!dendriteOrDomainContains(tmpParticle[0], tmpParticle[1]))
				{
					particle = tmpParticle;
				}
			}
		}
	}
}

bool Solver::collidesAnything(const int &j, const int &i)
{
	return neighbours4(j, i);
}

bool Solver::dendriteOrDomainContains(const int &j, const int &i)
{
	return dendrite[j][i] || domainMesh[j][i];
}

void Solver::exportDendrite(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			int val = dendrite[j][i];
			file << val << "\t";
		}
		file << std::endl;
	}
}

void Solver::exportComputationRegion(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			int val = domainMesh[j][i];
			file << val << "\t";
		}
		file << std::endl;
	}
}

void Solver::exportData(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			if (domainMesh[j][i])
			{ file << 1 << "\t"; }
			else if (dendrite[j][i])
			{ file << 2 << "\t"; }
			else
			{ file << 0 << "\t"; }
		}
		file << std::endl;
	}
}

vec2i Solver::randomShift(const int &j, const int &i)
{
	vec2d E = electricField[j][i];
	double coeff = h * Z * mu / (D + h * Z * mu * (abs(E[0]) + abs(E[1])));
	double Px = coeff * abs(E[0]);
	double Py = coeff * abs(E[1]);
	double basicProb = (1 - Px - Py) / 4;

	std::vector<double> transitionProbabilities(4, basicProb);
	if (E[0] >= 0)
	{
		transitionProbabilities[0] += Px;
	}
	else
	{ transitionProbabilities[2] += Px; }

	if (E[1] >= 0)
	{
		transitionProbabilities[1] += Py;
	}
	else
	{ transitionProbabilities[3] += Py; }


	int dI = 0, dJ = 0;
	double randNumber = urd(gen);
	for (int k = 0; k < transitionProbabilities.size(); ++k)
	{
		if (randNumber <= transitionProbabilities[k])
		{
			if (k == 0)
			{ ++dI; }
			if (k == 1)
			{ --dJ; }
			if (k == 2)
			{ --dI; }
			if (k == 3)
			{ ++dJ; }
			break;
		}
		else
		{
			randNumber -= transitionProbabilities[k];
		}
	}
	return {dJ, dI};
}

int Solver::neighbours4(const int &j, const int &i)
{
	int neighbours = 0;
	for (int k = -1; k <= 1; ++k)
		for (int l = -1; l <= 1; ++l)
			if (abs(k) xor abs(l))
			{
				int newJ = j + k, newI = i + l;

				//check whether (j, i) touches domain or dendrite
				if (computationAreaContains(newJ, newI) &&
					(dendrite[newJ][newI] || domainMesh[newJ][newI]))
				{
					++neighbours;
				}
			}
	return neighbours;
}

vec2i Solver::rectifyJI(const int &j, const int &i)
{
	int newJ = j, newI = i;
	if (NY <= j)
	{ newJ = NY - 1; }
	else if (j < 0)
	{ newJ = 0; }

	if (NX <= i)
	{ newI = NX - 1; }
	else if (i < 0)
	{ newI = 0; }
	return {newJ, newI};
}

void Solver::randomSeed(const double &fraction)
{
	int N = static_cast<int>(fraction * NX * NY);
	for (int i = 0; i < N; ++i)
	{
		vec2i randParticle = {randJ(gen), randI(gen)};
		while (dendriteOrDomainContains(randParticle[0], randParticle[1]))
			randParticle = {randJ(gen), randI(gen)};
		dendrite[randParticle[0]][randParticle[1]] = true;
		seedParticles.emplace_back(randParticle);
	}
}

void Solver::updateElectricPotential(const double &absError)
{
	std::cout << "Solving Laplace eq. over computation area" << std::endl;
	double **tmpPtr;


	int counter = 0;

	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			electricPotentialTemporary[j][i] = 0;
		}
	}

	double rNormSq0;
	rNormSq = 0;
	applyOperatorB(r, electricPotentialTemporary);

	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			r[j][i] *= -1;
			d[j][i] = r[j][i];
			rNormSq += r[j][i] * r[j][i];
		}
	}
	rNormSq0 = rNormSq;

	while (counter <= NX * NY && rNormSq > absError * absError * rNormSq0)
	{
		applyOperatorNoB(q, d);

		alpha = rNormSq / scalarProduct(d, q);

		for (int j = 0; j < NY; ++j)
		{
			for (int i = 0; i < NX; ++i)
			{
				electricPotential[j][i] = electricPotentialTemporary[j][i] + alpha * d[j][i];
			}
		}

		tmpPtr = electricPotentialTemporary;
		electricPotentialTemporary = electricPotential;
		electricPotential = tmpPtr;

		rNormSqPrev = rNormSq;
		rNormSq = 0;

		if (counter % 50 == 0)
		{
			applyOperatorB(r, electricPotentialTemporary);
			for (int j = 0; j < NY; ++j)
			{
				for (int i = 0; i < NX; ++i)
				{
					r[j][i] *= -1;
					rNormSq += r[j][i] * r[j][i];
				}
			}
		}
		else
		{
			for (int j = 0; j < NY; ++j)
			{
				for (int i = 0; i < NX; ++i)
				{
					r[j][i] = r[j][i] - alpha * q[j][i];
					rNormSq += r[j][i] * r[j][i];
				}
			}
		}

		beta = rNormSq / rNormSqPrev;

		for (int j = 0; j < NY; ++j)
		{
			for (int i = 0; i < NX; ++i)
			{
				d[j][i] = r[j][i] + beta * d[j][i];
			}
		}

		++counter;
	}
//	for (int j = 0; j < NY; ++j)
//	{
//		for (int i = 0; i < NX; ++i)
//		{
//			if (dendriteOrDomainContains(j, i))
//			{ electricPotentialTemporary[j][i] = 0; }
//		}
//	}

//		std::cout << rNormSq0 << std::endl;
//		std::cout << rNormSq << std::endl;
//		std::cout << counter << std::endl;
//	if (counter >= NX * NY)
//	{
//		std::fstream file;
//
//		file.open("data.txt", std::ios::out);
//		exportData(file);
//		file.close();
//
//		file.open("EP.txt", std::ios::out);
//		exportPotential(file);
//		file.close();
//		throw std::bad_alloc();
//	}
}


void Solver::updateElectricField()
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			electricField[j][i] = {(accessFuncIncludingBC(j, i + 1, electricPotentialTemporary) -
									accessFuncIncludingBC(j, i, electricPotentialTemporary)) / h,
								   (accessFuncIncludingBC(j - 1, i, electricPotentialTemporary) -
									accessFuncIncludingBC(j, i, electricPotentialTemporary)) / h};
		}
	}
}

void Solver::exportPotential(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			file << electricPotential[j][i] << "\t";
		}
		file << std::endl;
	}
}

void Solver::exportR(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			file << r[j][i] << "\t";
		}
		file << std::endl;
	}
}

bool Solver::computationAreaContains(const int &j, const int &i)
{
	return (0 <= j && j < NY) && (0 <= i && i < NX);
}

void Solver::exportField(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			file << electricField[j][i][0] << "\t" << electricField[j][i][1] << "\t";
		}
		file << std::endl;
	}
}


void Solver::applyOperatorB(double **result, double **x)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (!dendriteOrDomainContains(j, i))
			{
				for (int k = -1; k <= 1; k += 2)
					for (int l = -1; l <= 1; l += 2)
					{
						int newJ = j + k, newI = i + l;

						if (newJ < 0)
						{
							result[j][i] += (U - x[j][i]) * 2;
						}
						else if (!computationAreaContains(newJ, newI))
						{
							result[j][i] += 0;
						}
						else if (dendriteOrDomainContains(newJ, newI))
						{
							result[j][i] += (0 - x[j][i]) * 2;
						}
						else
						{
							result[j][i] += x[newJ][newI] - x[j][i];
						}
					}
			}
		}
	}
}

void Solver::applyOperatorNoB(double **result, double **x)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (!dendriteOrDomainContains(j, i))
			{
				for (int k = -1; k <= 1; k += 2)
					for (int l = -1; l <= 1; l += 2)
					{
						int newJ = j + k, newI = i + l;

						if (newJ < 0)
						{
							result[j][i] += (0 - x[j][i]) * 2;
						}
						else if (!computationAreaContains(newJ, newI))
						{
							result[j][i] += 0;
						}
						else if (dendriteOrDomainContains(newJ, newI))
						{
							result[j][i] += (0 - x[j][i]) * 2;
						}
						else
						{
							result[j][i] += x[newJ][newI] - x[j][i];
						}
					}
			}
		}
	}
}

double Solver::scalarProduct(double **x, double **y)
{
	double result = 0;
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result += x[j][i] * y[j][i];
		}
	}
	return result;
}

void Solver::printArray(double **arr)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			std::cout << arr[j][i] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Solver::printArray(bool **arr)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			std::cout << arr[j][i] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

bool Solver::bcf(const int &j, const int &i)
{
	return (!computationAreaContains(j, i) || dendriteOrDomainContains(j, i) || j < 0);
}

double Solver::accessFuncIncludingBC(const int &j, const int &i, double **x)
{
	if (j < 0)
	{
		return U;
	}

	if (!computationAreaContains(j, i))
	{
		vec2i rectified = rectifyJI(j, i);
		return x[rectified[0]][rectified[1]];
	}


	if (dendriteOrDomainContains(j, i))
	{
		return 0;
	}

	return x[j][i];
}

double Solver::accessFuncNotIncludingBC(const int &j, const int &i, double **x)
{
	if (j < 0)
	{
		return 0;
	}

	if (!computationAreaContains(j, i))
	{
//		vec2i rectified = rectifyJI(j, i);
//		return x[rectified[0]][rectified[1]];
	}


	if (dendriteOrDomainContains(j, i))
	{
		return 0;
	}

	return x[j][i];
}
