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
	// Delete dendrites treated as seed particles
	for (const auto &seedParticle : seedParticles)
	{
		dendrite[seedParticle[0]][seedParticle[1]] = false;
	}

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
			vec2i shift = randomShift();

			// If doesn't collide anything - step
			if (!collidesAnything(particle[0], particle[1]))
			{
				particle += shift;
				if (!boundingRect->contains(iToX(particle[1]), jToY(particle[0])))
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

					if (boundingRect->contains(iToX(particle[1] + 1), jToY(particle[0])) &&
						!dendrite[particle[0]][particle[1] + 1])
					{
						if (!dendrite[particle[0] + 1][particle[1] + 1] &&
							dendrite[particle[0] + 2][particle[1] + 1])
						{
							n1 = neighbours4(particle[0] + 1, particle[1] + 1);
						}
					}

					if (boundingRect->contains(iToX(particle[1] - 1), jToY(particle[0])) &&
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
				}
				break;
			}
			else
			{
				vec2i secondaryShift = randomShift();
				vec2i tmpParticle = particle + secondaryShift;
				if (!boundingRect->contains(iToX(tmpParticle[1]), jToY(tmpParticle[0])))
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

vec2i Solver::randomShift()
{
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
				if (boundingRect->contains(iToX(newI), jToY(newJ)) &&
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
	if (NY <= newJ)
	{ newJ = NY - 1; }
	if (newJ < 0)
	{ newJ = 0; }

	if (NX <= newI)
	{ newI = NY - 1; }
	if (newI < 0)
	{ newI = 0; }
	return {newJ, newI};
}

void Solver::randomSeed(const int &N)
{
	for (int i = 0; i < N; ++i)
	{
		vec2i randParticle = {randJ(gen), randI(gen)};
		while (dendriteOrDomainContains(randParticle[0], randParticle[1]))
			randParticle = {randJ(gen), randI(gen)};
		dendrite[randParticle[0]][randParticle[1]] = true;
		seedParticles.emplace_back(randParticle);
	}
}
