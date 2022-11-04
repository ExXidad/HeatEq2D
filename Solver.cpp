
#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const float &h, const float &U)
{
	this->U = U;
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
				domain.unionContains(iToX(i), jToY(j)) || j == NY - 1)
			{
				domainMesh[j][i] = true;
			}
		}
	}

	electricPotential = new float *[NY];
	for (int i = 0; i < NY; ++i)
		electricPotential[i] = new float[NX];

	electricPotentialTemporary = new float *[NY];
	for (int j = 0; j < NY; ++j)
	{
		electricPotentialTemporary[j] = new float[NX];
		for (int i = 0; i < NX; ++i)
		{
			electricPotentialTemporary[j][i] = U * (1 - static_cast<float>(j) / NY);
		}
	}

	r = new float *[NY];
	for (int i = 0; i < NY; ++i)
		r[i] = new float[NX];

	d = new float *[NY];
	for (int i = 0; i < NY; ++i)
		d[i] = new float[NX];

	q = new float *[NY];
	for (int i = 0; i < NY; ++i)
		q[i] = new float[NX];

	electricFieldI = new float *[NY];
	for (int i = 0; i < NY; ++i)
		electricFieldI[i] = new float[NX];

	electricFieldJ = new float *[NY];
	for (int i = 0; i < NY; ++i)
		electricFieldJ[i] = new float[NX];


	gen = std::mt19937(123);
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
		delete[] electricFieldJ[j];
	delete[] electricFieldJ;

	for (int j = 0; j < NY; ++j)
		delete[] electricFieldI[j];
	delete[] electricFieldI;

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

float Solver::iToX(const int &i)
{
	return boundingRect->getSize()[0][0] + h / 2 + i * h;
}

float Solver::jToY(const int &j)
{
	return boundingRect->getSize()[1][1] - (h / 2 + j * h);
}

void Solver::addNucleus(const int &j, const int &i)
{
	dendrite[j][i] = true;
}

void Solver::solve(const float &fraction, const float &reactionProbability)
{
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
		if (reachedTopEdgeFlag){
			break;
		}
		vec2i particle(2);
		if (seedParticles.empty())
		{
			// If no seed particles remained - generate random one
			particle = {0, randI(gen)};
		}
		else
		{
			// Otherwise pick one
			do
			{
				particle = seedParticles[seedParticles.size() - 1];
				seedParticles.erase(seedParticles.end());
				if (seedParticles.empty())
				{
					break;
				}
			} while (dendriteOrDomainContains(particle[0], particle[1]));
		}

		while (true)
		{
			randomShift(shift, particle[0], particle[1]);

			// If doesn't collide anything - step
			if (!collidesAnything(particle[0], particle[1]))
			{
				particle = {particle[0] + shift[0], particle[1] + shift[1]};
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
					int n0 = neighbours4(particle[0], particle[1]), n1 = -1, n2 = -1;

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

					if (n0 >= n1 && n0 >= n2)
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
					if (particle[0] == 0){
						std::cout << "Dendrite reached top edge" << std::endl;
						reachedTopEdgeFlag = true;
						break;
					}
				}

				// Progress output
				if (i % std::max(1, static_cast<int>(1. * N / 100)) == 0)
				{
					std::cout << std::endl << "Progress: " << 1. * i / N * 100 << "%" << std::endl;

					// Update EF
					updateElectricPotential(pow(10, -8));
					updateElectricField();

					if (saveProgressFlag)
					{
						// Export dendrite
						std::fstream file;
						file.open(std::to_string(i), std::ios::out);
						exportDendrite(file);
						file.close();
					}
				}
				break;
			}
				// If reaction failed continue moving
			else
			{
				randomShift(secondaryShift, particle[0], particle[1]);
				vec2i tmpParticle = {particle[0] + secondaryShift[0], particle[1] + secondaryShift[1]};
				// Delete particle if it left computation area
				if (!computationAreaContains(tmpParticle[0], tmpParticle[1]))
				{
					particle = {0, randI(gen)};
				}
					// Make a step
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

void Solver::randomShift(vec2i &shiftVar, const int &j, const int &i)
{
	float Ei = electricFieldI[j][i], Ej = electricFieldJ[j][i];
	float coeff = h * Z * mu / (D + h * Z * mu * (std::abs(Ei) + std::abs(Ej)));
	float Pi = coeff * std::abs(Ei);
	float Pj = coeff * std::abs(Ej);
	float basicProb = (1 - Pi - Pj) / 4;

	std::vector<float> transitionProbabilities(4, basicProb);
	if (Ei >= 0)
	{
		transitionProbabilities[0] += Pi;
	}
	else
	{ transitionProbabilities[2] += Pi; }

	if (Ej >= 0)
	{
		transitionProbabilities[3] += Pj;
	}
	else
	{ transitionProbabilities[1] += Pj; }


	int dI = 0, dJ = 0;
	float randNumber = urd(gen);
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
	shiftVar = {dJ, dI};
}

int Solver::neighbours4(const int &j, const int &i)
{
	int neighbours = 0;
	for (int k = -1; k <= 1; ++k)
		for (int l = -1; l <= 1; ++l)
		{
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
		}
	return neighbours;
}

void Solver::randomSeed(const float &fraction)
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

void Solver::updateElectricPotential(const float &absError)
{
	std::cout << "Solving Laplace eq. over computation area" << std::endl;
	float **tmpPtr;

	auto start = std::chrono::system_clock::now(); //start timer

	int counter = 0;

//#pragma omp parallel for
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			electricPotentialTemporary[j][i] = 0;
		}
	}
//#pragma omp barrier

	delNew = 0;
	applyOperatorB(r, electricPotentialTemporary);

//#pragma omp parallel for shared(r, d) reduction(+: delNew)
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			r[j][i] *= -1;
			d[j][i] = r[j][i] / 4;
			delNew += d[j][i] * r[j][i];
		}
	}
//#pragma omp barrier

	del0 = delNew;

	while (counter <= NX * NY && delNew > absError * absError * del0)
	{
		applyOperatorNoB(q, d);

		alpha = delNew / scalarProduct(d, q);

//#pragma omp parallel for
		for (int j = 0; j < NY; ++j)
		{
			for (int i = 0; i < NX; ++i)
			{
				electricPotential[j][i] = electricPotentialTemporary[j][i] + alpha * d[j][i];
			}
		}
//#pragma omp barrier

		tmpPtr = electricPotentialTemporary;
		electricPotentialTemporary = electricPotential;
		electricPotential = tmpPtr;

		delOld = delNew;
		delNew = 0;

		if (counter % 50 == 0)
		{
			applyOperatorB(r, electricPotentialTemporary);
//#pragma omp parallel for shared(r) reduction(+: delNew)
			for (int j = 0; j < NY; ++j)
			{
				for (int i = 0; i < NX; ++i)
				{
					r[j][i] *= -1;
					delNew += r[j][i] * r[j][i] / 4;
				}
			}
//#pragma omp barrier
		}
		else
		{
//#pragma omp parallel for shared(r, q) reduction(+: delNew)
			for (int j = 0; j < NY; ++j)
			{
				for (int i = 0; i < NX; ++i)
				{
					r[j][i] = r[j][i] - alpha * q[j][i];
					delNew += r[j][i] * r[j][i] / 4;
				}
			}
//#pragma omp barrier
		}

		beta = delNew / delOld;

//#pragma omp parallel for shared(r, d)
		for (int j = 0; j < NY; ++j)
		{
			for (int i = 0; i < NX; ++i)
			{
				d[j][i] = r[j][i] / 4 + beta * d[j][i];
			}
		}
//#pragma omp barrier

		++counter;
	}

	if (counter >= NX * NY)
	{
		std::cout << "Didn't converge within: " << NX * NY << " iterations" << std::endl;
	}
	else
	{ std::cout << "Converged within: " << counter << " iterations" << std::endl; }

	auto end = std::chrono::system_clock::now();//end timer
	std::chrono::duration<float> elapsed = end - start;
	std::cout << "Laplace took: " << elapsed.count() << "s" << std::endl;
}


void Solver::updateElectricField()
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			electricFieldJ[j][i] = 0;
			electricFieldI[j][i] = 0;
			if (!dendriteOrDomainContains(j, i))
			{
				float hx = 0, hy = 0;

				if (j - 1 < 0)
				{
					electricFieldJ[j][i] -= U;
					hy += 1. / 2;
				}
				else if (dendriteOrDomainContains(j - 1, i))
				{
					electricFieldJ[j][i] -= 0;
					hy += 1. / 2;
				}
				else
				{
					electricFieldJ[j][i] -= electricPotentialTemporary[j - 1][i];
					hy += 1;
				}


				if (j + 1 >= NY)
				{
					electricFieldJ[j][i] += electricPotentialTemporary[j][i];
					hy += 1. / 2;
				}
				else if (dendriteOrDomainContains(j + 1, i))
				{
					electricFieldJ[j][i] += 0;
					hy += 1. / 2;
				}
				else
				{
					electricFieldJ[j][i] += electricPotentialTemporary[j + 1][i];
					hy += 1;
				}


				if (i + 1 >= NX)
				{
					electricFieldI[j][i] += electricPotentialTemporary[j][i];
					hx += 1;
				}
				else if (dendriteOrDomainContains(j, i + 1))
				{
					electricFieldI[j][i] += 0;
					hx += 1. / 2;
				}
				else
				{
					electricFieldI[j][i] += electricPotentialTemporary[j][i + 1];
					hx += 1;
				}


				if (i - 1 < 0)
				{
					electricFieldI[j][i] -= electricPotentialTemporary[j][i];
					hx += 1. / 2;
				}
				else if (dendriteOrDomainContains(j, i - 1))
				{
					electricFieldI[j][i] -= 0;
					hx += 1. / 2;
				}
				else
				{
					electricFieldI[j][i] -= electricPotentialTemporary[j][i - 1];
					hx += 1;
				}


				electricFieldI[j][i] = -electricFieldI[j][i] / hx / h;
				electricFieldJ[j][i] = -electricFieldJ[j][i] / hy / h;
			}
		}
	}
}


void Solver::exportPotential(std::fstream &file)
{
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			file << iToX(i) << "\t" << jToY(j) << "\t" << electricPotential[j][i] << "\t" << std::endl;;
		}
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
			file << iToX(i) << "\t" << jToY(j) << "\t" << electricFieldI[j][i] << "\t" << electricFieldJ[j][i]
				 << std::endl;
		}
	}
}


void Solver::applyOperatorB(float **result, float **x)
{
//#pragma omp parallel for shared(x, result)
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (!dendriteOrDomainContains(j, i))
			{
				for (int k = -1; k <= 1; ++k)
					for (int l = -1; l <= 1; ++l)
					{
						if (abs(k) xor abs(l))
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
//#pragma omp barrier
}

void Solver::applyOperatorNoB(float **result, float **x)
{
//#pragma omp parallel for shared(x, result)
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (!dendriteOrDomainContains(j, i))
			{
				for (int k = -1; k <= 1; ++k)
					for (int l = -1; l <= 1; ++l)
					{
						if (abs(k) xor abs(l))
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
}

float Solver::scalarProduct(float **x, float **y)
{
	float result = 0;
//#pragma omp parallel for shared(r, d) reduction(+: result)
	for (int j = 0; j < NY; ++j)
	{
		for (int i = 0; i < NX; ++i)
		{
			result += x[j][i] * y[j][i];
		}
	}
	return result;
}

void Solver::printArray(float **arr)
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

float Solver::getH() const
{
	return h;
}

int Solver::getNx() const
{
	return NX;
}

int Solver::getNy() const
{
	return NY;
}

float Solver::getU() const
{
	return U;
}

float Solver::getMu() const
{
	return mu;
}

float Solver::getD() const
{
	return D;
}

void Solver::setSaveProgressFlag(bool saveProgressFlag)
{
	Solver::saveProgressFlag = saveProgressFlag;
}

//#pragma clang diagnostic pop