
#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const float &h,
			   double(&tempF)(const double &, const double &),
			   bool(&constTempFCond)(const double &, const double &),
			   double(&lambdaF)(const double &, const double &),
			   double(&sourceF)(const double &, const double &),
			   double T0
)
{
	this->constTempFCond = constTempFCond;
	this->tempF = tempF;
	this->lambdaF = lambdaF;
	this->sourceF = sourceF;
	this->T0 = T0;
	this->h = h;
	this->boundingRect = &boundingRect;
	this->domain = &domain;

	NY = static_cast<int>(boundingRect.getYSize() / h) + 1;
	NX = static_cast<int>(boundingRect.getXSize() / h) + 1;

	domainMesh = new bool *[NY];
	for (int j = 0; j < NY; ++j)
	{
		domainMesh[j] = new bool[NX];
		for (int i = 0; i < NX; ++i)
			if (boundingRect.contains(iToX(i), jToY(j))
				&&
				domain.unionContains(iToX(i), jToY(j)))
				domainMesh[j][i] = true;
	}


	T = new float *[NY];
	for (int i = 0; i < NY; ++i)
		T[i] = new float[NX];

	tmpT = new float *[NY];
	for (int j = 0; j < NY; ++j)
	{
		tmpT[j] = new float[NX];
		for (int i = 0; i < NX; ++i)
			tmpT[j][i] = tempF(iToX(i), jToY(j));
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
}

Solver::~Solver()
{
	for (int j = 0; j < NY; ++j)
		delete[] domainMesh[j];
	delete[] domainMesh;

	for (int j = 0; j < NY; ++j)
		delete[] r[j];
	delete[] r;

	for (int j = 0; j < NY; ++j)
		delete[] d[j];
	delete[] d;

	for (int j = 0; j < NY; ++j)
		delete[] T[j];
	delete[] T;

	for (int j = 0; j < NY; ++j)
		delete[] tmpT[j];
	delete[] tmpT;
}

double Solver::iToX(const int &i)
{
	return boundingRect->getSize()[0][0] + i * h;
}

double Solver::jToY(const int &j)
{
	return boundingRect->getSize()[1][1] - j * h;
}

void Solver::solve()
{
	updateElectricPotential(pow(10, -10));
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

void Solver::updateElectricPotential(const float &absError)
{
	std::cout << "Solving Heat eq. over computation area" << std::endl;
	float **tmpPtr;

	auto start = std::chrono::system_clock::now(); //start timer

	int counter = 0;

	for (int j = 0; j < NY; ++j)
		for (int i = 0; i < NX; ++i)
			tmpT[j][i] = 100;

	delNew = 0;
	applyOperatorB(r, tmpT);

	for (int j = 0; j < NY; ++j)
		for (int i = 0; i < NX; ++i)
			d[j][i] = r[j][i];

	delNew = scalarProduct(r, r);
	del0 = delNew;

	while (counter <= NX * NY && delNew > absError * absError * del0)
	{
		applyOperatorNoB(q, d);

		alpha = delNew / scalarProduct(d, q);

		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
				T[j][i] = tmpT[j][i] + alpha * d[j][i];

		tmpPtr = tmpT;
		tmpT = T;
		T = tmpPtr;

		if (counter % 50 == 0)
			applyOperatorB(r, tmpT);
		else
			for (int j = 0; j < NY; ++j)
				for (int i = 0; i < NX; ++i)
					r[j][i] = r[j][i] - alpha * q[j][i];

		delOld = delNew;
		delNew = scalarProduct(r, r);

		beta = delNew / delOld;

		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
				d[j][i] = r[j][i] + beta * d[j][i];

		++counter;
	}

	if (counter >= NX * NY)
	{
		std::cout << "Didn't converge within: " << NX * NY << " iterations" << std::endl;
	} else
	{ std::cout << "Converged within: " << counter << " iterations" << std::endl; }

	auto end = std::chrono::system_clock::now();//end timer
	std::chrono::duration<float> elapsed = end - start;
	std::cout << "Solver took: " << elapsed.count() << "s" << std::endl;
}

void Solver::exportTemp(std::fstream &file, const bool includeCoord)
{
	if (includeCoord)
		for (int j = 0; j < NY; ++j)
			for (int i = 0; i < NX; ++i)
				file << iToX(i) << "\t" << jToY(j) << "\t" << T[j][i] << "\n";
	else
		for (int j = 0; j < NY; ++j)
		{
			for (int i = 0; i < NX; ++i)
				file << T[j][i] << "\t";
			file << "\n";
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

void Solver::applyOperatorB(float **result, float **x)
{
	for (int j = 0; j < NY; ++j)
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (domainMesh[j][i])
			{
				for (int k = -1; k <= 1; ++k)
					for (int l = -1; l <= 1; ++l)

						if (abs(k) xor abs(l))
						{
							int newJ = j + k, newI = i + l;
							double newX = iToX(newI), newY = jToY(newJ);

							if (constTempFCond(newX, newY))
								result[j][i] += (tempF(newX, newY) - x[j][i]) * 2;
							else if (!computationAreaContains(newJ, newI) || !domainMesh[newJ][newI])
								result[j][i] += 0;
							else
								result[j][i] += x[newJ][newI] - x[j][i];
						}
			}
			result[j][i] = -sourceF(iToX(i), jToY(j)) - result[j][i] * 0.25 * lambdaF(iToX(i), jToY(j));
		}
}

void Solver::applyOperatorNoB(float **result, float **x)
{
	for (int j = 0; j < NY; ++j)
		for (int i = 0; i < NX; ++i)
		{
			result[j][i] = 0;

			if (domainMesh[j][i])
			{
				for (int k = -1; k <= 1; ++k)
					for (int l = -1; l <= 1; ++l)

						if (abs(k) xor abs(l))
						{
							int newJ = j + k, newI = i + l;

							if (constTempFCond(iToX(newI), jToY(newJ)))
								result[j][i] += (0. - x[j][i]) * 2;
							else if (!computationAreaContains(newJ, newI) || !domainMesh[newJ][newI])
								result[j][i] += 0;
							else
								result[j][i] += x[newJ][newI] - x[j][i];
						}
			}
			result[j][i] = result[j][i] * 0.25 * lambdaF(iToX(i), jToY(j));
		}
}

float Solver::scalarProduct(float **x, float **y)
{
	float result = 0;
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

void Solver::setSaveProgressFlag(bool saveProgressFlag)
{
	Solver::saveProgressFlag = saveProgressFlag;
}

//#pragma clang diagnostic pop