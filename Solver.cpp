
#include "Solver.h"

Solver::Solver(BoundingRect &boundingRect, Domain &domain, const double &h,
               double(&tempF)(const double &, const double &),
               bool(&constTempFCond)(const double &, const double &),
               double(&fluxF)(const double &, const double &),
               bool(&constFluxFCond)(const double &, const double &),
               double(&lambdaF)(const double &, const double &),
               pair(&sourceF)(const double &, const double &, const double &),
               double T0,
               const int &nmax
)
{
    this->constTempFCond = constTempFCond;
    this->tempF = tempF;
    this->constFluxFCond = constFluxFCond;
    this->fluxF = fluxF;
    this->lambdaF = lambdaF;
    this->sourceF = sourceF;
    this->T0 = T0;
    this->h = h;
    this->boundingRect = &boundingRect;
    this->domain = &domain;
    this->nmax = nmax;

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


    T = new double *[NY];
    for (int i = 0; i < NY; ++i)
        T[i] = new double[NX];

    tmpT = new double *[NY];
    for (int j = 0; j < NY; ++j)
        tmpT[j] = new double[NX];

    r = new double *[NY];
    for (int i = 0; i < NY; ++i)
        r[i] = new double[NX];

    d = new double *[NY];
    for (int i = 0; i < NY; ++i)
        d[i] = new double[NX];

    q = new double *[NY];
    for (int i = 0; i < NY; ++i)
        q[i] = new double[NX];
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

void Solver::solve(const double tol)
{
    solveTemperature(tol);
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

void Solver::solveTemperature(const double &absError)
{
    std::cout << "Solving Heat eq. over computation area" << std::endl;

    auto start = std::chrono::system_clock::now(); //start timer

    int counter = 0;

    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i)
            tmpT[j][i] = T0;

    delNew = 0;
    applyOperatorB(r, tmpT);

//	printArray(r);

    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i)
            d[j][i] = r[j][i];

    delNew = scalarProduct(r, r);
    del0 = delNew;

    std::cout << delNew << std::endl;
    system("rm r/*");

    while (counter <= NX * NY && delNew > absError * absError * NX * NY && counter < nmax)
    {
        if (counter % 20 == 0)
            std::cout << delNew << std::endl;
//		if (counter%1==0&&counter<300)
//		{
//			std::fstream file;
//			file.open("r/"+std::to_string(counter) + "_r.txt", std::ios::out);
//			exportR(file);
//			file.close();
//		}

        applyOperatorNoB(q, d);

        alpha = delNew / scalarProduct(d, q);

        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i)
                T[j][i] = tmpT[j][i] + alpha * d[j][i];

        std::swap(tmpT, T);

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
    } else { std::cout << "Converged within: " << counter << " iterations" << std::endl; }

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

void Solver::applyOperatorB(double **result, double **x)
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
                            double tmp = tempF(newX, newY);

                            if (constTempFCond(newX, newY))
                                result[j][i] += (tempF(newX, newY) - x[j][i]);
                            else if (!computationAreaContains(newJ, newI) || !domainMesh[newJ][newI])
                                result[j][i] += 0;
                            else if (constFluxFCond(newX, newY))
                                result[j][i] += fluxF(newX, newY) * h;
                            else
                                result[j][i] += (x[newJ][newI] - x[j][i]);
                        }
            }
            pair source = sourceF(iToX(i), jToY(j), x[j][i]);
            result[j][i] =
                    (source.first - source.second * x[j][i]) / lambdaF(iToX(i), jToY(j)) * h * h
                    + result[j][i];
        }
}

void Solver::applyOperatorNoB(double **result, double **x)
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
                            double tmp = tempF(newX, newY);

                            if (constTempFCond(newX, newY))
                                result[j][i] += (0. - x[j][i]);
                            else if (!computationAreaContains(newJ, newI) || !domainMesh[newJ][newI])
                                result[j][i] += 0;
                            else if (constFluxFCond(newX, newY))
                                result[j][i] += fluxF(newX, newY) * h;
                            else
                                result[j][i] += (x[newJ][newI] - x[j][i]);
                        }
            }
            pair source = sourceF(iToX(i), jToY(j), x[j][i]);
            result[j][i] = -(
                    (-source.second * x[j][i]) / lambdaF(iToX(i), jToY(j)) * h * h
                    + result[j][i]
            );
        }
}

double Solver::scalarProduct(double **x, double **y)
{
    double result = 0;
    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i)
            result += x[j][i] * y[j][i];
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