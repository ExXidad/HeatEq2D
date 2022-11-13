#include <iostream>
#include <random>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"
#include <algorithm>

constexpr double h = 1e-2;
constexpr double d = 5 * h;
constexpr double xmin = 0., xmax = 2.;
constexpr double ymin = 0., ymax = 1.;
constexpr double T1 = 100.;
constexpr double T2 = 150.;
BoundingRect boundingRect(xmin, xmax, ymin, ymax);

BoundingRect br1(xmin - d, xmin + 0.25, ymax - 0.25, ymax + d);
BoundingRect br2((xmin + xmax) * 0.5 - 0.25 / 2, (xmin + xmax) * 0.5 + 0.25 / 2,
                 (ymin + ymax) * 0.5 - 0.25 / 2, (ymin + ymax) * 0.5 + 0.25 / 2);
BoundingRect br3(xmax - 0.25, xmax + d, ymax - 0.25, ymax + d);
BoundingRect br4(xmax - 0.25, xmax + d, ymin - d, ymin + 0.25);
BoundingRect br5(xmin - d, xmin + 0.25, ymin - d, ymin + 0.25);

BoundingRect bottomRect(xmin - 0.5, xmax + 0.5, -0.5 + ymin, ymin);
BoundingRect topRect(xmin - 0.5, xmax + 0.5, ymax, ymax + 0.5);

bool domainFunction(const double &x, const double &y)
{
    return !(
            br3.contains(x, y) && !br3.innerShellContains(2. * h, x, y)
            ||
            br4.contains(x, y) && !br4.innerShellContains(2. * h, x, y)
            ||
            br5.contains(x, y) && !br5.innerShellContains(2. * h, x, y)
    );
}

bool constTempFCond(const double &x, const double &y)
{
    return (
            br1.innerShellContains(2. * h, x, y)
            ||
            br2.innerShellContains(2. * h, x, y)
            ||
            br5.innerShellContains(2. * h, x, y)
//            bottomRect.contains(x, y)
//            ||
//            topRect.contains(x, y)
//            ||
//    !boundingRect.contains(x, y)
    );
}

double tempF(const double &x, const double &y)
{
    if (br1.innerShellContains(2. * h, x, y))
        return T1;
    else if (br2.innerShellContains(2. * h, x, y))
        return T2;
    else if (br5.innerShellContains(2. * h, x, y))
        if (y > x)
            return T1;
        else return T2;
//	if (bottomRect.contains(x, y)) return 50.;
//    else if (topRect.contains(x, y)) return 150;
}

bool constFluxFCond(const double &x, const double &y)
{
    return (
            br3.innerShellContains(2. * h, x, y)
            ||
            br4.innerShellContains(2. * h, x, y)
    );
}

double fluxF(const double &x, const double &y)
{
    if (br3.innerShellContains(2. * h, x, y))
        return 0;
    if (br4.innerShellContains(2. * h, x, y))
        return 5*h;
//    else
    return 0;

//	if (bottomRect.contains(x, y)) return 50.;
//    else if (topRect.contains(x, y)) return 150;
}

double lambdaF(const double &x, const double &y)
{
//	if (std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25-5*h, 2.))
//		return 1e15;
//	else
    return 1.;
}

std::pair<double, double> sourceF(const double &x, const double &y, const double &T)
{
    double M = 1e1;
    if (br2.contains(x, y))
        return {M * T2, -M};
    else if (br4.outerShellContains(2. * h, x, y))
        return {M * T2 * 2, 0.};
    else if (br5.outerShellContains(2. * h, x, y) && x >= y)
        return {M * T2 / h, -M / h};

//    else
    return {0., 0.};
}


int main(int argc, char **argv)
{
    auto start = std::chrono::system_clock::now(); //start timer

    // Set geometry
    Domain domain;
    domain.addDomainFunction(domainFunction);

    // Initialize solver
    Solver solver(boundingRect, domain, h,
                  tempF, constTempFCond,
                  fluxF, constFluxFCond,
                  lambdaF, sourceF, T1);

    std::cout << "Run details:" << std::endl;
    std::cout << "h: " << solver.getH() << std::endl;
    std::cout << "[NX, NY]: [" << solver.getNx() << ", " << solver.getNy() << "]" << std::endl;
    std::cout << "Amount of cells: " << solver.getNy() * solver.getNx() << std::endl;

    // Start solving
    solver.solve(1e-2);

    std::fstream file;
    file.open("temp2.txt", std::ios::out);
    solver.exportTemp(file);
    file.close();

    auto end = std::chrono::system_clock::now();//end timer
    std::chrono::duration<double> elapsed = end - start;
    std::cout << std::endl << "Whole run took: " << elapsed.count() / 60 << "min" << std::endl;
    return 0;
}
