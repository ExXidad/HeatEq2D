#include <iostream>
#include <random>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"

constexpr double h = 1e-2;
constexpr double xmin = 0., xmax = 2.;
constexpr double ymin = 0., ymax = 1.;
constexpr double T0 = 100.;
BoundingRect boundingRect(xmin, xmax, ymin, ymax);
BoundingRect lambdaRect(xmin, xmin + 0.5, ymax - 0.5, ymax);
BoundingRect bottomRect(xmin - 0.5, xmax + 0.5, -0.5 + ymin, ymin);
BoundingRect topRect(xmin - 0.5, xmax + 0.5, ymax, ymax + 0.5);

bool domainFunction(const double &x, const double &y)
{
	return true;
}

bool constTempFCond(const double &x, const double &y)
{
//	return !boundingRect.contains(x, y);
//			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25, 2.)
//			&&
//			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) >= std::pow(0.25 - 5 * h, 2.)
	return (
			bottomRect.contains(x, y)
			||
			topRect.contains(x, y)
			||
			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25, 2.)
			&&
			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) >= std::pow(0.25 - 5 * h, 2.)
	);
}

double lambdaF(const double &x, const double &y)
{
	if (std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25-5*h, 2.))
		return 1e15;
	else
		return 1.;
}

double sourceF(const double &x, const double &y, const double &T)
{
//	return BoundingRect(0.5, 1.5, 0.25, 0.75).contains(x, y);
	double M = 1e20;
	if (std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25 - 5 * h, 2.))
		return M * (T0 - T);
	else
		return 0.;
}

double tempF(const double &x, const double &y)
{
//	return y < (ymin + ymax) * 0.5 ? 50 : 150;
//	if (constTempFCond(x, y)) return T0;
	if (bottomRect.contains(x, y)) return 50.;
	else if (topRect.contains(x, y)) return 150;
	else if (
			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) <= std::pow(0.25, 2.)
			&&
			std::pow(x - 1.5, 2.) + std::pow(y - 0.5, 2.) >= std::pow(0.25 - 5 * h, 2.)
			)
		return T0;
}


int main(int argc, char **argv)
{
	auto start = std::chrono::system_clock::now(); //start timer

	// Set geometry
	Domain domain;
	domain.addDomainFunction(domainFunction);

	// Initialize solver
	Solver solver(boundingRect, domain, h,
				  tempF, constTempFCond, lambdaF, sourceF, T0);

	std::cout << "Run details:" << std::endl;
	std::cout << "h: " << solver.getH() << std::endl;
	std::cout << "[NX, NY]: [" << solver.getNx() << ", " << solver.getNy() << "]" << std::endl;
	std::cout << "Amount of cells: " << solver.getNy() * solver.getNx() << std::endl;

	// Start solving
	solver.solve(1e-4);

	std::fstream file;
	file.open("temp2.txt", std::ios::out);
	solver.exportTemp(file, false);
	file.close();

	auto end = std::chrono::system_clock::now();//end timer
	std::chrono::duration<double> elapsed = end - start;
	std::cout << std::endl << "Whole run took: " << elapsed.count() / 60 << "min" << std::endl;
	return 0;
}
