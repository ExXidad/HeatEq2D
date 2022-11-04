#include <iostream>
#include <random>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"

constexpr double xmin = 0., xmax = 2.;
constexpr double ymin = 0., ymax = 1.;
constexpr double T0 = 1.;
BoundingRect bottomRect(xmin - 0.5, xmax + 0.5, -0.5 + ymin, ymin);
BoundingRect topRect(xmin - 0.5, xmax + 0.5, ymax, ymax + 0.5);

bool domainFunction(const double &x, const double &y)
{
	return true;
}

bool constTempFCond(const double &x, const double &y)
{
//	return !((xmin <= x && x <= xmax) && (ymin <= y && y <= ymax));
	return (bottomRect.contains(x, y)
			||
			topRect.contains(x, y));
}

double lambdaF(const double &x, const double &y)
{
	return 1.;
}

double sourceF(const double &x, const double &y)
{
	return BoundingRect(0.5, 1.5, 0.25, 0.75).contains(x, y);
	return 0.;
}

double tempF(const double &x, const double &y)
{
//	return y < (ymin + ymax) * 0.5 ? 50 : 150;
	if (bottomRect.contains(x, y))
		return 50.;
	else if (topRect.contains(x, y))
		return 150.;
}


int main(int argc, char **argv)
{
	auto start = std::chrono::system_clock::now(); //start timer

	// Set geometry
	BoundingRect boundingRect(xmin, xmax, ymin, ymax);
	Domain domain;
	domain.addDomainFunction(domainFunction);

	// Initialize solver
	Solver solver(boundingRect, domain, 1e-2,
				  tempF, constTempFCond, lambdaF, sourceF, 100.);

	std::cout << "Run details:" << std::endl;
	std::cout << "h: " << solver.getH() << std::endl;
	std::cout << "[NX, NY]: [" << solver.getNx() << ", " << solver.getNy() << "]" << std::endl;
	std::cout << "Amount of cells: " << solver.getNy() * solver.getNx() << std::endl;

	// Start solving
	solver.solve();

	std::fstream file;
	file.open("temp.txt", std::ios::out);
	solver.exportTemp(file, false);
	file.close();

	auto end = std::chrono::system_clock::now();//end timer
	std::chrono::duration<double> elapsed = end - start;
	std::cout << std::endl << "Whole run took: " << elapsed.count() / 60 << "min" << std::endl;
	return 0;
}
