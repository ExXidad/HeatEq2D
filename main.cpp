#include <iostream>
#include <random>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"
// smth to test branch
//DF - domain function
bool DF1(const double &x, const double &y)
{
	return y <= 2 * pow(10, -6);
}

bool DF2(const double &x, const double &y)
{
	return (3.33 <= x && x <= 6.66) && (1 <= y && y <= 2);
}

int main(int argc, char **argv)
{
	double filling_fraction = 0.35;

	if (argc != 3)
	{ throw std::runtime_error("Incorrect number of passed args"); }
	double n = atof(argv[1]);
	double U = atof(argv[2]);

	auto start = std::chrono::system_clock::now(); //start timer

	// Set geometry
	BoundingRect boundingRect(0, 0.025, 0, 0.016);
	Domain domain;
	domain.addDomainFunction(DF1);
//	domain.addDomainFunction(DF2);

	// Initialize solver
	Solver solver(boundingRect, domain, 9 * pow(10, -5), U);
	solver.setSaveProgressFlag(false);
	solver.randomSeed(n*filling_fraction);

	std::cout << "Run details:" << std::endl;
	std::cout << "h: " << solver.getH() << std::endl;
	std::cout << "[NX, NY]: [" << solver.getNx() << ", " << solver.getNy() << "]" << std::endl;
	std::cout << "Amount of cells: " << solver.getNy()*solver.getNx() << std::endl;
	std::cout << "Fraction: " << n << std::endl;
	std::cout << "U: " << solver.getU() << std::endl;
	std::cout << "D: " << solver.getD() << std::endl;
	std::cout << "mu: " << solver.getMu() << std::endl;


	// Start solving
	solver.solve(filling_fraction, 0.15);

	std::fstream file;
	file.open("domain.txt", std::ios::out);
	solver.exportComputationRegion(file);
	file.close();

	file.open("dendrite.txt", std::ios::out);
	solver.exportDendrite(file);
	file.close();

	file.open("data.txt", std::ios::out);
	solver.exportData(file);
	file.close();

	file.open("EP.txt", std::ios::out);
	solver.exportPotential(file);
	file.close();

	file.open("EF.txt", std::ios::out);
	solver.exportField(file);
	file.close();

	auto end = std::chrono::system_clock::now();//end timer
	std::chrono::duration<double> elapsed = end - start;
	std::cout << std::endl << "Whole run took: " << elapsed.count() / 60 << "min" << std::endl;
	return 0;
}
