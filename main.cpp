#include <iostream>
#include <random>
#include <cmath>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"

//DF - domain function
bool DF1(const double &x, const double &y)
{
	return y <= 1;
}

bool DF2(const double &x, const double &y)
{
	return (3.33 <= x && x <= 6.66) && (1 <= y && y <= 2);
}

//ICF - initial condition function
double ICF1(const double &x, const double &y)
{
	return exp(-((x - 5) * (x - 5) + (y - 5) * (y - 5)));
}

int main()
{
	BoundingRect boundingRect(0, 10, 0, 10);

	Domain domain;
	domain.addDomainFunction(DF1);
	domain.addDomainFunction(DF2);
	domain.setDFInteractionType(UNION);


	Solver solver(boundingRect, domain, 0.05, 10);
	solver.solve(ICF1, MINMOD, 0.01, 1);
	return 0;
}
