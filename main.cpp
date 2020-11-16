#include <iostream>
#include <random>
#include <cmath>

#include "Domain.h"
#include "BoundingRect.h"
#include "Solver.h"

//ICF - initial condition function
double ICF1(const double &x, const double &y)
{
	if (-3 <= x && x <= 3) return 1;
	return 0;
}

double ICF2(const double &x, const double &y)
{
	if (-3 <= x && x <= 3) return pow(2.71,-x*x/2)-pow(2.71,-9./2);
	return 0;
}

int main()
{
	// Bounding rectangle x-dimensions [xmin,xmax], y-dimensions [0,tMax]
	BoundingRect boundingRect(-10, 10, 0, 10);

	Solver solver(boundingRect, 0.05, 1, 0.1, 0.1);
	solver.solve(ICF2, MINMOD);

	solver.save("data");
	return 0;
}
