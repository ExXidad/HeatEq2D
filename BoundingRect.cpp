#include "BoundingRect.h"

BoundingRect::BoundingRect(const double &xmin, const double &xmax, const double &ymin, const double &ymax)
{
	if (xmax <= xmin || ymax <= ymin)
	{ throw std::invalid_argument("Error: wrong bounding rectangle dimensions"); }
	size = {{xmin, xmax},
			{ymin, ymax}};
}

const std::vector<std::vector<double>> &BoundingRect::getSize() const
{
	return size;
}

void BoundingRect::setSize(const std::vector<std::vector<double>> &size)
{
	BoundingRect::size = size;
}

double BoundingRect::getXSize()
{
	return size[0][1] - size[0][0];
}

double BoundingRect::getYSize()
{
	return size[1][1] - size[1][0];
}

bool BoundingRect::contains(const double &x, const double &y)
{
	return (size[0][0] <= x && x <= size[0][1]) && (size[1][0] <= y && y <= size[1][1]);
}
