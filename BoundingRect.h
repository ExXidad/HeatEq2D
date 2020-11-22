#ifndef DENDRIVEV3_BOUNDINGRECT_H
#define DENDRIVEV3_BOUNDINGRECT_H

#include <vector>
#include <stdexcept>

class BoundingRect
{
private:
	std::vector<std::vector<double>> size = std::vector<std::vector<double>>(2, std::vector<double>(2, 0));
public:
	const std::vector<std::vector<double>> &getSize() const;

	void setSize(const std::vector<std::vector<double>> &size);

	double getXSize();

	double getYSize();


public:
	BoundingRect(const double &xmin, const double &xmax, const double &ymin, const double &ymax);

	bool contains(const double &x, const double &y);
};


#endif //DENDRIVEV3_BOUNDINGRECT_H
