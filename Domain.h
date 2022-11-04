#ifndef DENDRIVEV3_DOMAIN_H
#define DENDRIVEV3_DOMAIN_H

#include <vector>

enum class DFInteractionType
{
	UNION, INTERSECTION, EXCLUSION
};

class Domain
{
private:
	std::vector<bool (*)(const double &, const double &)> domainFunctions;
	DFInteractionType type = DFInteractionType::UNION;

	bool (Domain::*containsFunction)(const double &, const double &) = &Domain::intersectionContains;

public:
	Domain();

	void setDFInteractionType(const DFInteractionType &type);

	bool contains(const double &x, const double &y);

	void addDomainFunction(bool(&domainFunction)(const double &, const double &));

	bool intersectionContains(const double &x, const double &y);

	bool unionContains(const double &x, const double &y);
};


#endif //DENDRIVEV3_DOMAIN_H
