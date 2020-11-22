#include "Domain.h"

void Domain::addDomainFunction(bool (&domainFunction)(const double &, const double &))
{
	domainFunctions.emplace_back(&domainFunction);
}

bool Domain::intersectionContains(const double &x, const double &y)
{
	for (auto &domainFunction : domainFunctions) {
		if (!domainFunction(x, y)) return false;
	}
	return true;
}

bool Domain::unionContains(const double &x, const double &y)
{
	for (auto &domainFunction : domainFunctions) {
		if (domainFunction(x, y)) return true;
	}
	return false;
}

Domain::Domain()
{

}

void Domain::setDFInteractionType(const DFInteractionType &type)
{
	this->type = type;
	switch (type) {
		case UNION:
			containsFunction = &Domain::unionContains;
			break;
		case INTERSECTION:
			containsFunction = &Domain::intersectionContains;
			break;
	}
}

bool Domain::contains(const double &x, const double &y)
{
	return (this->*containsFunction)(x, y);
}
