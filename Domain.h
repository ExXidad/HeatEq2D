//
// Created by xidad on 26.10.2020.
//

#ifndef DENDRIVEV3_DOMAIN_H
#define DENDRIVEV3_DOMAIN_H

#include <vector>

class Domain
{
private:
	std::vector<bool (*)(const double &, const double &)> domainFunctions;

public:
	Domain();

	void addDomainFunction(bool(&domainFunction)(const double &, const double &));

	bool intersectionContains(const double &x, const double &y);

	bool unionContains(const double &x, const double &y);
};


#endif //DENDRIVEV3_DOMAIN_H
