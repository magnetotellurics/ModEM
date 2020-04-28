#include <string>
#include <vector>
#include <iostream>
//
//
#ifndef PERIOD_H_
#define PERIOD_H_
//
class Period {
public:
	int index;
	double period;
	std::vector< double > _pol_1_residuals;
	std::vector< double > _pol_2_residuals;
};
//
#endif

