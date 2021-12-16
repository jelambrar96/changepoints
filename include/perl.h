#ifndef PERL_H
#define PERL_H

#include <cmath>
#include <vector>

#include "cost.h"

std::vector<int> perl(AbstractCost * cost, double pen = NULL);

#endif
