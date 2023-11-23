#ifndef TDLCR_UTILS_H
#define TDLCR_UTILS_H

#include <iostream>
#include <stack>
#include <vector>
#include "Constants.h"


inline bool le(const double &x, const double &y) { return x <= y + EPSILON; }

inline bool lt(const double &x, const double &y) { return x + EPSILON < y; }

inline bool eq(const double &x, const double &y) { return fabs(x - y) <= EPSILON; }

inline bool neq(const double &x, const double &y) { return !eq(x, y); }

inline bool gt(const double &x, const double &y) { return lt(y, x); }

inline bool ge(const double &x, const double &y) { return le(y, x); }




#endif 
