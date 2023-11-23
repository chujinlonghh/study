#ifndef TDLCR_POINT_H
#define TDLCR_POINT_H

#include <climits>
#include <iostream>
#include "misc.h"
#include<limits.h>

struct Segment {
	double t;   //departure time
	double w;   //weight, i.e., travel time
	int intv;   // intermediate vertex of current segment or connenction info of intv belong to(negative)
	Segment() : t(0), w(INT_MAX), intv(INT_MIN) {}

	Segment(double _t, double _w) noexcept: t(_t), w(_w), intv(INT_MIN) {}

	Segment(double _t, double _w, int _intv) noexcept : t(_t), w(_w), intv(_intv) {}

	friend std::ostream &operator<<(std::ostream &os, const Segment &seg) {
		os << "(" << seg.t << ", " << seg.w << ", " << seg.intv << ")";
		return os;
	}

	inline bool operator==(const Segment& rhs) const {
		return eq(t, rhs.t) && eq(w, rhs.w) ;
	}

	inline bool operator!=(const Segment& rhs) const {
		return !(*this == rhs);
	}
};

#endif 
