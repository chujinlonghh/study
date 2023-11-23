#ifndef TDLCR_CONSTANTS_H
#define TDLCR_CONSTANTS_H

#include <string>
#include <map>
#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <cassert>


#define DE_INTV -2147483648//32位机器能表示最小有符号整数
#define INTV_CNTED -2147483647
#define DE_W 2147483647//32位机器能表示最大有符号整数


double EPSILON = 1e-5; 
int TMAX = 86400;
unsigned long LEAF_SIZE = 64;
unsigned long FANOUT = 4;


#endif 
