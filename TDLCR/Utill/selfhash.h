#ifndef SELFHASH_H
#define SELFHASH_H
#include <vector>
#include <unordered_set>
#include<set>
#include <unordered_map>
#include<iostream>
#include<queue>
#include <functional>

using namespace std;
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& p) const {
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);
		return h1 ^ h2;
	}
};

struct pair_equal {
	template<class T1, class T2>
	bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
		return lhs.first == rhs.first && lhs.second == rhs.second;
	}
};


#endif
