#ifndef TDP_TWOQUERY_H
#define TDP_TWOQUERY_H

#include<iostream>
#include <functional>
#include <string>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <queue>
#include <vector>
#include<algorithm>
#include <unordered_map>
#include"tdp2h.h"
#include"pathassit.h"

using namespace std;

bool tdptwohquery(int s, int d, double td, double y, LCR_index &lcrindex, int language)
{
	int index = 0;
	vector<int> path = convert(language);
	//检查s是否在Iin(d)中
	if (lcrindex.I_in[d].size() > 0) {
		for (const auto& pair : lcrindex.I_in[d]) {
			int currentIndex = index; // 保存当前的 index 值
			if (pair.node == s) {
				double arr = pair.travetime(td);
				if (arr <= y && checkispath(path, index, pair.labels)) {
						std::cout << " 到达喽！" << endl;
						return true;
				}
			}
			index = currentIndex;
		}
	}
	// 检查d是否在Iout(s)中
	if (lcrindex.I_out[s].size() > 0) {
		for (const auto& pair : lcrindex.I_out[s]) {
			int currentIndex = index; // 保存当前的 index 值
			if (pair.node == d) {
				double arr = pair.travetime(td);
				if (arr <=y && checkispath(path, index, pair.labels)) {
						std::cout  << " 到达喽！" << endl;
						return true;
				}
			}
			index = currentIndex;
		}
	}
	// 检查Iout(s)和Iin(d)的交集
	if (lcrindex.I_out[s].size() > 0 && lcrindex.I_in[d].size() > 0) {
		unordered_set<int> commonNodes;
		// 找到存在于Iout(s)和Iin(d)中的公共顶点
		for (const auto& outPair : lcrindex.I_out[s]) {
			for (const auto& inPair : lcrindex.I_in[d]) {
				if (outPair.node == inPair.node) {
					commonNodes.insert(outPair.node);
				}
			}
		}
		// 检查公共顶点的可达性
		for (int node : commonNodes) {
			double sToNodeTravelTime = 0.0;
			double TravelTime = 0.0;
			// 获取s到公共顶点的行程时间
			for (const auto& outPair : lcrindex.I_out[s]) {
				int currentIndex = index; // 保存当前的 index 值
				if (outPair.node == node && checkispath(path, index, outPair.labels)) {
					sToNodeTravelTime = outPair.travetime(td);
					for (const auto& inPair : lcrindex.I_in[d]) {
						int curIndex = index; // 保存当前的 index 值
						if (inPair.node == node && checkispath(path, index, inPair.labels)) {
							TravelTime = inPair.travetime(sToNodeTravelTime);
							if (TravelTime <= y) {
								std::cout << " 到达喽！" << endl;
								return true;
							}
						}
						index = curIndex;
					}
				}
				index = currentIndex;
			}	
		}
	}
	std::cout << "无法到达";
	return false;
}












#endif