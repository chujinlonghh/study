#ifndef TDLCR_ONEINDEX_H
#define TDLCR_ONEINDEX_H
#include <vector>
#include <unordered_set>
#include<set>
#include <unordered_map>
#include<iostream>
#include<queue>
#include <functional>
#include"TGraph.h"
#include "PLF.h"
#include "Constants.h"
#include"selfhash.h"

using namespace std;

class Index
{
private:
	unordered_map<pair<int, int>, queue<pair<int, function<double(double)>>>, pair_hash, pair_equal> index;
	int vertexnum;
public:
	unordered_map<pair<int, int>, vector<PLF>, pair_hash> preplfs;
	Index() {};
	void buildindex(TGraph& TG)
	{
		this->vertexnum = TG.n;
		for (int i = 0; i < TG.n; i++)
		{
			for (int j = TG.head[i]; j != -1; j = TG.next[j]) {
				int nei = TG.adjv[j];
				int label = TG.label[j];
				PLF prePLF = TG.weights[j];
				auto travetime_ = [j, &TG](double td) -> double {
					return TG.weights[j].dpt2arr(td);
				};
				index[make_pair(i, label)].push(make_pair(nei, travetime_));
				buildsame(TG, i, label, nei, travetime_, prePLF);
			}
			for (auto it = preplfs.begin(); it != preplfs.end();) {
				if (it->first.first == i) {
					it = preplfs.erase(it);
				}
				else {
					++it;
				}
			}
		}
	}
	void buildsame(TGraph& TG, int src, int prelabel, int cur, function<double(double)> traveltime, PLF &prePLF) {
		for (int j = TG.head[cur]; j != -1; j = TG.next[j]) {
			PLF curt = PLF();
			int nei = TG.adjv[j];
			int label = TG.label[j];
			TG.weights[j].compound(prePLF, curt, cur);
		
			if (!preplfs[make_pair(src, cur)].empty()) {
				PLF plf_ = preplfs[make_pair(src, cur)].front();
				curt.minimize(plf_);
				preplfs[make_pair(src, cur)].pop_back();
				preplfs[make_pair(src, cur)].push_back(curt);
			}
			else {
				preplfs[make_pair(src, cur)].emplace_back(curt);
			}
			auto travetime = [curt](double td) mutable -> double {
				return curt.dpt2arr(td);
			};
			if (label == prelabel) {
				auto key = make_pair(src, label);

				if (index.find(key) != index.end())
				{
					bool found = false;
					auto key = make_pair(src, label);
					if (index.find(key) != index.end())
					{
						// 判断队列中是否已存在相同的元素
						auto queue = index[key];
						while (!queue.empty())
						{
							auto elem = queue.front();
							queue.pop();
							if (elem.first == nei)
							{
								found = true;
								break;
							}
						}
						if (!found) {
							index[make_pair(src, prelabel)].push(make_pair(nei, travetime));
							buildsame(TG, src, label, nei, travetime, curt);
						}
					}
				}
			}
		}
	}
	vector<pair<int, function<double(double)>>> getData(int vertex, int label)
	{
		vector<pair<int, function<double(double)>>> result;

		auto it = index.find(make_pair(vertex, label));
		if (it != index.end())
		{
			queue<pair<int, function<double(double)>>> tempQueue; // 临时队列
			while (!it->second.empty())
			{
				result.push_back(it->second.front());
				tempQueue.push(it->second.front());
				it->second.pop();
			}
			while (!tempQueue.empty())
			{
				it->second.push(tempQueue.front());
				tempQueue.pop();
			}
		}

		return result;
	}
	int getvertexnum()
	{
		return vertexnum;
	}
	//根据顶点获取外向边标签
	vector<int> getlabelFromVertex(int vertex)
	{
		vector<int> result;
		set<int> labelSet;

		for (const auto& entry : index)
		{
			if (entry.first.first == vertex && labelSet.find(entry.first.second) == labelSet.end())
			{
				result.push_back(entry.first.second);
				labelSet.insert(entry.first.second);
			}
		}

		return result;
	}
};


#endif
