#ifndef TDLCR_TDP2H_H
#define TDLCR_TDP2H_H
#include<iostream>
#include"TGraph.h"
#include <functional>
#include <string>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <queue>
#include <vector>
#include<algorithm>
#include <unordered_map>
#include"selfhash.h"

using namespace std;



struct ppl_LCR_pair {
	int node;
	int labels;
	function<double(double)> travetime;
	ppl_LCR_pair() {};
};
class LCR_index
{
public:
	unordered_map<int, vector<ppl_LCR_pair>> I_in;
	unordered_map<int, vector<ppl_LCR_pair>> I_out;
	unordered_set<int> prunrule1;
	LCR_index() {};
	unordered_map<pair<int, int>, vector<PLF>, pair_hash> preplfs;
	std::vector<int> getvertexorder(TGraph& TG)
	{
		std::vector<int> vertexorder;

		std::vector<int> outCount(TG.n);
		std::vector<int> countvex(TG.n);  

		for (int i = 0; i < TG.n; i++) {
			double count = 0.0;
			countvex[i] = 0;  // 初始化 countvex
			int timevex = 0;
			for (int j = TG.head[i]; j != -1; j = TG.next[j]) {
				countvex[i]++;  // 更新 countvex
				timevex += (TG.turnpoint[j] - 1);
			}
			for (int j = TG.rhead[i]; j != -1; j = TG.rnext[j]) {
				countvex[i]++;  // 更新 countvex
			}
			count = timevex;
			outCount[i] = count;
		}

		vertexorder.clear();
		for (int i = 0; i < TG.n; i++) {
			vertexorder.push_back(i);
		}

		// 修改排序条件
		std::sort(vertexorder.begin(), vertexorder.end(), [&](int a, int b) {
			if (outCount[a] == outCount[b]) {
				// 如果 timevex 相同，按照 countvex 进行比较
				if (countvex[a] == countvex[b]) {
					// 如果 countvex 也相同，按照顶点序号降序排列
					return a > b;
				}
				return countvex[a] > countvex[b];
			}
			return outCount[a] > outCount[b];
		});
		return vertexorder;
	}

	void insert_I_in(int index_node, int node, int label, function<double(double)> travetime)
	{
		ppl_LCR_pair temp;
		temp.node = node;
		temp.labels = label;
		temp.travetime = travetime;
		I_in[index_node].push_back(temp);
	}
	void insert_I_out(int index_node, int node, int label, function<double(double)> travetime)
	{
		ppl_LCR_pair temp;
		temp.node = node;
		temp.labels = label;
		temp.travetime = travetime;
		I_out[index_node].push_back(temp);
	}

	void construction(TGraph& TG) {
		unordered_set<int> visited;
		unordered_set<int> pathSet;
		vector<int> vertexorder(TG.n);
		vertexorder = getvertexorder(TG);
		int i = 0;
		for (auto it : vertexorder) {
			PLF prePLF = PLF();
			dfs(TG, it, it, visited, pathSet, 0, prePLF);
			dfsReverse(TG, it, it, visited, 0, prePLF);
			for (auto it1 = preplfs.begin(); it1 != preplfs.end();) {
				if (it1->first.second == it || it1->first.first == it) {
					it1 = preplfs.erase(it1);
				}
				else {
					++it1;
				}
			}
			prunrule1.insert(it);
			pathSet.clear();
		}
	}
	int reverseInteger(int num) {
		int reversedNum = 0;

		while (num > 0) {
			int remainder = num % 10;
			reversedNum = reversedNum * 10 + remainder;
			num = num / 10;
		}

		return reversedNum;
	}
	bool isDigitInInteger(int digit, int num) {
		while (num > 0) {
			int lastDigit = num % 10;
			if (lastDigit == digit) {
				return true;
			}
			num = num / 10;
		}
		return false;
	}
	// 递归遍历函数
	void dfs( TGraph& g, int src, int v, unordered_set<int>& visited, unordered_set<int>& pathSet, int path, PLF &prePLF)
	{
		visited.insert(v);
		int savedPath = path; // 保存当前的path值
		for (int j = g.head[v]; j != -1; j = g.next[j]) {
			int u = g.adjv[j];
			//当前访问点已经在以前遍历过
			PLF cur = PLF();
			if (prunrule1.find(u) != prunrule1.end() || u == src) {

				continue;
			}
			g.weights[j].compound(prePLF, cur, v);
			int label = g.label[j];
			if (isDigitInInteger(label, path / 10))
			{
				continue;
			}
			auto travetime = [cur](double td) mutable-> double {
				return cur.dpt2arr(td);
			};
			PLF curt = cur;
			if (!preplfs[make_pair(src, u)].empty()) {
				PLF plf_ = preplfs[make_pair(src, u)].front();
				curt.minimize(plf_);
				preplfs[make_pair(src, u)].pop_back();
				preplfs[make_pair(src, u)].push_back(curt);
			}
			else {
				preplfs[make_pair(src, u)].emplace_back(curt);
			}

			pathSet.insert(label);
			if (label != path % 10) {
				path = path * 10 + label;
			}
			insert_I_in(u, src, path, travetime);
			if (visited.find(u) == visited.end()) {
				dfs(g, src, u, visited, pathSet, path, curt);
			}
			pathSet.erase(label);
			path = savedPath; // 恢复path的值
		}
		visited.erase(v);
	}
	void dfsReverse( TGraph& g, int src, int v, unordered_set<int>& visited, int path, PLF &prePLF) {
		visited.insert(v);
		int savedPath = path; // 保存当前的path值
		for (int j = g.rhead[v]; j != -1; j = g.rnext[j]) {
			PLF cur = PLF();
			int u = g.radjv[j];
			if (prunrule1.find(u) != prunrule1.end() || u == src) {
				continue;
			}

			prePLF.compound(g.weights[j], cur, v);

			if (visited.find(u) == visited.end()) {

				int label = g.label[j];
				if (isDigitInInteger(label, path / 10))
				{
					continue;
				}
				auto travetime = [cur](double td)mutable -> double {

					return cur.dpt2arr(td);
				};
				PLF curt = cur;
				if (!preplfs[make_pair(u, src)].empty()) {
					PLF plf_ = preplfs[make_pair(u, src)].front();
					curt.minimize(plf_);
					preplfs[make_pair(u, src)].pop_back();
					preplfs[make_pair(u, src)].push_back(curt);
				}
				else {
					preplfs[make_pair(u, src)].emplace_back(curt);
				}
				if (label != path % 10) {
					path = path * 10 + label;
				}
				insert_I_out(u, src, reverseInteger(path), travetime);

				dfsReverse(g, src, u, visited, path, curt);
				
			}
			path = savedPath; // 恢复path的值
		}
		visited.erase(v);
	}
};



#endif