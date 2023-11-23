#ifndef TDLCR_BII_H
#define TDLCR_BII_H

#include <iostream>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include"TGraph.h"


using namespace std;

struct biibn_pair {
	int node;
	int labels;
	double weight;
	biibn_pair() {};
};

struct biibf_pair {
	int node;
	double weight;
	function<double(double)> mintravetime;
	unordered_map<pair<double, double>, int, pair_hash, pair_equal> labels;
	biibf_pair() {};
}; 

struct BNin_out {
	vector<biibn_pair> in;
	vector<biibn_pair> out;
};

struct BFin_out {
	vector<biibf_pair> in;
	vector<biibf_pair> out;
};

class BII_index {
public:
	unordered_map<int, BNin_out> BN;
	unordered_map<int, BFin_out> BF;
	unordered_set<int> visited;
	unordered_map<pair<int, int>, vector<PLF>, pair_hash> preplfs;

	void biiconstruct(TGraph&TG, double delta)
	{
		buildBN(TG, delta);
		buildBF(TG, delta);
	}

	void insertBN(int src,int dst,int path,double weight) {
		biibn_pair outpair;
		outpair.node = dst;
		outpair.labels = path;
		outpair.weight = weight;
		BN[src].out.push_back(outpair);

		biibn_pair inpair;
		inpair.node = src;
		inpair.labels = path;
		inpair.weight = weight;
		BN[dst].in.push_back(inpair);
	}

	void insertBF(int src, int dst, double weight, function<double(double)> mintravetime, unordered_map<pair<double, double>, int, pair_hash, pair_equal> labels) {
		biibf_pair outpair;
		outpair.node = dst;
		outpair.mintravetime = mintravetime;
		outpair.weight = weight;
		outpair.labels = labels;
		BF[src].out.push_back(outpair);

		biibf_pair inpair;
		inpair.node =src;
		inpair.mintravetime = mintravetime;
		inpair.weight = weight;
		inpair.labels = labels;
		BF[dst].in.push_back(inpair);
	}

	void buildBN(TGraph& G, double delta) {
		for (int v = 0; v < G.n; v++) {
			unordered_set<int> visited;
			extendBN(G, v,v, delta,0.0, visited,0);
		}
	}

	void extendBN(TGraph& G, int src,int dst, double delta,double travelt, unordered_set<int>& visited, int path) {
		if (G.head[dst] == -1) {
			return;
		}
		visited.insert(dst);
		int savedPath = path; // 保存当前的path值
		
		for (int j = G.head[dst]; j != -1; j = G.next[j]) {
			int u = G.adjv[j];
			if (src == u) {
				continue;
			}
			int lab = G.label[j];
			if (isDigitInInteger(lab, path / 10))
			{
				continue;
			}
			double travelTime = G.weights[j].getMaxW()+travelt;
			if (lab != path % 10) {
				path = path * 10 + lab;
			}
			if (travelTime <=delta) {
				insertBN(src, u, path, travelTime);
				if (visited.find(u) == visited.end()) {
					extendBN(G, src, u, delta, travelTime, visited, path);
				}
				path = savedPath;
			}
		}
		visited.erase(dst);
	}
	void buildBF(TGraph& G, double delta) {
		for (int v = 0; v < G.n; v++) {
			PLF prePLF = PLF();
			extendBF(G, v,v, delta,0.0, visited, 0,prePLF);
			for (auto it = preplfs.begin(); it != preplfs.end();) {
				if (it->first.first == v) {
					it = preplfs.erase(it);
				}
				else {
					++it;
				}
			}
			//对preves清空；
		}
	}
	
	void extendBF(TGraph& g, int src, int v, double delta, double travelt, unordered_set<int>& visited, int path, PLF &prePLF)
	{
		if (g.head[v] == -1) {
			return;
		}
		visited.insert(v);
		int savedPath = path; // 保存当前的path值
		for (int j = g.head[v]; j != -1; j = g.next[j]){
			int u = g.adjv[j];
			PLF cur=PLF();
			if (src == u) {
				continue;
			}
			int label = g.label[j];
			if (isDigitInInteger(label, path / 10))
			{
				continue;
			}
			if (label != path % 10) {
				path = path * 10 + label;
			}
			g.weights[j].compound(prePLF,cur, v);
			assert(cur.f->front().w != DE_W);
			cur.constructpath(path);
			if (!preplfs[make_pair(src, u)].empty()) {
				PLF plf_ = preplfs[make_pair(src, u)].front();
				cur.minimize(plf_);
				preplfs[make_pair(src, u)].pop_back();
				preplfs[make_pair(src, u)].push_back(cur);
			}
			else {
 				preplfs[make_pair(src, u)].emplace_back(cur);
			}
			auto travetime = [cur](double td) mutable -> double {
				return cur.dpt2arr(td);
			};
			double statictravelTime = g.weights[j].getMinW() + travelt;
			if (statictravelTime <= delta) {
				bool existed = true;
				for (auto &it : BF[src].out)
				{
					if (it.node == u)
					{
						if (it.weight > statictravelTime)
						{
							it.weight = statictravelTime;
						}
						it.mintravetime = travetime;
						it.labels = cur.labels;
						for (auto &it : BF[u].in)
						{
							if (it.node == src)
							{
								if (it.weight > statictravelTime)
								{
									it.weight = statictravelTime;
								}
								it.mintravetime = travetime;
								it.labels = cur.labels;
								break;
							}
						}
						existed = false;
						break;
					}
				}
				if (existed) {
					insertBF(src, u, statictravelTime, travetime,cur.labels);
				}
				if(visited.find(u) == visited.end()) {
					extendBF(g, src, u, delta, statictravelTime, visited, path, cur);
				}
				path = savedPath;
			}
		}
		visited.erase(v);
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
	
};

#endif