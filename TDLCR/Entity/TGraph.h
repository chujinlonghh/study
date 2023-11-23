#ifndef TDLCR_TGRAPH_H
#define TDLCR_TGRAPH_H

#include <memory>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include "PLF.h"
#include "Segment.h"
#include "misc.h"

struct TGraph {
	unsigned long n, m;
	int mCnt; //count current added edges
	std::vector<int> id, head, next, adjv, turnpoint;  //id maintain global id of vertices
	std::vector<int> rhead, rnext, radjv;  // reverse adjacency list
	std::vector<PLF> weights;   //time dependent edge weights
	std::unordered_map<int, int> label;

	TGraph() {
		n = 0;
		m = 0;
		mCnt = 0;

	}

	~TGraph() {
		head.clear();
		id.clear();
		next.clear();
		adjv.clear();
		rhead.clear();
		rnext.clear();
		radjv.clear();
		weights.clear();
		turnpoint.clear();
	}

	void init(unsigned long _n, unsigned long _m) {
		this->n = _n;
		this->m = _m;
		head = std::vector<int>(n, -1);
		id = std::vector<int>(n, -1);
		next = std::vector<int>(m, -1);
		adjv = std::vector<int>(m, -1);
		rhead = std::vector<int>(n, -1);
		rnext = std::vector<int>(m, -1);
		radjv = std::vector<int>(m, -1);
		turnpoint = std::vector<int>(m, -1);
		weights = std::vector<PLF>(m);
	}

	void readGraph(const std::string& path) {
		std::ifstream in(path);
		if (!in.is_open()) {
			std::cout << "Error opening file: " << path << std::endl;
			return;
		}
		assert(in.is_open());
		unsigned int n, m;
		in >> n >> m;
		init(n, m);
		int vs, vt, weight_piece_num;
		int lab;
		while (in >> vs >> vt >> weight_piece_num) { //input edges
			auto f = std::make_shared<std::vector<Segment>>(weight_piece_num);
			double t, w;
			for (int i = 0; i < weight_piece_num; i++) {
				in >> t >> w;
				(*f)[i] = { t, w, INTV_CNTED };
			}
			in >> lab;
			addEdge(vs, vt, f, lab, weight_piece_num);
		}
	}

	void addEdge(int s, int t, std::shared_ptr<std::vector<Segment>>& f, int lab, int num) {
		adjv[mCnt] = t;
		next[mCnt] = head[s];
		head[s] = mCnt;
		weights[mCnt].f = f;
		label[mCnt] = lab;
		turnpoint[mCnt] = num;

		radjv[mCnt] = s;
		rnext[mCnt] = rhead[t];
		rhead[t] = mCnt;

		mCnt++;
	}

};

#endif