#ifndef TDLCR_PL_H
#define TDLCR_PL_H

#include <vector>
#include <utility>//实用函数
#include <cmath>
#include<unordered_map>
#include <memory>//智能指针
#include <initializer_list>//用于表示初始化列表，可以方便地传递一组值给函数或构造函数
#include "Segment.h"
#include "Constants.h"
#include"selfhash.h"

using namespace std;

struct PLF {
	typedef shared_ptr<vector<Segment>> FPtr;//智能指针shared_ptr通过跟踪对象的引用计数，管理动态分配对象，避免内存泄漏
	typedef vector<Segment>::size_type size_type;

	FPtr f;
	unordered_map<pair<double, double>, int, pair_hash, pair_equal> labels;
	PLF() : f(make_shared<vector<Segment>>()) {}//make_shared用于创建shared_ptr对象

	PLF(initializer_list<Segment> segList) : f(make_shared<vector<Segment>>(segList)) {}
	/*用于初始化容器、构造函数和其他支持初始化列表*/
	PLF(const PLF &other) {
		if (other.f)
			f = make_shared<vector<Segment>>(*other.f);
		if (!other.labels.empty())
		{
			labels = other.labels;
		}

	}

	PLF(PLF &&rhs) noexcept : f(std::move(rhs.f)) {}//移动构造函数，将rhs内容移动到本对象中，无需创建副本

//copy assignments
	PLF &operator=(const PLF &other) {
		if (other.f)
			f = make_shared<vector<Segment>>(*other.f);
		if (!other.labels.empty())
		{
			labels = other.labels;
		}

		return *this;
	}

	PLF &operator=(PLF &&rhs) noexcept {
		f = rhs.f;
		return *this;
	}//移动赋值

	Segment &operator[](size_type id) {
		return (*f)[id];
	}
	//返回容器中id位置segment元素  
	friend ostream &operator<<(ostream &out, const PLF &plf) {
		for (auto &e : *plf.f)
			out << e << ",  ";
		return out;
	}

	//find segment for time t，找到出发时间t矩阵中对应点，返回元素迭代器
	vector<Segment>::iterator dpt2seg(double t) {
		auto beg = (*f).begin(), end = (*f).end();
		if (beg + 1 == end) { // 当前对象只有一个segment
			return beg;
		}
		auto mid = beg + (end - beg) / 2;
		while (mid != end) {// find the first t* > t, or end
			if (ge(t, mid->t)) // t* in [mid+1,end]
				beg = mid + 1;
			else    //t* in [begin, mid]
				end = mid;
			mid = beg + (end - beg) / 2;
		}
		return mid - 1;

	}
	vector<Segment>::iterator dpt2seg(double t, vector<Segment>& f) {
		auto beg = f.begin(), end = f.end();
		if (beg + 1 == f.end()) {
			return beg;
		}
		auto mid = beg + (end - beg) / 2;
		while (mid != end) {
			if (ge(t, mid->t))
				beg = mid + 1;
			else
				end = mid;
			mid = beg + (end - beg) / 2;
		}
		return mid - 1;
	}
	inline double dpt2wgt(double t, vector<Segment>::iterator s1) {  //weight func
		auto s2 = s1 + 1;
		if (s2 != f->end()) {// adopt s2 to get slop
			//assure s1->t≤t≤s2->t (right equal happen for minimize with same next point)
			assert(lt(s1->t, s2->t) && ge(t, s1->t) && le(t, s2->t));
			double k = (s2->w - s1->w) / (s2->t - s1->t);
			assert(gt(k, -1)); // assure FIFO property
			return k * (t - s1->t) + s1->w;
		}
		else { //last point
			assert(ge(t, s1->t));//&& le(t, TMAX)
			return s1->w;
		}
	}

	inline double dpt2wgt(double t) {
		vector<Segment>::iterator pseg = dpt2seg(t);
		return dpt2wgt(t, pseg);
	}

	// get dpt time for an arrival time
	inline double arr2dpt(double arr, std::vector<Segment>::iterator s1) {//根据到达时间，寻找出发时间
		auto s2 = s1 + 1;
		if (s2 != f->end()) {
			assert(lt(s1->t, s2->t) && ge(arr, s1->t + s1->w) && le(arr, s2->t + s2->w));  //assure s2 is not redundant
			double k = (s2->w - s1->w) / (s2->t - s1->t);
			assert(gt(k, -1)); // assure FIFO property
			double arr1 = s1->w + s1->t;
			return s1->t + (arr - arr1) / (k + 1);//逆插值运算找到出发点
		}
		else {
			assert(ge(arr, s1->t + s1->w));
			return arr - s1->w;
		}
	}

	// get arrival time from dpt time
	inline double dpt2arr(double t, std::vector<Segment>::iterator s1) {//根据出发时间计算到达时间
		auto s2 = s1 + 1;
		if (s2 != f->end()) {
			assert(lt(s1->t, s2->t) && ge(t, s1->t) && le(t, s2->t));  //assure s2 is not redundant
			double k = (s2->w - s1->w) / (s2->t - s1->t);

			return k * (t - s1->t) + s1->w + t;
		}
		else {    //last point
			assert(ge(t, s1->t));//&& le(t, TMAX)
			return s1->w + t;
		}
	}

	void constructpath(int path) {
		double tj, tj1;
		auto pw0 = f->begin();
		auto pw1 = pw0 + 1;
		while (pw1 != f->end())
		{
			tj = pw0->t;
			tj1 = pw1->t;
			labels[make_pair(tj, tj1)] = path;
			pw0++;
			pw1++;
		}
		tj = pw0->t;
		tj1 = TMAX;
		labels[make_pair(tj, tj1)] = path;
	}

	inline double dpt2arr(double td) {

		vector<Segment>::iterator pseg = dpt2seg(td);
		return dpt2arr(td, pseg);
	}

	// push a new segment point, initial  value: (ts,INT_MAX), intv=INTV_CNTED (directly connected),
	inline void push(double &pre_k, double dpt, double wgt, int intv) {
		auto &segs = *f;
		assert(le(dpt, TMAX) && le(wgt, INT_MAX) && gt(wgt, 0));
		if (segs.empty())
			segs.emplace_back(dpt, wgt, intv);
		else if (neq(segs.back().t, dpt)) {
			double new_k = (wgt - segs.back().w) / (dpt - segs.back().t);
			assert(gt(dpt, segs.back().t) && gt(new_k, -1)); // assure FIFO property
			// the next point wo pw0 may have different intv
			if (eq(pre_k, new_k)) {
				segs.back().t = dpt;
				segs.back().w = wgt;
				segs.back().intv = intv;
			}
			else {
				segs.emplace_back(dpt, wgt, intv);
				pre_k = new_k;
			}
		}
	}
	void compound(PLF &PLF0, PLF &PLFc, int intv) {    //suppose wgts belongs to (vs,ve),arr0 starts at v0
		double pre_k = INT_MAX;
		auto pw0 = PLF0.f->begin();
		if (pw0 == PLF0.f->end()) {
			PLFc.f = f;
			return;
		}
		auto pw3 = f->begin();
		if (pw3 == f->end()) {
			PLFc.f = PLF0.f;
			return;
		}
		double dpt = pw0->t; //dpt from start vertex v0
		double wgt0 = pw0->w;
		double arr0 = dpt + wgt0;
		auto pw1 = dpt2seg(arr0);// for edge (v_s, v_e), find dpt@v_s for dpt t at start v
		double wgt1 = dpt2wgt(arr0, pw1);
		PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
		while (true) {
			bool valid0 = pw0 + 1 < PLF0.f->end();
			bool valid1 = pw1 + 1 < f->end();
			int turn; // next point from base or from cur plf(1)
			if (valid0 && valid1) {
				double next_arr0 = (pw0 + 1)->t + (pw0 + 1)->w;
				if (lt((pw1 + 1)->t, next_arr0)) turn = 1;
				else turn = 0;
			}
			else if (valid1) turn = 1;
			else if (valid0) turn = 0;
			else break;
			if (turn) {
				pw1++;
				dpt = PLF0.arr2dpt(pw1->t, pw0);
			
				if (valid0) assert(lt(dpt, (pw0 + 1)->t));
				wgt0 = pw1->t - dpt;
				wgt1 = pw1->w;
			}
			else {
				pw0++;
				dpt = pw0->t;
				wgt0 = pw0->w;
				arr0 = dpt + wgt0;
				wgt1 = dpt2wgt(arr0, pw1);
			}
			PLFc.push(pre_k, dpt, wgt0 + wgt1, intv);
		}
	}

	double getMinW() {
		double minW = INT_MAX;

		for (const auto& segment : *f) {
			if (segment.w < minW) {
				minW = segment.w;
			}
		}

		return minW;
	}

	double getMaxW() {
		double maxW = INT_MIN;

		for (const auto& segment : *f) {
			if (segment.w > maxW) {
				maxW = segment.w;
			}
		}

		return maxW;
	}

	inline void minseg(double &pre_k, double t1, double t2, double f1, double f2,
		double g1, double g2, int intvf, int intvg) {
		assert(ge(t2, t1)); //t2≥t1
		if (eq(t1, t2)) return; //not interval
		if (eq(f1, g1)) {// f1=g1
			if (le(f2, g2)) push(pre_k, t1, f1, intvf);// f1=g1, f2≤g2
			else push(pre_k, t1, g1, intvg);// f1=g1, f2>g2
		}
		else { // f1!=g1
			double delta1 = g1 - f1;
			double delta2 = f2 - g2;//
			int flag_cross = 0;//
			if (gt(delta1, 0)) { // f1<g1
				push(pre_k, t1, f1, intvf);
				if (gt(delta2, 0)) flag_cross = -1; //f2>g2, crossed
			}
			else { //f1>g1
				push(pre_k, t1, g1, intvg);
				if (lt(delta2, 0)) flag_cross = 1; //f2<g2, crossed
			}
			if (flag_cross) {
				double denominator = delta1 + delta2;
				double inter_t = (delta1 * t2 + delta2 * t1) / denominator;
				double inter_w = (delta1 * g2 + delta2 * g1) / denominator;
				if (flag_cross == 1) push(pre_k, inter_t, inter_w, intvf);//f2<g2
				else push(pre_k, inter_t, inter_w, intvg);//f2>g2
			}
		}
	}


	//min{f,PLFc->f}
	void minimize(PLF &PLFc) {
		unordered_map<pair<double, double>, int, pair_hash, pair_equal> templabels = labels;
		labels.clear();
		double tj, tj1;
		double pre_k = INT_MAX;
		assert(!PLFc.f->empty());
		PLF PL_new;
		auto pw0 = f->begin(), pwc = PLFc.f->begin();
		assert(pw0->t == pwc->t); // same start time
		double t1 = pw0->t, t2, f1 = pw0->w, f2, g1 = pwc->w, g2;
		tj = t1;
		while (true) {
			int turn; //get next ptr to expand
			bool valid0 = pw0 + 1 < f->end();
			bool valid1 = pwc + 1 < PLFc.f->end();
			if (valid0 && valid1) {
				if (lt((pw0 + 1)->t, (pwc + 1)->t)) turn = 0;
				else turn = 1;
			}
			else if (valid0) turn = 0;
			else if (valid1) turn = 1;
			else break;
			if (turn) { // pwc turn
				t2 = (pwc + 1)->t;
				f2 = dpt2wgt(t2, pw0);
				g2 = (pwc + 1)->w;
			}
			else { //pw0 turn
				t2 = (pw0 + 1)->t;
				f2 = (pw0 + 1)->w;
				g2 = PLFc.dpt2wgt(t2, pwc);
			}
			PL_new.minseg(pre_k, t1, t2, f1, f2, g1, g2,
				pw0->intv, pwc->intv);
			//构造对应最小旅行时间函数路径
			assert(ge(t2, t1)); //t2≥t1

			double delta1 = g1 - f1;
			double delta2 = f2 - g2;//两个起始点之间的差异值
			int flag_cross = 0;//表示两个线段不相交
			if (gt(delta1, 0)) { // f1<g1

				if (gt(delta2, 0)) flag_cross = -1; //f2>g2, crossed
			}
			else { //f1>g1
				if (lt(delta2, 0)) flag_cross = 1; //f2<g2, crossed
			}
			if (flag_cross) {
				double denominator = delta1 + delta2;
				double inter_t = (delta1 * t2 + delta2 * t1) / denominator;
				tj1 = inter_t;
				if (flag_cross == 1) {//f2<g2
					for (const auto& pair : PLFc.labels) {
						if (inter_t > pair.first.first && inter_t <= pair.first.second) {
							labels[make_pair(tj, tj1)] = pair.second;
						}
					}
				}
				else { //f2>g2
					for (const auto& pair : templabels) {
						if (inter_t > pair.first.first && inter_t <= pair.first.second) {
							labels[make_pair(tj, tj1)] = pair.second;
						}
					}
				}
			}
			else {
				tj1 = t2;
				if (f1 <= g1) {
					for (const auto& pair : templabels) {
						if (tj1 > pair.first.first && tj1 <= pair.first.second) {
							labels[make_pair(tj, tj1)] = pair.second;
						}
					}
				}
				else {
					for (const auto& pair : PLFc.labels) {
						if (tj1 > pair.first.first && tj1 <= pair.first.second) {
							labels[make_pair(tj, tj1)] = pair.second;
						}
					}
				}
			}
			assert(lt(t2, TMAX));
			t1 = t2;
			tj = tj1;
			f1 = f2;
			g1 = g2;
			if (turn) pwc++; else pw0++;
		}
		PL_new.minseg(pre_k, t1, TMAX, f1, pw0->w, g1, pwc->w, // last point ,k = 0
			pw0->intv, pwc->intv);
		f = PL_new.f;


	}

};



#endif 
