#ifndef BII_QUERY_H
#define BII_QUERY_H

#include"biiindex.h"
#include"pathassit.h"
#include<cmath>
struct biiqueryf {
	int node;
	double weight;
	function<double(double)> mintravetime;
	function<int(double)> labels;
	biiqueryf(){};
};

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
int reverseInteger(int num) {
	int reversedNum = 0;

	while (num > 0) {
		int remainder = num % 10;
		reversedNum = reversedNum * 10 + remainder;
		num = num / 10;
	}

	return reversedNum;
}
bool isFirstDigitEqualTo(int num, int digit) {
	// 获取首位数
	int firstDigit = num;
	while (firstDigit >= 10) {
		firstDigit /= 10;
	}

	// 判断首位数与给定数字是否相等
	return (firstDigit == digit);
}
vector<biibn_pair>FISN(int s, int k, BII_index& biindex, double td, int path, double y) {
	// 如果已经达到最大深度，返回空集合
	if (k == 0) {
		return vector<biibn_pair>();
	}
	// 在索引中查找起点 s 在 delta 时间内可以到达的点的集合
	vector<biibn_pair> result;
	for (auto it : biindex.BN[s].out)
	{
		double curtime = it.weight + td;
		if (isDigitInInteger(it.labels, path / 10))
		{
			continue;
		}
		int curpath = path;
		if (it.labels != path % 10) {
			curpath =path * 10 + it.labels;
		}
		if (curtime < y)
		{
			biibn_pair temp;
			temp.labels = curpath;
			temp.node = it.node;
			temp.weight = curtime;
			result.push_back(temp);
			// 递归调用，对每条路径都进行递归调用，并将结果合并到最终的结果集中
			vector<biibn_pair> subResult = FISN(it.node, k - 1, biindex, curtime, curpath, y);
			result.insert(result.end(), subResult.begin(), subResult.end());
		}
	}
	return result;
}
vector<biibn_pair>BISN(int d, int k, BII_index& biindex, double td, int path, double y) {
	// 如果已经达到最大深度，返回空集合
	if (k == 0) {
		return vector<biibn_pair>();
	}
	// 在索引中查找起点 s 在 delta 时间内可以到达的点的集合
	vector<biibn_pair> result;
	for (auto it : biindex.BN[d].in)
	{
		double curtime = it.weight + td;
		if (isDigitInInteger(it.labels, path / 10))
		{
			continue;
		}
		int curpath = path;
		if (it.labels != path % 10) {
			curpath = path * 10 + it.labels;
		}
		curpath = reverseInteger(curpath);
		if (curtime < y)
		{
			biibn_pair temp;
			temp.labels = curpath;
			temp.node = it.node;
			temp.weight = curtime;
			result.push_back(temp);
			// 递归调用，对每条路径都进行递归调用，并将结果合并到最终的结果集中
			vector<biibn_pair> subResult = FISN(it.node, k - 1, biindex, curtime, curpath, y);
			result.insert(result.end(), subResult.begin(), subResult.end());
		}
	}
	return result;
}

vector<biiqueryf> FISF(int s, int k, BII_index& biindex, double td, double y, function<double(double)> pretraveltime, function<int(double)> preminpath) {
	// 如果已经达到最大深度，返回空集合
	if (k == 0) {
		return vector<biiqueryf>();
	}
	// 在索引中查找起点 s 在 delta 时间内可以到达的点的集合
	vector<biiqueryf> result;
	for (auto it : biindex.BF[s].out)
	{
		double curtime = td + it.weight;
		if (curtime < y)
		{
			biiqueryf temp;
			temp.node = it.node;
			temp.weight = curtime;
			auto  travetime = [it, pretraveltime](double t)->double {
				double result = pretraveltime(t);
				return it.mintravetime(result);
			};
			auto  minpath = [it, preminpath, pretraveltime](double t)->int {
				int curp;
				bool findpath = false;
				double td = pretraveltime(t);
				for (const auto &curit : it.labels)
				{
					if (td > curit.first.first && td <= curit.first.second) {
						curp = curit.second;
						findpath = true;
						break;
					}
				}
				if (!findpath)
				{
					pair<double, double> maxkey{ 0.0, 0.0 };
					for (const auto& entry : it.labels) {
						const pair<double, double>& key = entry.first;
						const int& value = entry.second;
						if (key.first > maxkey.first || (key.first == maxkey.first && key.second > maxkey.second)) {
							maxkey = key;
							curp = value;
						}
					}
				}
				int s = preminpath(t);
				if (curp != s % 10) {
					s = s * 10 + curp;
				}
				return s;
			};
			temp.mintravetime = travetime;
			temp.labels = minpath;
			result.push_back(temp);
			vector<biiqueryf> subResult = FISF(it.node, k - 1, biindex, curtime, y, travetime, minpath);
			result.insert(result.end(), subResult.begin(), subResult.end());
		}
	}
	return result;
}
vector<biiqueryf> BISF(int d, int k, BII_index& biindex, double td, double y, function<double(double)> pretraveltime, function<int(double)> preminpath) {
	// 如果已经达到最大深度，返回空集合
	if (k == 0) {
		return vector<biiqueryf>();
	}
	// 在索引中查找起点 s 在 delta 时间内可以到达的点的集合
	vector<biiqueryf> result;
	for (auto it : biindex.BF[d].in)
	{
		double curtime = td + it.weight;
		if (curtime < y)
		{
			biiqueryf temp;
			temp.node = it.node;
			temp.weight = curtime;
			auto  travetime = [it, pretraveltime](double t)->double {
				double result = it.mintravetime(t);
				double s = pretraveltime(result);
				return pretraveltime(result);
			};
			auto  minpath = [it, preminpath](double t)->int {
				int curp;
				bool findpath = false;
				double result = it.mintravetime(t);
				for (const auto &curit : it.labels)
				{
					if (t > curit.first.first && t <= curit.first.second) {
						curp = curit.second;
						findpath = true;
						break;
					}
				}
				if (!findpath)
				{
					pair<double, double> maxkey{ 0.0, 0.0 };
					for (const auto& entry : it.labels) {
						const pair<double, double>& key = entry.first;
						const int& value = entry.second;
						if (key.first > maxkey.first || (key.first == maxkey.first && key.second > maxkey.second)) {
							maxkey = key;
							curp = value;
						}
					}
				}
				int s = preminpath(result);
				if (!isFirstDigitEqualTo(s,curp)&&s!=0) {
					s = curp * 10 +s;
				}
				else if (s == 0) {
					s = curp;
				}
				return s;
			};
			temp.mintravetime = travetime;
			temp.labels = minpath;
			result.push_back(temp);
			vector<biiqueryf> subResult = BISF(it.node, k - 1, biindex, curtime, y, travetime, minpath);
			result.insert(result.end(), subResult.begin(), subResult.end());
		}
	}
	return result;
}
bool check(vector<biibn_pair> &Ns, vector<biibn_pair> &Nd, vector<biiqueryf> &Fs, vector<biiqueryf> &Fd, int s, int d, double td, double y, int language) {
	int index = 0;
	vector<int> path = convert(language);
	vector<int> FintersectNodes;
	for (const auto& fs : Fs) {
		int fsNode = fs.node;
		for (const auto& fd : Fd) {
			int fdNode = fd.node;
			if (fsNode == fdNode || fsNode == d) {
				FintersectNodes.push_back(fsNode);
				break;
			}
		}
	}
	vector<int> NintersectNodes;
	for (const auto& ns : Ns) {
		int nsNode = ns.node;
		for (const auto& nd : Nd) {
			int ndNode = nd.node;
			if (nsNode == ndNode || nsNode == d) {
				NintersectNodes.push_back(nsNode);
				break;
			}
		}
	}
	if (!FintersectNodes.empty()) {
		if (!NintersectNodes.empty()) {
			int path1;
			int path2;
			double curtime;
			for (const auto& node : NintersectNodes) {
				int currentindex = index;
				for (const auto& ns : Ns)
				{
					if (ns.node == node) {
						path1 = ns.labels;
						curtime = ns.weight;
						break;
					}
				}
				for (const auto& nd : Nd) {
					if (nd.node == node) {
						path2 = nd.labels;
						curtime += nd.weight;
						break;
					}
				}
				if (checkispath(path, index, path1)) {
					if (checkispath(path, index, path2)) {
						if (curtime - td < y) {
							cout << "到达喽！";
							return true;
						}
					}
				}
				index = currentindex;
			}

		}
		for (const auto& node : FintersectNodes) {
			int curindex = index;
			double time;
			int path1;
			int path2;
			for (const auto& it : Fs)
			{
				if (node == it.node)
				{
					time = it.mintravetime(td);
					path1 = it.labels(td);
					if (node == d) {
						if (time < y)
						{
							if (checkispath(path, index, path1)) {
								cout << " 到达喽！";
								return true;
							}
						}
					}
					break;
				}
			}
			for (const auto& it : Fd)
			{
				if (node == it.node)
				{
					time = it.mintravetime(time);
					path2 = it.labels(time);
					break;
				}
			}
			if (time < y)
			{
				if (checkispath(path, index, path1)) {
					if (checkispath(path, index, path2)) {
						cout << " 到达喽！";
						return true;
					}
				}

			}
			index = curindex;
		}
	}
	return false;
}


bool biiquery(BII_index& biindex, int s, int d, double td, double y, double delta, int language) {

	int m = ceil(y / delta);
	for (int k = 1; k <= m; k++)
	{
		vector<biibn_pair> Ns;
		vector<biibn_pair> Nd;
		vector<biiqueryf> Fs;
		vector<biiqueryf> Fd;
		auto pretravetime = [](double td)->double {
			return td;
		};
		auto preminpath = [](double td)->int {
			return 0;
		};
		Ns = FISN(s, k,biindex,td,0,y);
		Nd = BISN(d, k, biindex, td, 0, y);
		Fs = FISF(s, k, biindex,td,y,pretravetime,preminpath);
		Fd = BISF(d, k, biindex, td, y, pretravetime, preminpath);
		if (check(Ns, Nd, Fs, Fd, s, d,td,y,language))
		{
			return true;
		}
	}
	cout << "无法到达";
	return false;
}



#endif
