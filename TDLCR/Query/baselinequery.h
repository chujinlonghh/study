#ifndef LCR_QUERY_H
#define LCR_QUERY_H
#include"oneindex.h"
#include<iostream>
#include <string>
#include<queue>
#include <vector>
#include <functional> 
#include"TGraph.h"
#include"pathassit.h"


using namespace std;


bool baselinequery(int s, int d, double td, double y, Index &IG,int language)
{	
	vector<bool> visited(IG.getvertexnum(), false);
	queue<pair<int, double>> S;
	S.emplace(s, td);
	visited[s] = true;
	unordered_map<int, vector<int>> indexMap;
	indexMap[s] = { 0 };
	vector<int> path = convert(language);
	while (!S.empty())
	{
		auto vdptime = S.front();
		S.pop();
		int vexsrc = vdptime.first;
		double temp = vdptime.second;
		vector<int> labels = IG.getlabelFromVertex(vexsrc);
		if (!labels.empty()) {
			for (const int& label : labels) {
				vector<pair<int, function<double(double)>>> result = IG.getData(vexsrc, label);
				for ( int itindex : indexMap[vexsrc])
				{
					
					if (!result.empty() && checkispath(path, itindex, label))
					{

						for (auto it = result.begin(); it != result.end(); ++it) {
							pair<int, function<double(double)>> tempres = *it;

							if (tempres.first == d && tempres.second(temp) <= y)
							{
								cout << " 到达喽！" << endl;
								return true;
							}
							else if (tempres.second(temp) > y || visited[tempres.first] == true)
							{
								continue;
							}
							else {
								indexMap[tempres.first].push_back(itindex + 1);
								S.emplace(tempres.first, tempres.second(temp));
								visited[tempres.first] = true;
							}
						}
					}
					result.clear();
				}
			}
		}

	}
	cout << "无法到达";
	return false;

}

#endif 
