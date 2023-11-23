#ifndef SUB_H
#define SUB_H

#include"TGraph.h"
#include <unordered_set>
#include <random>
#include <ctime>
#include <vector>
#include <fstream>
#include <sstream>

void getfilesubgraph(string path, TGraph& graph, int targetSiz) {

	// 打开输出文件
	std::ofstream outputFile(path);
	if (!outputFile.is_open()) {
		std::cout << "Failed to open output file." << std::endl;
	}
	outputFile << targetSiz << " " << 2 * (targetSiz - 1) << endl;
	std::queue<int> bfsQueue;
	std::unordered_set<int> visited;
	int s = 1;
	int m = 0;
	bfsQueue.push(135698);
	visited.insert(135698);
	while (!bfsQueue.empty() && s < targetSiz) {
		int current = bfsQueue.front();
		bfsQueue.pop();
		for (int j = graph.head[current]; j != -1; j = graph.next[j]) {

			int neighbor = graph.adjv[j];
			auto f = graph.weights[j].f;;
			int la = graph.label[j];
			int num = graph.turnpoint[j];

			if (visited.count(neighbor) == 0) {
				s++;
				visited.insert(neighbor);
				bfsQueue.push(neighbor);
			}
			
			outputFile << current << " " << neighbor << " " << num << endl;
			m++;
			for (int i = 0; i < num; i++)
			{
				outputFile << (*f)[i].t << " " << (*f)[i].w << " ";
			}
			outputFile << endl;
			outputFile << la << endl;
		}
			
		
	}
	outputFile << m << endl;
	outputFile.close();

}




#endif