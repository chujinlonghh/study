#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unistd.h>
#include"oneindex.h"
#include"tdp2h.h"
#include"tdptwohquery.h"
#include"baselinequery.h"
#include"biiindex.h"
#include"biiquery.h"
#include <malloc.h>
#include <sys/resource.h>
#include"getsubgraph.h"


int main()
{
	TGraph graph;
	graph.readGraph("input.txt");
	Index IG;
	IG.buildindex(graph);
	LCR_index lcrindex;
	lcrindex.construction(graph);
    BII_index biiindex;
	biiindex.biiconstruct(graph,15);
	ifstream inputFile("queryNY.txt");
	if (!inputFile.is_open()) {
		std::cout << "Failed to open input file." << std::endl;
		return 1;
	}
	string line;
	while (getline(inputFile, line)) {
		std::istringstream iss(line);
		int s, d;
		double td, y;
		if (!(iss >> s >> d >> td >> y)) {
			continue; // 忽略格式不正确的行
		}
		cout << s << "->" << d << "出发时间："<<td<<"gamma:"<<y<<"开始查询：" << endl;
		
		baselinequery(s, d, td, y, IG, 123);
		tdptwohquery(s, d, td, y, lcrindex, 123);
		biiquery(biiindex, s, d, td, y, 10, 123);
	}
	inputFile.close();
	return 0;
}
