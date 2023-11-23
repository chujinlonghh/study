#ifndef LABELL_H
#define LABELL_H

#include<iostream>
#include<vector>

vector<int> convert(int number) {
	vector<int> result;

	while (number > 0) {
		result.push_back(number % 10);
		number /= 10;
	}

	reverse(result.begin(), result.end());

	return result;
}

bool checkPath(const vector<int>& path, int& index, int input) {
	if (index> path.size()) {
		return false;
	}
	else if ((index - 1) >= 0 && path[index - 1] == input)
	{
		return true;
	}
	else if (index == path.size()) {
		return false;
	}
	else if(path[index] == input){
		index++;
		return true;
	}
	else {
		return false;
	}
}
bool checkispath(const vector<int>& path, int& index, int input){
vector<int> query = convert(input);
	for (int i = 0; i < query.size(); i++) {
		if (checkPath(path, index, query[i])) {
			continue;
		}
		else {
			return false;
		}
	}
	return true;
}
}
#endif