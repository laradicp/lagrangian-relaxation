#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>

using namespace std;

class Kruskal{
public:
	Kruskal(vector<vector<double>> dist);

	double MST(int nodes);
	vector<pair<int, int>> getEdges();

private:
	priority_queue<pair<double, pair<int, int>>> graph;
	vector <int> pset;
	vector<pair<int, int>> edges;

	void initDisjoint(int n);
	int findSet(int i);
	void unionSet(int i, int j);
	bool isSameSet(int i, int j);
};

#endif
