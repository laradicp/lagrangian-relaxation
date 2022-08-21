#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include "Kruskal.h"

using namespace std;

class LagrangianRelaxation{
public:
	struct Solution {
		double bound;
		vector<double> u;
		vector<pair<int, int>> edges;
		vector<pair<int, int>> forbiddenEdges;

		Solution(double bound, vector<double>& u, vector<pair<int, int>> edges, vector<pair<int, int>> forbiddenEdges)
		{
			this->bound = bound;
			this->u = u;
			this->edges = edges;
			this->forbiddenEdges = forbiddenEdges;
		}

		bool operator<(const Solution &sol) const
		{
			return bound < sol.bound;
		}
	};

	LagrangianRelaxation(int dim, double**& costMatrix);
	Solution subgradient(vector<pair<int, int>>& forbiddenEdges, vector<double> u, int ub);
	int getDimension();
	double getCost(int i, int j);
	vector<vector<double>> computeCostsMatrix();

private:
	int dimension;
	double** originalCosts;
};

#endif
