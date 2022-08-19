#include "LagrangianRelaxation.h"

LagrangianRelaxation::LagrangianRelaxation(int dim, double**& costMatrix)
{
    dimension = dim;
    originalCosts = costMatrix;
}

LagrangianRelaxation::Solution LagrangianRelaxation::subgradient(
    vector<pair<int, int>>& forbiddenEdges, vector<double> u, int ub
)
{
    int pace = 2, updatedPace = 0, iter = 0, lb = 0;
    while(1)
    {
        // update cost matrix
        vector<vector<double>> updatedCostMatrix;
        for(int i = 1; i < dimension; i++)
        {
            vector<double> updatedCostVector(originalCosts[i] + 1, originalCosts[i] + dimension);
            for(int j = i + 1; j < dimension; j++)
            {
                updatedCostVector[j - 1] -= u[i] + u[j];
            }
            updatedCostMatrix.push_back(updatedCostVector);
        }
        for(int i = 0; i < forbiddenEdges.size(); i++)
        {
            updatedCostMatrix[forbiddenEdges[i].first][forbiddenEdges[i].second] = __FLT_MAX__;
        }
        for(int i = 0; i < dimension - 2; i++)
        {
            for(int j = i + 1; j < dimension - 1; j++)
            {
                updatedCostMatrix[j][i] = updatedCostMatrix[i][j];
            }
        }

        // solve MST
        Kruskal solver(updatedCostMatrix);
        double bound = solver.MST(dimension - 1);
        Solution solution(bound, u, solver.getEdges(), forbiddenEdges);
        if(solution.bound > lb)
        {
            lb = solution.bound;
            iter = 0;
        }

        // fix numeration of vertices in edges (updatedCostMatrix ignores depot)
        for(int i = 0; i < solution.edges.size(); i++)
        {
            solution.edges[i].first++;
            solution.edges[i].second++;
        }

        // insert depot
        // calculate the cost of edges with depot
        priority_queue<pair<double, int>> bestVertices;
        if(dimension > 2)
        {
            bestVertices.push(make_pair(originalCosts[0][1] - u[0] - u[1], 1));
            bestVertices.push(make_pair(originalCosts[0][2] - u[0] - u[2], 2));
        }
        for(int i = 3; i < dimension; i++)
        {
            if(originalCosts[0][i] - u[0] - u[i] < bestVertices.top().first)
            {
                bestVertices.pop();
                bestVertices.push(make_pair(originalCosts[0][i] - u[0] - u[i], i));
            }
        }
        // add edges to solution
        solution.edges.push_back(make_pair(0, bestVertices.top().second));
        solution.bound += bestVertices.top().first;
        bestVertices.pop();
        solution.edges.push_back(make_pair(bestVertices.top().second, 0));
        solution.bound += bestVertices.top().first;

        // update u
        double delta = pace*(ub - lb)/(double)((2*solution.edges.size() - 2)*(2*solution.edges.size() - 2)); // check
        for(int i = 0; i < dimension; i++)
        {
            u[i] -= delta;

            if(u[i] < 0)
            {
                u[i] = 0;
            }
        }

        // update pace
        if(++iter >= 30)
        {
            pace /= 2;
            iter = 0;
        }

        // break criteria
        if(updatedPace++ >= 30) // pace == 2^-30
        {
            return solution;
        }
    }
}

int LagrangianRelaxation::getDimension()
{
    return dimension;
}

double LagrangianRelaxation::getCost(int i, int j)
{
    return originalCosts[i][j];
}

vector<vector<double>> LagrangianRelaxation::computeCostsMatrix()
{
    vector<vector<double>> costsMatrix;
    for(int i = 0; i < dimension; i++)
    {
        vector<double> costsVector;
        for(int j = 0; j < dimension; j++)
        {
            costsVector.push_back(originalCosts[i][j]);
        }
        costsMatrix.push_back(costsVector);
    }
    return costsMatrix;
}
