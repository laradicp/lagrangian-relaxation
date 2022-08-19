#include "BranchAndBound.h"

BranchAndBound::BranchAndBound(int argc, char** argv)
{
    int dimension;
    double** originalCosts;
    readData(argc, argv, &dimension, &originalCosts);
    relax = new LagrangianRelaxation(dimension, originalCosts);
}

LagrangianRelaxation::Solution BranchAndBound::run(int ub)
{
    priority_queue<LagrangianRelaxation::Solution> tree;
    vector<double> u(relax->getDimension(), 0);
    
    // initialize tree
    vector<vector<double>> originalCosts = relax->computeCostsMatrix();
    vector<pair<int, int>> forbiddenEdges;
    tree.push(relax->subgradient(forbiddenEdges, u, ub));

    LagrangianRelaxation::Solution bestSolution = tree.top();

    while(!tree.empty())
    {
        LagrangianRelaxation::Solution solution = tree.top();
        tree.pop();
        cout << tree.size() << endl;

        // update best solution
        if(solution.bound < bestSolution.bound)
        {
            bestSolution = solution;
        }

        // find vertex with the maximum number of edges
        int branchVertex = maxEdgesVertex(solution);

        // found a feasible solution for TSP (vertex 0 always has exactly 2 edges)
        if(branchVertex == 0)
        {
            continue;
        }
        
        for(int i = 1; i < solution.edges.size(); i++)
        {
            if(solution.edges[i].first == branchVertex || solution.edges[i].second == branchVertex)
            {
                if(solution.edges[i].first < solution.edges[i].second)
                {
                    solution.forbiddenEdges.push_back(make_pair(solution.edges[i].first,
                        solution.edges[i].second));
                }
                else
                {
                    solution.forbiddenEdges.push_back(make_pair(solution.edges[i].second,
                        solution.edges[i].first));
                }
                tree.push(relax->subgradient(solution.forbiddenEdges, solution.u, ub));
                solution.forbiddenEdges.pop_back();
            }
        }
    }

    return bestSolution;
}

int BranchAndBound::maxEdgesVertex(LagrangianRelaxation::Solution& solution)
{
    vector<int> edgesCounter(relax->getDimension(), 0);
    pair<int, int> maxVertex = make_pair(0, 2); // vertex 0 always has 2 edges
    for(int i = 1; i < solution.edges.size(); i++)
    {
        if(++edgesCounter[solution.edges[i].first] > maxVertex.second)
        {
            maxVertex.first = solution.edges[i].first;
            maxVertex.second = edgesCounter[solution.edges[i].first];
        }
        if(++edgesCounter[solution.edges[i].second] > maxVertex.second)
        {
            maxVertex.first = solution.edges[i].second;
            maxVertex.second = edgesCounter[solution.edges[i].second];
        }
    }
    return maxVertex.first;
}
