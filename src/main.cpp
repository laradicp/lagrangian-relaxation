#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "BranchAndBound.h"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        cout << "Usage: ./lagragian-relax.exe <instance>" << endl;
        exit(1);
    }

    ifstream ubFile;
    ubFile.open("ubs.txt");

    if(!ubFile.is_open())
    {
        cout << "Problem opening ubs.txt for reading." << endl;
        exit(1);
    }

    // get upper bound for instance
    string instanceName;
    int ub = 0;
    while(getline(ubFile, instanceName))
    {
        ubFile >> ub;
        if(instanceName == argv[1])
        {
            break;
        }
        ubFile.ignore();
    }

    if(ub == 0)
    {
        cout << "Upper bound for instance " << argv[1] << " not found in file ubs.txt." << endl;
        exit(1);
    }

    BranchAndBound bb = BranchAndBound(argc, argv);
    LagrangianRelaxation::Solution solution = bb.run(++ub);

    cout << "Lagrangian bound: " << solution.bound << endl;
    for(int i = 0; i < solution.edges.size(); i++)
    {
        cout << "(" << solution.edges[i].first << "," << solution.edges[i].second << ") ";
    }
    cout << endl;

    return 0;
}