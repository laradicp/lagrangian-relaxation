#ifndef BRANCHANDBOUND_H
#define BRANCHANDBOUND_H

#include <iostream>
#include <vector>
#include "LagrangianRelaxation.h"
#include "readData.h"

using namespace std;

class BranchAndBound{
public:
    LagrangianRelaxation* relax;

    BranchAndBound(int argc, char** argv);
    LagrangianRelaxation::Solution run(int ub);
	
private:
    int maxEdgesVertex(LagrangianRelaxation::Solution& solution);
};

#endif
