#ifndef FUNCTIONSTSP_H
#define FUNCTIONSTSP_H

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "MatrixGraph.h"
#include <algorithm>
#include "Timer.h"

using namespace std;

int* bruteForce(MatrixGraph &);
int* branchAndBound(MatrixGraph &);



#endif // !FUNCTIONSTSP_H
