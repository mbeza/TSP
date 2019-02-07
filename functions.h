#ifndef FUNCTIONSTSP_H
#define FUNCTIONSTSP_H

#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "MatrixGraph.h"
#include <algorithm>
#include "Timer.h"
#include <stack>

using namespace std;

int* bruteForce(MatrixGraph &);
int* branchAndBound(MatrixGraph &);
int* simulatedAnnealing(MatrixGraph &, int, double, double, double*, int*);
bool checkCost(MatrixGraph &, int *, int*, int*, int*);
int* startPermutation(MatrixGraph &, int);
double P(double, double, double);
int getCost(MatrixGraph &, int*);



#endif // !FUNCTIONSTSP_H
