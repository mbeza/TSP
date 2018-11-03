#ifndef MATRIXGRAPH_H
#define MATRIXGRAPH_H

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>

using namespace std;


class MatrixGraph
{
	int **Matrix;
	int numOfVertex;

public:
	MatrixGraph();
	~MatrixGraph();
	void loadFromFile(string filename);
	void showMatrix();
	int getNumOfVertex();
	int getCost(int v1, int v2);
};


#endif // !MATRIXGRAPH_H


