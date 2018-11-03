#include "pch.h"
#include "MatrixGraph.h"


MatrixGraph::MatrixGraph()
{
	Matrix = new int *[0];
	Matrix[0] = new int[0];
	numOfVertex = 0;
}


MatrixGraph::~MatrixGraph()
{
	if (numOfVertex != 0)
	{
		for (int i = 0; i < numOfVertex; i++)
		{
			delete[] Matrix[i];
		}
		delete[] Matrix;
	}
	numOfVertex = 0;

}

int MatrixGraph::getNumOfVertex()
{
	return numOfVertex;
}

int MatrixGraph::getCost(int v1, int v2)
{
	return Matrix[v1][v2];
}

void MatrixGraph::loadFromFile(string filename)
{
	ifstream fileIn;
	fileIn.open(filename.c_str());

	this->~MatrixGraph();

	fileIn >> numOfVertex;
	Matrix = new int *[numOfVertex];

	for (int i = 0; i < numOfVertex; ++i)
	{

		Matrix[i] = new int[numOfVertex];
		for (int j = 0; j < numOfVertex; ++j)
		{
			fileIn >> Matrix[i][j];
		}
	}
	fileIn.close();
}

void MatrixGraph::showMatrix()
{
	for (int i = 0; i < numOfVertex; i++)
	{
		for (int j = 0; j < numOfVertex; j++)
		{
			cout << setw(3) << Matrix[i][j] << " ";
		}
		cout << endl;
	}
}
