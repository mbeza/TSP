#include "pch.h"
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <Windows.h>

using namespace std;

int main()
{
	MatrixGraph graph1;
	cout << "Hello\n";
	graph1.loadFromFile("14.txt");

	int* result;
	int* result2;
	




	
	//result = bruteForce(graph1);
	//for (int i = 0; i <= graph1.getNumOfVertex(); ++i)
	//	cout << result[i] << " ";

	//cout << "\nkoszt: " << result[graph1.getNumOfVertex() + 1] << endl<<endl;

	result2 = branchAndBound(graph1);

	for (int i = 0; i <= graph1.getNumOfVertex(); ++i)
		cout << result2[i] << " ";
	cout << "\nkoszt: " << result2[graph1.getNumOfVertex() + 1]<<endl;

	
	
	


	system("pause");
	return 0;

}


