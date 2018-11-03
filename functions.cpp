#include "functions.h"
#include "pch.h"

using namespace std;




int* bruteForce(MatrixGraph & graph)
{
	double czas;
	Timer time;
	cout << "\nBRUTEFORCE\n";

	int numOfVertex = graph.getNumOfVertex();
	int *permutationArray = new int[numOfVertex];
	int *resultArray = new int[numOfVertex+2];
	for (int i = 0; i < numOfVertex; ++i)
	{
		permutationArray[i] = i;
		resultArray[i] = -1;
	}
	int minCost = 2147483111;
	int licznik = 0;

	time.timerStart();
	do {
		int totalCost = 0;
		for (int i = 0; i < numOfVertex - 1; ++i)
			totalCost += graph.getCost(permutationArray[i], permutationArray[i + 1]);

		totalCost += graph.getCost(permutationArray[numOfVertex - 1], permutationArray[0]);
		if (totalCost < minCost)
		{
			minCost = totalCost;
			for (int j = 0; j < numOfVertex; ++j)
				resultArray[j] = permutationArray[j];
		}

		next_permutation(permutationArray, permutationArray + numOfVertex);
		licznik++;
	} while (permutationArray[0] == 0);
	time.timerStop();

	czas = time.calculateTime() / (1000000.0);
	cout <<"czas: "<< czas << endl;
	
	resultArray[numOfVertex] = resultArray[0];
	resultArray[numOfVertex + 1] = minCost;


	
	return resultArray;
}


int* branchAndBound(MatrixGraph & graph)
{
	cout << "\nBRANCH AND BOUND\n";
	double czas;
	Timer time;
	int numOfVertex = graph.getNumOfVertex();
	
	int **Matrix = new int *[numOfVertex + 1]; //tworzymy macierz o wielkosci (n+1)x(n+1)
	
	bool *activeRow = new bool[numOfVertex];
	bool *activeColumn = new bool[numOfVertex];
	
	int *result = new int[numOfVertex];
	
	int activeNumOfVertex = numOfVertex;

	for (int i = 0; i <= numOfVertex; ++i)
		Matrix[i] = new int[numOfVertex + 1];

	for (int i = 0; i < numOfVertex; ++i)	//uzupelniamy macierz danymi z grafu
	{
		for (int j = 0; j < numOfVertex; ++j)
		{
			Matrix[i][j] = graph.getCost(i, j);
			if (i == j) Matrix[i][j] = 2147483111;
		}
	}

	for (int i = 0; i <= numOfVertex; i++)		//brzegowe wartosci (na koncu kolumny i na koncu wiersza uzupelniamy duzymi liczbami
	{
		Matrix[i][numOfVertex] = 2147483111;
		Matrix[numOfVertex][i] = 2147483111;
	}

	for (int i = 0; i < numOfVertex; ++i)
	{
		activeColumn[i] = true;
		activeRow[i] = true;
	}


	int lowerLimit = 0;
	time.timerStart();
	while (activeNumOfVertex >= 2)
	{
		for (int i = 0; i < numOfVertex; ++i)
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[i][j] != -1 && Matrix[i][j] < Matrix[i][numOfVertex]) Matrix[i][numOfVertex] = Matrix[i][j];
			}
		}

		for (int i = 0; i < numOfVertex; ++i)
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[i][j] != -1) Matrix[i][j] -= Matrix[i][numOfVertex];
			}
		}


		

		for (int i = 0; i < numOfVertex; ++i)
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[j][i] != -1 && Matrix[j][i] < Matrix[numOfVertex][i]) Matrix[numOfVertex][i] = Matrix[j][i];
			}
		}

		for (int i = 0; i < numOfVertex; ++i)
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[j][i] != -1) Matrix[j][i] -= Matrix[numOfVertex][i];
			}
		}

		

		int minRow;
		int minColumn;
		int countZeroRow;
		int countZeroColumn;
		int maxElem = -2;
		int row = -1;
		int column = -1;

		for (int i = 0; i < numOfVertex; ++i)
		{
			minRow = 2147483111;
			minColumn = 2147483111;
			countZeroRow = 0;
			countZeroColumn = 0;
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[i][j] == 0) countZeroRow++;
				else
				{
					if (Matrix[i][j] != -1 && Matrix[i][j] < minRow) minRow = Matrix[i][j];
				}

				if (Matrix[j][i] == 0) countZeroColumn++;
				else
				{
					if (Matrix[j][i] != -1 && Matrix[j][i] < minColumn) minColumn = Matrix[j][i];
				}
			}

			if (activeRow[i])
			{
				if (countZeroRow > 1) Matrix[i][numOfVertex] = 0;
				else
				{
					Matrix[i][numOfVertex] = minRow;
					if (minRow >= maxElem)
					{
						maxElem = minRow;
						row = i;
						column = -1;
					}
				}
			}

			if (activeColumn[i])
			{
				if (countZeroColumn > 1) Matrix[numOfVertex][i] = 0;
				else
				{
					Matrix[numOfVertex][i] = minColumn;
					if (minColumn >= maxElem)
					{
						maxElem = minColumn;
						column = i;
						row = -1;
					}
				}
			}


		}

		

		if (row < 0)
		{
			for (int i = 0; i < numOfVertex; ++i)
			{
				if (Matrix[i][column] == 0) row = i;
			}
		}
		else
		{
			for (int i = 0; i < numOfVertex; ++i)
			{
				if (Matrix[row][i] == 0) column = i;
			}
		}

	


		if (activeRow[column] && activeColumn[row]) Matrix[column][row] = 2147483111;

		activeColumn[column] = false;
		activeRow[row] = false;
		result[row] = column;

		

		for (int i = 0; i < numOfVertex; ++i)
		{
			Matrix[row][i] = -1;
			Matrix[i][column] = -1;
			Matrix[i][numOfVertex] = 2147483111;
			Matrix[numOfVertex][i] = 2147483111;
		}


		activeNumOfVertex--;

	}

	int lastRow;
	int lastColumn;

	for (int i = 0; i < numOfVertex; ++i)
	{
		if (activeRow[i]) lastRow = i;
		if (activeColumn[i]) lastColumn = i;
	}

	result[lastRow] = lastColumn;

	time.timerStop();

	czas = time.calculateTime() / (1000000.0);
	cout <<"czas "<< czas << endl;

	int *resultab = new int[numOfVertex + 2];
	int resultiter = 1;

	int iterator = numOfVertex;
	int i = 0;

	int totalCost = 0;

	resultab[0] = i;
	while (iterator > 0)
	{
		resultab[resultiter] = result[i];
		totalCost += graph.getCost(i, result[i]);
		i = result[i];
		iterator--;
		resultiter++;
	}

	resultab[resultiter] = totalCost;

	return resultab;

}