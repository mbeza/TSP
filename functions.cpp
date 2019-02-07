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
	int minTotalCost = 2147483111;

	stack <Variant> variants;
	
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

	int countMax;
	int maxElem;
	int row;
	int column;
	int countEqualsMax;
	
	time.timerStart();
	while (activeNumOfVertex >= 2)
	{
		for (int i = 0; i < numOfVertex; ++i)	//szukamy minimum w wierszach
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[i][j] != -1 && Matrix[i][j] < Matrix[i][numOfVertex]) Matrix[i][numOfVertex] = Matrix[i][j];
			}
		}

		for (int i = 0; i < numOfVertex; ++i)	//odejmujemy znalezione minimum z od kazdego elementu w wierszach
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[i][j] != -1) Matrix[i][j] -= Matrix[i][numOfVertex];
			}
		}


		

		for (int i = 0; i < numOfVertex; ++i)	//szukamy minimum w kolumnach
		{
			for (int j = 0; j < numOfVertex; ++j)
			{
				if (Matrix[j][i] != -1 && Matrix[j][i] < Matrix[numOfVertex][i]) Matrix[numOfVertex][i] = Matrix[j][i];
			}
		}

		for (int i = 0; i < numOfVertex; ++i)	//odejmujemy minimum z kolumn
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
		maxElem = -2;
		row = -1;
		column = -1;
		countMax = 0;
		countEqualsMax = 0;

		for (int i = 0; i < numOfVertex; ++i)
		{
			minRow = 2147480111;
			minColumn = 2147480111;
			countZeroRow = 0;
			countZeroColumn = 0;
			for (int j = 0; j < numOfVertex; ++j)	
			{
				if (Matrix[i][j] == 0) countZeroRow++; //liczymy ile zer jest w danym wierszu
				else
				{
					if (Matrix[i][j] != -1 && Matrix[i][j] < minRow) minRow = Matrix[i][j];	//jeœli to nie jest zero, a jest mniejszy ni¿ najmniejszy element minRow, to wstawiamy go jako minimum
				}

				if (Matrix[j][i] == 0) countZeroColumn++;	//liczymy ile zer w kolumnie
				else
				{
					if (Matrix[j][i] != -1 && Matrix[j][i] < minColumn) minColumn = Matrix[j][i]; //analogicznie
				}
			}


			if (activeRow[i])	//jesli wiersz nie jest usuniety
			{
				if (countZeroRow > 1)
				{
					Matrix[i][numOfVertex] = 0;	//sprawdzamy czy sa co najmneij dwa zera, jesli tak to minimum w wierszu jest 0
					minRow = 0;
					if (minRow >= maxElem)		//jesli mniej niz 2 zera, to jako najwieksze minimum oznaczamy najmniejszy element z wiersza, poza zerem
					{
						if (minRow == maxElem && minRow < 2147380111)
						{
							Variant var;
							var.Copy(Matrix, result, numOfVertex, row, column, activeNumOfVertex, activeRow, activeColumn);
							variants.push(var);
								cout << "\n dorzucam na stos " << maxElem<<" z wiersza\n";
							maxElem = minRow;
							row = i;			//tutaj zaznaczamy, ze wybieramy wiersz do usuniecia
							column = -1;		//znak, ze wybralismy wiersz, bo do kolumny wstawiamy -1
						//	countEqualsMax++;

						}
						else
						{
							maxElem = minRow;
							row = i;			//tutaj zaznaczamy, ze wybieramy wiersz do usuniecia
							column = -1;		//znak, ze wybralismy wiersz, bo do kolumny wstawiamy -1

							while (countEqualsMax > 0)
							{
								variants.pop();
								countEqualsMax--;
							}
						}

					}
				}
				else
				{
					Matrix[i][numOfVertex] = minRow; 
					if (minRow >= maxElem)		//jesli mniej niz 2 zera, to jako najwieksze minimum oznaczamy najmniejszy element z wiersza, poza zerem
					{
						if (minRow == maxElem && minRow < 2147380111)
						{
							Variant var;
							var.Copy(Matrix,result,numOfVertex,row,column,activeNumOfVertex,activeRow,activeColumn);
							variants.push(var);
							cout << "\n dorzucam na stos " << maxElem<<" z wiersza\n";
							maxElem = minRow;
							row = i;			//tutaj zaznaczamy, ze wybieramy wiersz do usuniecia
							column = -1;		//znak, ze wybralismy wiersz, bo do kolumny wstawiamy -1
						//	countEqualsMax++;

						}
						else
						{
							maxElem = minRow;
							row = i;			//tutaj zaznaczamy, ze wybieramy wiersz do usuniecia
							column = -1;		//znak, ze wybralismy wiersz, bo do kolumny wstawiamy -1

							while (countEqualsMax > 0)
							{
								variants.pop();
								countEqualsMax--;
							}
						}

					}
				}
			}

			if (activeColumn[i])	//ten sam przypadek z kolumna
			{
				if (countZeroColumn > 1)
				{
					Matrix[numOfVertex][i] = 0;
					minColumn = 0;
					if (minColumn >= maxElem)
					{
						if (minColumn == maxElem && minColumn < 2147380111)
						{
							Variant var;
							var.Copy(Matrix, result, numOfVertex, row, column, activeNumOfVertex, activeRow, activeColumn);
							variants.push(var);
							cout << "\n dorzucam na stos " << maxElem <<" z kolumny\n";
							maxElem = minColumn;
							column = i;		//jesli wybralismy kolumne to wstawiamy jej numer
							row = -1;		//i ustawiamy ze nie wybralismy wiersza
						//	countEqualsMax++;
						}
						else
						{
							maxElem = minColumn;
							column = i;		//jesli wybralismy kolumne to wstawiamy jej numer
							row = -1;		//i ustawiamy ze nie wybralismy wiersza
							while (countEqualsMax > 0)
							{
								variants.pop();
								countEqualsMax--;
							}
						}

					}
				}
				
				else
				{
					Matrix[numOfVertex][i] = minColumn;
					if (minColumn >= maxElem)
					{
						if (minColumn == maxElem && minColumn < 2147380111)
						{
							Variant var;
							var.Copy(Matrix, result, numOfVertex, row, column, activeNumOfVertex, activeRow, activeColumn);
							variants.push(var);
							cout << "\n dorzucam na stos " << maxElem <<" z kolumny\n";
							maxElem = minColumn;
							column = i;		//jesli wybralismy kolumne to wstawiamy jej numer
							row = -1;		//i ustawiamy ze nie wybralismy wiersza
						//	countEqualsMax++;
						}
						else
						{
							maxElem = minColumn;
							column = i;		//jesli wybralismy kolumne to wstawiamy jej numer
							row = -1;		//i ustawiamy ze nie wybralismy wiersza
							while (countEqualsMax > 0)
							{
								variants.pop();
								countEqualsMax--;
							}
						}

					}
				}
			}


		}

		//cout << endl;
		//for (int i = 0; i <= numOfVertex; ++i)
		//{
		//	for (int j = 0; j <= numOfVertex; ++j)
		//	{
		//		if (Matrix[i][j] > 2147380111)
		//		{
		//			cout << setw(3) << " # ";
		//		}
		//		else
		//		{
		//			cout << setw(3) << Matrix[i][j] << " ";
		//		}

		//	}
		//	cout << endl;
		//}
		//cout << endl<<" "<<maxElem<<endl;



		if (false) //miejsce do którego wraca algorytm, gdy s¹ rozwa¿ane inne drogi
		{
		nextvariant:
			Variant newvar;
			newvar = variants.top();
			variants.pop();
			activeNumOfVertex = newvar.activeNumOfVertex;
			numOfVertex = newvar.numOfVertex;
			row = newvar.row;
			column = newvar.column;

			for (int i = 0; i <= numOfVertex; ++i)
			{
				for (int j = 0; j <= numOfVertex; ++j)
				{
					Matrix[i][j] = newvar.Matrix[i][j];
				}
			}

			for (int i = 0; i < numOfVertex; ++i)
			{
				result[i] = newvar.result[i];
				activeColumn[i] = newvar.activeColumn[i];
				activeRow[i] = newvar.activeRow[i];
			}
			
		}

		if (row < 0)	//jesli wybralismy kolumne, czyli row==-1
		{
			for (int i = 0; i < numOfVertex; ++i)
			{
				if (Matrix[i][column] == 0) row = i; //to szukamy do niej wiersza, w ktorym na przecieciu z kolumna ma 0
			}
		}
		else //gdy wybralismy wiersz
		{
			for (int i = 0; i < numOfVertex; ++i)
			{
				if (Matrix[row][i] == 0) column = i;	// to szukamy do niej kolumny
			}
		}

		//cout << column << " " << row << endl;
			if (activeRow[column] && activeColumn[row]) Matrix[column][row] = 2147483111; //blokujemy element "symetryczny"



			activeColumn[column] = false;
			activeRow[row] = false;
			result[row] = column; //do wyniku dajemy usuniety element
			cout << "usuwamy " << row << " " << column << endl;



			for (int i = 0; i < numOfVertex; ++i) //usuwamy wiersz lub kolumne
			{
				Matrix[row][i] = -1;
				Matrix[i][column] = -1;
				Matrix[i][numOfVertex] = 2147483111;
				Matrix[numOfVertex][i] = 2147483111;
			}


			activeNumOfVertex--;
		
	}

	//cout << endl;
	//for (int i = 0; i <= numOfVertex; ++i)
	//{
	//	for (int j = 0; j <= numOfVertex; ++j)
	//	{
	//		if (Matrix[i][j] > 2147380111)
	//		{
	//			cout << setw(3) << " # ";
	//		}
	//		else
	//		{
	//			cout << setw(3) << Matrix[i][j] << " ";
	//		}

	//	}
	//	cout << endl;
	//}

	int lastRow;
	int lastColumn;

	for (int i = 0; i < numOfVertex; ++i)
	{
		if (activeRow[i]) lastRow = i;
		if (activeColumn[i]) lastColumn = i;
	}

	result[lastRow] = lastColumn;


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
	//cout <<endl<< variants.size()<<endl;
	if (totalCost < minTotalCost)
		minTotalCost = totalCost;
	if (!variants.empty())
		goto nextvariant;

	resultab[resultiter] = minTotalCost;
	time.timerStop();

	czas = time.calculateTime() / (1000000.0);
	cout << "czas " << czas << endl;


	return resultab;

}

int* simulatedAnnealing(MatrixGraph & graph, int mode, double Tstart, double deltaT, double* time, int* totalcost)
{
	int numOfVertex = graph.getNumOfVertex();
	int *resultArray = new int[numOfVertex];

	int *actualPermutation = new int[numOfVertex];
	int *previousPermutation;
	int* bestPermutation = new int[numOfVertex];

	int actualCost;
	int prevCost;

	int* actualCostPtr = &actualCost;
	int* prevCostPtr = &prevCost;

	int bestCost;

	Timer countedTime;
	countedTime.timerStart();

	previousPermutation = startPermutation(graph, mode); //pocz¹tkowe rozwi¹zanie
	bestCost = getCost(graph, previousPermutation);

	double T = Tstart;

	while (T > 0.0001)
	{
		//cout << T << endl;
		for (int i = 0; i < numOfVertex; ++i)
		{
			actualPermutation[i] = previousPermutation[i];	//jako aktualne rozwi¹zanie ustawiamy poprzednie rozwi¹zanie
		}

		int candidate1 = rand() % numOfVertex;
		int candidate2;
		do
		{
			candidate2 = rand() % numOfVertex;	//losujemy elementy do zmiany
		} while (candidate1 == candidate2);

		swap(actualPermutation[candidate1], actualPermutation[candidate2]); //zmieniamy dwa losowe elementy

		if (checkCost(graph, actualPermutation, previousPermutation, actualCostPtr, prevCostPtr))	//sprawdzamy, czy koszt aktualnej permutacji(tej po swap) jest mniejszy ni¿ poprzedniej
		{
			for (int i = 0; i < numOfVertex; ++i)	//jesli jest mniejszy, to zamieniamy nasze poprzednie rozwi¹zanie na to, które analizujemy
			{
				previousPermutation[i] = actualPermutation[i];
			}

			int tmpcost = getCost(graph, previousPermutation);
			if (tmpcost < bestCost)
			{
				bestCost = tmpcost;
				for (int i = 0; i < numOfVertex; ++i)
				{
					bestPermutation[i] = previousPermutation[i];	//jako najlepsze rozwi¹zanie ustawiamy aktualne rozwi¹zanie
				}
			}

		}
		else
		{
			double rnd = (double)rand() / RAND_MAX;		//jesli nie jest, to losujemy liczbe z przedzialu (0,1)

			if (rnd < P(T, actualCost, prevCost))		//i sprawdzamy czy jest mniejsza od prawdopodobienstwa P
			{
				for (int i = 0; i < numOfVertex; ++i)	//jesli jest mniejsza, to zamieniamy nasze poprzednie rozwi¹zanie na to, które analizujemy
				{
					previousPermutation[i] = actualPermutation[i];
				}										//a jesli nie jest, to nie zamieniamy naszego poprzedniego rozwiazania, zeby moc sie do niego cofnac i wygenerowac inny ruch
			}

		}
		T *= deltaT;
	}
	countedTime.timerStop();

	*time = countedTime.calculateTime() / (1000000.0);
	*totalcost = bestCost;

	return bestPermutation;
}

bool checkCost(MatrixGraph & graph, int *tab1, int* tab2, int* actual, int* prev)
{
	int numOfVertex = graph.getNumOfVertex();

	int cost1 = 0, cost2 = 0;

	for (int i = 0; i < numOfVertex - 1; ++i)
	{
		cost1 += graph.getCost(tab1[i], tab1[i + 1]);
		cost2 += graph.getCost(tab2[i], tab2[i + 1]);
	}

	cost1 += graph.getCost(tab1[numOfVertex - 1], tab1[0]);
	cost2 += graph.getCost(tab2[numOfVertex - 1], tab2[0]);


	*actual = cost1;
	*prev = cost2;

	if (cost1 < cost2) return true;
	else return false;

}

int* startPermutation(MatrixGraph &graph, int mode)
{
	int numOfVertex = graph.getNumOfVertex();
	int* permutation = new int[numOfVertex];

	if (mode == 1) //simply 0,1,2,3,...,n
	{
		for (int i = 0; i < numOfVertex; ++i)
			permutation[i] = i;
	}
	else
	{
		int* numbers = new int[numOfVertex];
		for (int i = 0; i < numOfVertex; ++i)
			numbers[i] = i;

		int existnumbers = numOfVertex;
		int iterator = 0;
		int randNumber;

		if (mode == 2) //rand
		{
			while (existnumbers > 0)
			{
				randNumber = rand() % existnumbers;
				permutation[iterator] = numbers[randNumber];
				iterator++;
				swap(numbers[randNumber], numbers[existnumbers - 1]);
				existnumbers--;
			}

		}
		if (mode == 3) //greedy algorithm
		{
			bool* isVisited = new bool[numOfVertex];
			for (int i = 0; i < numOfVertex; ++i)
				isVisited[i] = false;

			randNumber = rand() % numOfVertex;
			permutation[iterator] = randNumber;
			isVisited[randNumber] = true;


			for (int i = 0; i < numOfVertex - 1; ++i)
			{
				do
				{
					randNumber = rand() % numOfVertex;
				} while (isVisited[randNumber]);


				int mincost = graph.getCost(permutation[iterator], randNumber);
				int cost;
				int numOfElement = randNumber;

				for (int j = 0; j < numOfVertex; ++j)
				{
					if (isVisited[j]) continue;
					cost = graph.getCost(permutation[iterator], j);
					if (cost < mincost)
					{
						mincost = cost;
						numOfElement = j;
					}
				}
				iterator++;
				permutation[iterator] = numOfElement;
				isVisited[numOfElement] = true;
			}

		}
	}


	return permutation;
}

double P(double T, double actual, double previous)
{
	return exp(((previous - actual) / T));
}

int getCost(MatrixGraph &graph, int* tab)
{
	int numOfVertex = graph.getNumOfVertex();
	int cost = 0;
	for (int i = 0; i < numOfVertex - 1; ++i)
	{
		cost += graph.getCost(tab[i], tab[i + 1]);
	}

	cost += graph.getCost(tab[numOfVertex - 1], tab[0]);

	return cost;
}