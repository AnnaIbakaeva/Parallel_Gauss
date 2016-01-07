#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int** initMatrix(int n, int m)
{
	int **a = new int*[n];
	for (int i = 0; i < n; i++)
		a[i] = new int[m];
	return a;
}

int** fillInMatrix(int** a, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			a[i][j] = i;// rand() % 100;
		}
	}
	return a;
}

void forwardSubstitution(int** a, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{

		}
	}
}

int main(int argc, char* argv[])
{
	int rank, size;
	double t1, t2;
	int rowAmount, columnAmount;

	cout << "Enter row amount:";
	cin >> rowAmount;
	cout << "\n\nEnter column amount:";
	cin >> columnAmount;

	int **a = initMatrix(rowAmount, columnAmount);
	a = fillInMatrix(a, rowAmount, columnAmount);


}