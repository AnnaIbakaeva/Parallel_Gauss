#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdlib>

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
	srand(time(NULL));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			a[i][j] = rand() % 10;
		}
	}
	return a;
}

void forwardSubstitution(int** a, int n, int m)
{
	int ** b = initMatrix(n, m);

	for (int j = 0; j < m; j++)
		b[0][j] = a[0][j];

	for (int k = 1; k < n; k++)
	{
		for (int i = k; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				
				
			}
		}

		/*cout << "\nb:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << b[i][j] << " ";
			}
			cout << "\n";
		}*/

		for (int l = 0; l < n; l++)
		{
			for (int t = 0; t < m; t++)
			{
				a[l][t] = b[l][t];
			}
		}
	}	
}

int main(int argc, char* argv[])
{
	int rank, size;
	double t1, t2;
	int rowAmount, columnAmount;

	cout << "Enter row amount:\n";
	cin >> rowAmount;
	cout << "\nEnter column amount:\n";
	cin >> columnAmount;

	int **a = initMatrix(rowAmount, columnAmount);
	a = fillInMatrix(a, rowAmount, columnAmount);

	cout << "a:\n";
	for (int i = 0; i < rowAmount; i++)
	{
		for (int j = 0; j < columnAmount; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << "\n";
	}

	forwardSubstitution(a, rowAmount, columnAmount);
	cin >> rank;
}