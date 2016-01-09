#include <stdio.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdlib>

using namespace std;

float** initMatrix(int n, int m)
{
	float **a = new float*[n];
	for (int i = 0; i < n; i++)
		a[i] = new float[m];
	return a;
}

float** fillInMatrix(float** a, int n, int m)
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

float** forwardSubstitution(float** a, int n, int m)
{
	float ** b = initMatrix(n, m);

	for (int j = 0; j < m; j++)
		b[0][j] = a[0][j];

	for (int k = 1; k < n; k++)
	{
		for (int i = k; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				b[i][j] = a[i][j] - ((a[k - 1][j] * a[i][k - 1]) / a[k - 1][k - 1]);
				if (b[i][j] < 0.0001)
					b[i][j] = 0;
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
	return b;
}


float * backSubstitution(float** a, int n, int m)
{
	float * result = new float[m - 1];
	//Обратный ход метода Гаусса
	result[m - 2] = a[n - 1][m - 1] / a[n - 1][m - 2];
	for (int i = m - 2; i >= 0; i--)
	{
		int buf = 0;
		for (int j = i + 1; j < m - 1; j++)
			buf += a[i][j] * result[j];
		result[i] = (a[i][n] - buf) / a[i][i];
	}
	return result;
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

	float **a = initMatrix(rowAmount, columnAmount);
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

	a = forwardSubstitution(a, rowAmount, columnAmount);

	cout << "\nresult a:\n";
	for (int i = 0; i < rowAmount; i++)
	{
		for (int j = 0; j < columnAmount; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << "\n";
	}

	float * r = backSubstitution(a, rowAmount, columnAmount);

	cout << "\nResult vector: \n";
	for (int i = 0; i < columnAmount - 1; i++)
		cout << r[i]<<" ";

	cin >> rank;
}