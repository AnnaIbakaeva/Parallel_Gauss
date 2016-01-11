#include <stdio.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdlib>

using namespace std;

double** initMatrix(int n, int m)
{
	double **a = new double*[n];
	for (int i = 0; i < n; i++)
		a[i] = new double[m];
	return a;
}

double** fillInMatrix(double** a, int n, int m)
{
	srand(time(NULL));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			a[i][j] = 1;// rand() % 10 + 1;
		}
		a[i][i] = 2;
		a[i][m - 1] = m;
	}
	return a;
}

double** forwardSubstitution(double** a, int n, int m)
{
	double ** b = initMatrix(n, m);

	for (int j = 0; j < m; j++)
		b[0][j] = a[0][j];

	for (int k = 1; k < n; k++)
	{
		for (int i = k; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				b[i][j] = a[i][j] - ((a[k - 1][j] * a[i][k - 1]) / a[k - 1][k - 1]);
				if (fabs(b[i][j]) < 0.0000001)
					b[i][j] = 0;
			}
		}

		cout << "\nb:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << b[i][j] << " ";
			}
			cout << "\n";
		}

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


double * backSubstitution(double** a, int n, int m)
{
	double * result = new double[m - 1];
	//Обратный ход метода Гаусса
	result[m - 2] = a[n - 1][m - 1] / a[n - 1][m - 2];
	int k = 0;
	for (int i = m - 3; i >= 0; i--)
	{
		int buf = 0;
		for (int j = i + 1; j < m - 1; j++)
			buf += a[n - 2 - k][j] * result[j];
		result[i] = (a[n - 2 - k][m - 1] - buf) / a[i][i];
		k++;
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

	double **a = initMatrix(rowAmount, columnAmount);
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

	double * r = backSubstitution(a, rowAmount, columnAmount);

	cout << "\nResult vector: \n";
	for (int i = 0; i < columnAmount - 1; i++)
		cout << r[i] << " ";

	cin >> rank;
}