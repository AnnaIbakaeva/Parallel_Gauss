#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdlib>

using namespace std;

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

double** initMatrix(int n, int m)
{
	double **a = new double*[n];
	for (int i = 0; i < n; i++)
		a[i] = new double[m];
	a = fillInMatrix(a, n, m);

	if (n < 21)
	{
		cout << "a:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				cout << a[i][j] << " ";
			}
			cout << "\n";
		}
	}
	return a;
}

double * backSubstitution(double** a, int n, int m)
{
	double * result = new double[m - 1];
	//Обратный ход метода Гаусса
	result[m - 2] = a[n - 1][m - 1] / a[n - 1][m - 2];
	int k = 0;
	for (int i = m - 3; i >= 0; i--)
	{
		float buf = 0;
		for (int j = i + 1; j < m - 1; j++)
		{
			buf += a[n - 2 - k][j] * result[j];
		}
		result[i] = (a[n - 2 - k][m - 1] - buf) / a[n - 2 - k][n - 2 - k];
		k++;
	}
	return result;
}

int main(int argc, char* argv[])
{
	int rank, size;
	double t1, t2;
	int matrixSize;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int rowAmount, columnAmount;
	if (rank == 0)
	{
		cout << "Enter row amount:\n";
		cin >> rowAmount;
		cout << "\nEnter column amount:\n";
		cin >> columnAmount;
	}

	MPI_Bcast(&rowAmount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&columnAmount, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int residue = rowAmount % size;
	int rankRowAmount = rowAmount / size;

	int headPosition = rankRowAmount * rank;

	if (rank == 0)
	{
		double **a = initMatrix(rowAmount, columnAmount);

		t1 = MPI_Wtime();

		int r = 0;
		//отправляем потокам их строки матрицы
		for (int m = 1; m < size; m++)
		{
			int amount = rankRowAmount;
			if (residue > 0 && size - m <= residue)
			{
				amount++;
			}
			double* sendingRow = new double[columnAmount * amount];
			for (int i = 0; i < amount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					sendingRow[j + i*columnAmount] = a[i + m*rankRowAmount + r][j];
				}
			}
			MPI_Send(sendingRow, columnAmount * amount, MPI_DOUBLE, m, 0, MPI_COMM_WORLD);
			if (amount > rankRowAmount)
				r++;
		}

		for (int r = 0; r < rankRowAmount; r++)
		{
			//определяем активную строку
			double* activeRow = new double[columnAmount];
			for (int j = 0; j < columnAmount; j++)
			{
				activeRow[j] = a[r][j];
			}
			//отправляем активную строку остальным потокам
			for (int m = 1; m < size; m++)
			{
				MPI_Send(activeRow, columnAmount, MPI_DOUBLE, m, m, MPI_COMM_WORLD);
			}
			double ** b = new double*[rankRowAmount];
			for (int i = 0; i < rankRowAmount; i++)
			{
				b[i] = new double[columnAmount];
				for (int j = 0; j < columnAmount; j++)
					b[i][j] = a[i][j];
			}
			//пересчитываем остальные свои строки для текущей активной строки
			for (int k = r + 1; k < rankRowAmount; k++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					b[k][j] = a[k][j] - ((a[r][j] * a[k][r]) / a[r][r]);
					if (fabs(b[k][j]) < 0.00001)
						b[k][j] = 0;
				}
			}
			for (int i = r + 1; i < rankRowAmount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
					a[i][j] = b[i][j];
			}
		}

		//принимаем строки от остальных потоков
		//собираем строки в одну матрицу
		int l = 0;
		for (int m = 1; m < size; m++)
		{
			int amount = rankRowAmount;
			if (residue > 0 && size - m <= residue)
			{
				amount++;
			}
			double* reseivedRow = new double[columnAmount * amount];
			MPI_Status status;
			MPI_Recv(reseivedRow, columnAmount * amount, MPI_DOUBLE, m, 1, MPI_COMM_WORLD, &status);

			for (int i = 0; i < amount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					a[i + m*rankRowAmount + l][j] = reseivedRow[j + i*columnAmount];
				}
			}
			if (amount > rankRowAmount)
				l++;
		}
		
		if (rowAmount <= 20 && columnAmount <= 20)
		{
			cout << "\nResult matrix:\n";
			for (int i = 0; i < rowAmount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
					cout << a[i][j] << " ";
				cout << "\n";
			}
		}	

		//cout << "\nResult vector: ";
		double * res = backSubstitution(a, rowAmount, columnAmount);
		t2 = MPI_Wtime();
		/*for (int i = 0; i < columnAmount - 1; i++)
			cout << res[i] << " ";
		cout << "\n";*/
		cout << "\nTime: " << (t2 - t1);
	}
	else
	{
		//принимаем от 0 потока свои строки матрицы
		MPI_Status status;
		int amount = rankRowAmount;
		if (residue > 0 && size - rank <= residue)
		{
			amount++;
		}
		double * rankRow = new double[columnAmount * amount];
		MPI_Recv(rankRow, columnAmount * amount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		//пока не дошли до своей активной строки
		int h = 0;
		int counter = 0;
		int source = 0;
		bool first = true;
		while (h < headPosition)
		{
			//принимаем активную строку
			MPI_Status status;
			double* activeRow = new double[columnAmount];
			int tmp = rankRowAmount;
			if (residue > 0 && size - source <= residue)
			{				
				tmp++;
			}
			if (counter >= tmp)
			{
				source++;
				counter = 0;
				first = true;
			}
			if (residue > 0 && size - source <= residue)
			{
				if (first)
				{
					headPosition++;
					first = false;
				}
			}
			MPI_Recv(activeRow, columnAmount, MPI_DOUBLE, source, rank, MPI_COMM_WORLD, &status);

			double *tempRow = new double[columnAmount*amount];
			for (int i = 0; i < columnAmount*amount; i++)
				tempRow[i] = rankRow[i];
			//считаем свои строки для этой активной строки
			//cout << "\nrank " << rank << " row:\n";
			for (int i = 0; i < amount; i++)
			{
				for (int j = h; j < columnAmount; j++)
				{
					int c = j + (i * columnAmount);
					tempRow[c] = rankRow[c] - ((activeRow[j] * rankRow[i * columnAmount + h]) / activeRow[h]);
					if (fabs(tempRow[c]) < 0.00001)
						tempRow[c] = 0;
					//cout << tempRow[c] << " ";
				}
				//cout << "\n";
			}

			for (int i = 0; i < columnAmount*amount; i++)
			{
				rankRow[i] = tempRow[i];
			}
			h++;
			counter++;
		}

		for (int r = 0; r < amount; r++)
		{
			//определяем активную строку
			double* activeRow = new double[columnAmount];
			for (int j = 0; j < columnAmount; j++)
			{
				activeRow[j] = rankRow[j + r*columnAmount];
			}
			//отправляем активную строку остальным потокам
			for (int m = rank + 1; m < size; m++)
			{
				MPI_Send(activeRow, columnAmount, MPI_DOUBLE, m, m, MPI_COMM_WORLD);
			}
			double *tempRow = new double[columnAmount*amount];
			for (int i = 0; i < columnAmount*amount; i++)
				tempRow[i] = rankRow[i];
			//пересчитываем остальные свои строки для текущей активной строки
			for (int k = r + 1; k < amount; k++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					int c = j + (k * columnAmount);
					tempRow[c] = rankRow[c] - ((activeRow[j] * rankRow[k * columnAmount + headPosition + r]) / activeRow[headPosition + r]);
					if (fabs(tempRow[c]) < 0.00001)
						tempRow[c] = 0;
				}
			}
			for (int i = 0; i < columnAmount*amount; i++)
			{
				rankRow[i] = tempRow[i];
			}
		}

		//Обратный ход Гаусса
		//if (rank == size - 1)
		//{
		//	for (int k = 0; k < rankRowAmount; k++)
		//	{
		//		float * actRow = new float[columnAmount];
		//		for (int j = 0; j < columnAmount; j++)
		//		{
		//			actRow[j] = rankRow[(rankRowAmount - (k + 1))*columnAmount + j];
		//		}
		//		for (int r = 0; r < rank; r++)
		//		{
		//			//MPI_Send(actRow, columnAmount, MPI_FLOAT, r, r, MPI_COMM_WORLD);
		//		}
		//		float *tempRow = new float[columnAmount*rankRowAmount];
		//		for (int i = 0; i < columnAmount*rankRowAmount; i++)
		//			tempRow[i] = rankRow[i];
		//		//cout << "\nBack Rank row: \n";
		//		for (int i = 0; i < rankRowAmount - (k + 1); i++)
		//		{
		//			for (int j = 0; j < columnAmount; j++)
		//			{
		//				int c = j + (i * columnAmount);
		//				tempRow[c] = rankRow[c] - ((actRow[j] * rankRow[i * columnAmount + (headPosition + rankRowAmount - 1) - k]) / actRow[(headPosition + rankRowAmount - 1) - k]);
		//				if (fabs(tempRow[c]) < 0.00001)
		//					tempRow[c] = 0;
		//				//cout << tempRow[c] << " ";
		//			}
		//			//cout << "\n";
		//		}
		//		//cout << "\n";
		//		for (int i = 0; i < columnAmount*rankRowAmount; i++)
		//		{
		//			rankRow[i] = tempRow[i];
		//		}
		//	}
		//}
		//else
		//{
		//}

		MPI_Send(rankRow, columnAmount*amount, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}

