#include <stdio.h>
#include <mpi.h>
#include <math.h> 
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdlib>

using namespace std;

float** fillInMatrix(float** a, int n, int m)
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

float** initMatrix(int n, int m)
{
	float **a = new float*[n];
	for (int i = 0; i < n; i++)
		a[i] = new float[m];
	a = fillInMatrix(a, n, m);

	cout << "a:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << "\n";
	}

	return a;
}

float * backSubstitution(float** a, int n, int m)
{
	float * result = new float[m - 1];
	//Обратный ход метода Гаусса
	result[m - 2] = a[n - 1][m - 1] / a[n - 1][m - 2];
	for (int i = m - 3; i >= 0; i--)
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
	float * rankRow = new float[columnAmount * rankRowAmount];
	int headPosition = rankRowAmount * rank;

	if (rank == 0)
	{
		float **a = initMatrix(rowAmount, columnAmount);

		t1 = MPI_Wtime();

		//отправляем потокам их строки матрицы
		for (int m = 1; m < size; m++)
		{
			float* sendingRow = new float[columnAmount * rankRowAmount];
			for (int i = 0; i < rankRowAmount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					sendingRow[j + i*columnAmount] = a[i + m*rankRowAmount][j];
				}
			}
			MPI_Send(sendingRow, columnAmount * rankRowAmount, MPI_FLOAT, m, 0, MPI_COMM_WORLD);
		}

		for (int r = 0; r < rankRowAmount; r++)
		{
			//определяем активную строку
			float* activeRow = new float[columnAmount];
			for (int j = 0; j < columnAmount; j++)
			{
				activeRow[j] = a[r][j];
			}
			//отправляем активную строку остальным потокам
			for (int m = 1; m < size; m++)
			{
				MPI_Send(activeRow, columnAmount, MPI_FLOAT, m, m, MPI_COMM_WORLD);
			}
			float ** b = new float*[rankRowAmount];
			for (int i = 0; i < rankRowAmount; i++)
			{
				b[i] = new float[columnAmount];
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
		for (int m = 1; m < size; m++)
		{
			float* reseivedRow = new float[columnAmount * rankRowAmount];
			MPI_Status status;
			MPI_Recv(reseivedRow, columnAmount * rankRowAmount, MPI_FLOAT, m, 1, MPI_COMM_WORLD, &status);

			for (int i = 0; i < rankRowAmount; i++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					a[i + m*rankRowAmount][j] = reseivedRow[j + i*columnAmount];
				}
			}
		}
		t2 = MPI_Wtime();
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
		cout << "\nTime: " << (t2 - t1);

		cout << "\nResult vector: ";
		float * res = backSubstitution(a, rowAmount, columnAmount);
		for (int i = 0; i < columnAmount-1; i++)
			cout << res[i] << " ";
		cout << "\n";
	}
	else
	{
		//принимаем от 0 потока свои строки матрицы
		MPI_Status status;
		MPI_Recv(rankRow, columnAmount * rankRowAmount, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);

		//пока не дошли до своей активной строки
		int h = 0;
		cout << "\nHead position: " << headPosition << "\n";
		int counter = 0;
		int source = 0;
		while (h < headPosition)
		{
			//принимаем активную строку
			MPI_Status status;
			float* activeRow = new float[columnAmount];
			if (counter >= rankRowAmount)
			{
				source++;
				counter = 0;
			}
			MPI_Recv(activeRow, columnAmount, MPI_FLOAT, source, rank, MPI_COMM_WORLD, &status);

			float *tempRow = new float[columnAmount*rankRowAmount];
			for (int i = 0; i < columnAmount*rankRowAmount; i++)
				tempRow[i] = rankRow[i];
			//считаем свои строки для этой активной строки
			cout << "\nrank " << rank << " row:\n";
			for (int i = 0; i < rankRowAmount; i++)
			{
				for (int j = h; j < columnAmount; j++)
				{
					int c = j + (i * columnAmount);
					tempRow[c] = rankRow[c] - ((activeRow[j] * rankRow[i * columnAmount + h]) / activeRow[h]);
					if (fabs(tempRow[c]) < 0.00001)
						tempRow[c] = 0;
					cout << tempRow[c] << " ";
				}
				cout << "\n";
			}

			for (int i = 0; i < columnAmount*rankRowAmount; i++)
			{
				rankRow[i] = tempRow[i];
			}
			h++;
			counter++;
		}

		for (int r = 0; r < rankRowAmount; r++)
		{
			//определяем активную строку
			float* activeRow = new float[columnAmount];
			for (int j = 0; j < columnAmount; j++)
			{
				activeRow[j] = rankRow[j + r*columnAmount];
			}
			//отправляем активную строку остальным потокам
			for (int m = rank + 1; m < size; m++)
			{
				MPI_Send(activeRow, columnAmount, MPI_FLOAT, m, m, MPI_COMM_WORLD);
			}
			float *tempRow = new float[columnAmount*rankRowAmount];
			for (int i = 0; i < columnAmount*rankRowAmount; i++)
				tempRow[i] = rankRow[i];
			//пересчитываем остальные свои строки для текущей активной строки
			for (int k = r + 1; k < rankRowAmount; k++)
			{
				for (int j = 0; j < columnAmount; j++)
				{
					int c = j + (k * columnAmount);
					tempRow[c] = rankRow[c] - ((activeRow[j] * rankRow[k * columnAmount + headPosition + r]) / activeRow[headPosition + r]);
					if (fabs(tempRow[c]) < 0.00001)
						tempRow[c] = 0;
				}
			}
			for (int i = 0; i < columnAmount*rankRowAmount; i++)
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

		MPI_Send(rankRow, columnAmount*rankRowAmount, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
		//если есть потоки с большим рангом, то
		//отправил активную строку потокам с большим рангом
		//посчитал остальные свои строки
		//отправил следующую активную строку
		//посчитал остальные строки, если они остались
		//отправляем 0 потоку свои строки

		//если это последний поток, то считаем свои строки для текущей активной строки
		//меняем активную строки и считаем строки
		//отправляем 0 потоку свои строки
	}

	MPI_Finalize();
}

