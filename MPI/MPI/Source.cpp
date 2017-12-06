#include<mpi.h>
#include<iostream>
#include <windows.h>
using namespace std;
int ProcRank; //Номер процесса
int CountOfProcs; //Число процессов
int pointcount; //Количество известных точек
double* function; //Массив значений функции
double arguments_begin; //Первый аргумент
double step; // Шаг
int algorithm_step; // Шаг подсчета разделенных разностей
int block_size; // Размер блока
int last_process; //Номер процесса содержащий последний элемент
double* deference; //Массив разделенных разностей
double* koef_mas; //Массив коэффициентов
double firsfunct;
void Inizialize(char * path)
{
	if (ProcRank == 0)
	{
		errno_t err;
		FILE* f;
		err = fopen_s(&f, path, "r");
		fscanf_s(f, "%i", &pointcount);
		fscanf_s(f, "%lf", &arguments_begin);
		fscanf_s(f, "%lf", &step);
		function = new double[pointcount];
		for (int i = 0; i < pointcount; i++)
		{
			fscanf_s(f, "%lf", &function[i]);
		}
		fclose(f);
		deference = new double[pointcount];
		algorithm_step = 0;
	}
	MPI_Bcast(&arguments_begin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
double fact(int N)
{
	if (N < 0) // если пользователь ввел отрицательное число
		return 0; // возвращаем ноль
	if (N == 0) // если пользователь ввел ноль,
		return 1; // возвращаем факториал от нуля - не удивляетесь, но это 1 =)
	else // Во всех остальных случаях
		return N * fact(N - 1); // делаем рекурсию.
}
void Data_Sharing()
{
	last_process = CountOfProcs - 1;
	MPI_Bcast(&pointcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		int SimpleBlock = pointcount / CountOfProcs;
		int RestBlock = pointcount - SimpleBlock*(CountOfProcs-1);
		block_size = SimpleBlock;
		double* Send = new double[SimpleBlock];
		for (int i = 1; i < CountOfProcs - 1; i++)
		{
			Send = function + i*SimpleBlock;
			MPI_Send(Send, SimpleBlock, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
		}
		Send = function + SimpleBlock*(CountOfProcs - 1);
		MPI_Send(Send, RestBlock, MPI_DOUBLE, CountOfProcs - 1, CountOfProcs - 1, MPI_COMM_WORLD);
		double* buf = new double[SimpleBlock];
		for (int i = 0; i < SimpleBlock; i++)
		{
			buf[i] = function[i];
		}
		delete function;
		function = buf;
	}
	else
	{
		if (ProcRank != (CountOfProcs - 1))
		{
			block_size= pointcount / CountOfProcs;
		}
		else
		{
			block_size= pointcount - pointcount / CountOfProcs*(CountOfProcs - 1);
		}
		function = new double[block_size];
		MPI_Recv(function, block_size, MPI_DOUBLE, 0, ProcRank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
}
void Show()
{
	if (ProcRank == 0)
	{
		for (int j = 0; j < pointcount; j++)
		{
			printf("%lf \n",deference[j]);
		}
	}
}
void Step_Difference()
{
	double transfer;
	int flag = 0;
	if ((ProcRank != 0)&&(ProcRank<=last_process))
	{
		transfer = function[0];
		MPI_Send(&transfer, 1, MPI_DOUBLE, ProcRank - 1, ProcRank-1, MPI_COMM_WORLD);
	}
	if (ProcRank < last_process)
	{
		MPI_Recv(&transfer, 1, MPI_DOUBLE, ProcRank + 1, ProcRank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
	for (int i = 0; i < block_size-1; i++)
	{
		function[i] = function[i + 1] - function[i];
	}
	if (ProcRank < last_process)
	{
		function[block_size - 1] = transfer - function[block_size - 1];
	}
	if (ProcRank == last_process)
	{
		block_size--;
		if (block_size == 0)
		{
			flag = 1;
		}
	}
	MPI_Bcast(&flag, 1,MPI_INT,last_process, MPI_COMM_WORLD);
	if (flag)
	{
		last_process--;
	}
	if (ProcRank == 0)
	{
		deference[algorithm_step] = function[0];
		algorithm_step = algorithm_step + 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void Shift(double * mas)
{
	for (int i = 0; i < pointcount-1; i++)
	{
		mas[i] = mas[i + 1];
	}
	mas[pointcount - 1] = 0;
}
void CopyMas(double* tmp, double* new_mas)
{
	for (int i = 0; i < pointcount; i++)
	{
		new_mas[i] = tmp[i];
	}
}
void MultyMas(double* mas, double koef)
{
	for (int i = 0; i < pointcount; i++)
	{
		mas[i] = mas[i] * koef;
	}
}
double* SumMas(double* frs, double* sec)
{
	for (int i = 0; i < pointcount; i++)
	{
		frs[i] += sec[i];
	}
	return frs;
}
int FirstChild()
{
	int res = 0;
	res = pointcount / CountOfProcs*ProcRank;
	return res;
}
int LastChild()
{
	int res = 0;
	if (ProcRank != CountOfProcs - 1)
		res = pointcount / CountOfProcs*(ProcRank + 1)-1;
	else
		res = pointcount;
	return res;
}
double * Recurs_Coeff_Arguments()
{
	double* mas = new double[pointcount];
	double* buf = new double[pointcount];
	double* copy = new double[pointcount];
	double firstkoef = 1 / step;
	double secondkoef = -arguments_begin / step;
	for (int i = 0; i < pointcount; i++)
	{
		mas[0] = 0;
		buf[0] = 0;
		copy[0] = 0;
	}
	int n = LastChild();
	int k = FirstChild();
	mas[pointcount-2] = firstkoef;
	mas[pointcount-1] = secondkoef;
	if (k == 1)
	{
		CopyMas(mas, koef_mas);
		MultyMas(koef_mas, deference[k - 1]);
	}
	else
	{
		for (int i = 0; i < pointcount; i++)
		{
			koef_mas[i] = 0;
		}
	}
	if (n > 1)
	{
		for (int i = 2; i <= n; i++)
		{
			CopyMas(mas, buf);
			Shift(mas);
			secondkoef=secondkoef-1;
			MultyMas(mas, firstkoef);
			MultyMas(buf, secondkoef);
			mas = SumMas(mas, buf);
			if (i >= k)
			{
				CopyMas(mas, copy);
				double koef = deference[i - 1] / (fact(i));
				MultyMas(copy, koef);
				koef_mas = SumMas(koef_mas, copy);
			}
		}
	}
	return mas;
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &CountOfProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	Inizialize("input.txt");
	Data_Sharing();
	if (ProcRank == 0)
	{
		firsfunct= function[0];
	}
	for (int i = 0; i < pointcount-1; i++)
	{
		Step_Difference();
	}
	if (ProcRank != 0)
	{
		deference = new double[pointcount - 1];
	}
	MPI_Bcast(deference, pointcount - 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	koef_mas = new double[pointcount];
	Recurs_Coeff_Arguments();
	MPI_Reduce(koef_mas, deference, pointcount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		deference[pointcount-1] = firsfunct;
	}
	Show();
	MPI_Finalize();
	return 0;
}