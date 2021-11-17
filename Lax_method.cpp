#define _CRT_SECURE_NO_WARNINGS 
#define _USE_MATH_DEFINES
#include <stdio.h> 
#include <conio.h> 
#include <iostream>
#include <fstream>
#include <locale.h> 
#include <math.h> 
#include <cmath>
using namespace std;


#define pi 3.1415926535
#define l 1. //длина по x
#define N 30 //разбиение по x
#define M 20 //разбиение по времени
#define a 1.
#define k 1.0 // коэффициент запаса для tao

double h = l / N;
double tao = k * h / a;
double T = tao * M;

//значение при x=0
double mu(double t)
{
	return 0;
	
}

//значение при t = 0;
double u_0(double x) {
	return sin(2.*pi*x);
	}

//точное решение
double u_ex(double x, double t)
{
	if (t <= x) {
		return sin(2 * pi*(x - t));
	}
	else {
		return 0;
	}

}

//разбиение по пространству
double grid(double *x) {
	for (int i = 0;i <= N;i++) {
		x[i] = i*h;
	}
	return *x;
}

//разбиение по времени
double time(double *t) {
	for (int i = 0; i <= M; i++)
	{
		t[i] = i*tao;
	}
	return *t;
}

//схема Лакса
double *Lax_method()
{

	double *x = (double*)malloc((N + 1) * sizeof(double));
	grid(x);

	double *t = (double*)malloc((M + 1) * sizeof(double));
	time(t);


	double *u0 = (double*)malloc((N + 1) * sizeof(double));
	double *u = (double*)malloc((N + 1) * sizeof(double));

	for (int j = 0; j <= N; j++)
		u0[j] = u_0(x[j]);

	double eps = 0;
	for (int n = 0; n <= M; n++)
	{
		u[0] = mu(t[n]);
		for (int j = 1; j < N; j++)
		{
			u[j] = 0.5*(u0[j + 1] + u0[j - 1]) - ((a*tao) / (2 * h))*(u0[j + 1] - u0[j - 1]);
		}
		u[N] = u0[N] - ((a*tao) / (2 * h))*(3 * u0[N] - 4 * u0[N - 1] + u0[N - 2]) + ((a*a*tao*tao) / (2 * h*h))*(u0[N] - 2 * u0[N - 1] + u0[N - 2]);
		for (int j = 0; j <= N; j++)
		{
			u0[j] = u[j];
			if (fabs(u[j] - u_ex(x[j], t[n])) > eps)
				eps = fabs(u[j] - u_ex(x[j], t[n]));
		}
	}

	printf("\n");
	printf(" Число шагов N =  %d\n", N);
	printf(" Шаг h = %f\n", h);
	printf(" Коэффицикет запаса k = %f\n", k);
	printf(" Шаг по времени tao %f\n", tao);
	printf(" Погрешность eps = %f\n", eps);
	printf("\n");

	return u;
}

//вывод данных в файл 
void get_results(double *x, double *u, double T) {
	setlocale(LC_NUMERIC, "C");
	FILE *name;
	name = fopen("D:\\new.txt", "w");

	for (int j = 0; j <= N; j++)
	{
		fprintf(name, "%f %f %f\n", x[j], u[j], u_ex(x[j], T));
	}

}

int main()
{
	setlocale(LC_NUMERIC, "C");
	setlocale(LC_ALL, "rus");

	double *x = (double*)malloc((N + 1) * sizeof(double));
	grid(x);

	double *u = (double*)malloc((N + 1) * sizeof(double));
	u = Lax_method();

	get_results(x, u, T);

	getchar();
	return 0;
}




