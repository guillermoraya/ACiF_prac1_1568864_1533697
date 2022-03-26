#ifndef __INTEGRACIO_H__
#define __INTEGRACIO_H__

#include <math.h>
#include <stdio.h>
double integrar_trapezi_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);
double integrar_simpson_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);

//Utils to implement Gauss-legendre
double P(double x, int n);
double dP(double x, int n);
double newton(int n, int i);
double weights(int n, double xi);

double integrar_gauss_legendre(double (*f)(double), double a, double b, int n);
double integrar_gauss_chebyshev(double (*f)(double), int n);


#endif
