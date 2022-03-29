#ifndef __INTEGRACIO_H__
#define __INTEGRACIO_H__

#include <math.h>
#include <stdio.h>
double integrar_trapezi_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);
double integrar_simpson_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);

//Utils to implement Gauss-legendre
double P(double x, int n); //Legendre polynomials
double dP(double x, int n); //Derivative of legendre polynomials
double newton(int n, int i, double (*f)(double*, double), double (*df)(double*, double), double* args); // Newton method
double weight(int n, double xi); // Calculates the gauss nodes

double integrar_gauss_legendre(double (*f)(double*, double), double (*df)(double*, double), double *args, int numArgs, double a, double b, int n);
double integrar_gauss_chebyshev(double (*f)(double*, double), double *args, int numArgs, double a, double b, int n);


#endif
