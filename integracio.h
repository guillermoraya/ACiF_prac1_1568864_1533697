#ifndef __INTEGRACIO_H__
#define __INTEGRACIO_H__

#include <math.h>
#include <stdio.h>
double integrar_trapezi_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);
double integrar_simpson_compost(double (*f)(double*,double),double* args,int numArgs, double interval_min, double interval_max, int n);
double integrar_gauss_legendre(double (*f)(double), double a, double b, int n);
double integrar_gauss_chebyshev(double (*f)(double), int n);


#endif
