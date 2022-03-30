#ifndef __INTEGRACIO_H__
#define __INTEGRACIO_H__

#include <math.h>
#include <stdio.h>
double integrar_trapezi_compost(double (*f)(double*,double),double* args, double interval_min, double interval_max, int n);
double integrar_simpson_compost(double (*f)(double*,double),double* args, double interval_min, double interval_max, int n);
double integrar_gauss_legendre(double (*f)(double*, double), double* args,double a, double b, int n);
double integrar_gauss_chebyshev(double (*f)(double*, double), double *args, double a, double b, int n);
double chebyshev_node(int i,int n);


#endif
