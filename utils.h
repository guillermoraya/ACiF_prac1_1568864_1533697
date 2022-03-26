//Library of the prototypes for some extra functions that we need for the integration methods in integracio.c
#ifndef __UTILS_H__
#define __UTILS_H__
#include <math.h>

//To calculate the ith gauss node ( or ith root of the Lengendre Polynomials)
double newton(int n, int i);
double weights(int n, double xi);

#endif