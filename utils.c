#define _USE_MATH_DEFINES
#include <math.h>
#include "func.h"
#include "utils.h"

double newton(int n, int i) {

    double x0, xn, tol;
    x0 = cos(M_PI * (i-1/4)/(n+1/4));
    xn = INFINITY;
    tol = 1e-8;

    while(fabs(xn - x0) < tol) {
        xn = x0 - P(x0,n)/dP(x0,n);
        x0 = xn;
    }

    return xn;
}

double weight(int n, double xi) {
    
    double dPi;
    dPi = dP(xi, n);

    return 2/((1-xi*xi)*(dPi*dPi));
}