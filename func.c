#include "func.h"

double P(double x, int n) {

    if (n == 0) return 1;
    else if (n == 1) return x;
    else {
        return ((2*n-1)*x*P(x,n-1) - (n-1)*P(x,n-2))/n;
    }
}

double dP(double x, int n) {
    return (n/(x*x-1)) * (x*P(x,n) - P(x,n-1));
}