#include "integracio.h"
#include <string.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

double integrar_trapezi_compost(double (*f)(double*,double),double* args, double a, double b, int n)
{
	//This function implements the numerical integration of a function using
	// the composite trapezoidal rule.
	//Input:
		// double (*f)(double*, double): function to integrate. It must be called 
		//   with first a string of arguments, and then the variable to integrate.
		// double* args: arguments to pass to the integrated function.
		// double a: "leftmost" (most negative) limit of the integration interval.
		// double b: "rightmost" (most positive) limit of the integration interval.
		// int n: Number of intervals to use for calculations.
	//Output:
		// Result of the numerical calculation of the integral of the desired 
		// function, with the established arguments, over the given interval, using 
		// 'n' intervals.

	if (n < 2)                                           // Control d'errors.
	{
		fprintf(stderr, "ERROR: integrar_trapezi_compost necessita valor de 'n' major o igual a 2.\n\n");
		return -1;
   }
   
   double pas=(b-a)/n;                                 //Càlcul del pas.
   
   double sumatori = (f(args,a)+f(args,b))/2;          //Inicialització del sumatori.
   int  i=0;
   for(i=0; i<n; i++)                              //Càlcul en bucle del sumatori.
   {
   	sumatori+=f(args,a+i*pas);
   }
   return sumatori*pas;                                //Multipliquem el sumatori pel pas i en retornem els resultats.
}

double integrar_simpson_compost(double (*f)(double*,double),double* args, double a, double b, int n)
{
	//This function implements the numerical integration of a function using
	// the composite Simpson rule.
	//Input:
		// double (*f)(double*, double): function to integrate. It must be called 
		//   with first a string of arguments, and then the variable to integrate.
		// double* args: arguments to pass to the integrated function.
		// double a: "leftmost" (most negative) limit of the integration interval.
		// double b: "rightmost" (most positive) limit of the integration interval.
		// int n: Number of intervals to use for calculations.
	//Output:
		// Result of the numerical calculation of the integral of the desired 
		// function, with the established arguments, over the given interval, using 
		// 'n' intervals.
	if (n < 2)                                           // Control d'errors.
	{
		fprintf(stderr, "ERROR: integrar_simpson_compost necessita valor de 'n' major o igual a 2.\n\n");
		return -1;
   }
   
   double pas=(b-a)/n;                                 //Càlcul del pas.
   
   double sumatori_senars = 0;                         //Inicialització dels sumatoris.
   double sumatori_parells = 0;
   int i;
   for(i=1; i<(n/2)-1; i++)                            //Càlcul en bucle dels sumatoris.
   {
   	sumatori_senars+=f(args,a+(2*i-1)*pas);
   	sumatori_parells+=f(args,a+2*i*pas);
   }
   sumatori_senars+=f(args,a+(2*i)*pas);
   
   double total = f(args,a)+4*sumatori_senars+2*sumatori_parells+f(args,b);
   
   return pas*total/3;                                //Multipliquem el total pel pas/3 i en retornem els resultats.
}

double P(double x, int n) {

    if (n == 0) return 1;
    else if (n == 1) return x;
    else {
        return ((2*n-1)*x*P(x,n-1) - (n-1)*P(x,n-2))/n;
    }
}

double dP(double x, int n) {

    if(n==0) return 0;

    else if (n==1) return 1;

    else return (n/(x*x-1)) * (x*P(x,n) - P(x,n-1));
}

double newton(int n, int i, double (*f)(double*, double), double (*df)(double*, double), double* args){

    double x0, xn, tol;
    x0 = cos(M_PI * (i-1/4)/(n+1/4));
    xn = 99999;
    tol = 1e-8;

    while(fabs(xn - x0) > tol) {
        xn = x0 - f(args,x0)/df(args,x0);
        x0 = xn;
    }

    return xn;
}

double weight(int n, double xi) {
    
    double dPi;
    if(n==1) return 2;
    else if(n==2) return 1;
    else {
        dPi = dP(xi, n);
        return 2/((1-xi*xi)*(dPi*dPi));
    }
}

double integrar_gauss_legendre(double (*f)(double*, double), double (*df)(double*, double), double* args,double a, double b, int n){   
    /*if(n != 2 && n != 5 && n != 10) {
        fprintf(stderr, "The n value is incorrect, it only accepts 2, 5 or 10. \n");
        return -1;
    }*/

    double 	fx,xi,wi,mid_length, mid;
    fx = 0;
    mid_length = (b-a)/2;
    mid = (a+b)/2;

    for (int j=1; j<=n; j++){
        xi = newton(n,j,f, df, args);
        wi = weight(n, xi);

        fx += wi*f(args, mid_length * xi + mid);
    } 

    return mid_length * fx;
}

double integrar_gauss_chebyshev(double (*f)(double*, double), double* args, double a, double b, int n)
{
	//This function implements the numerical integration of a function through
	// the gauss_chebyshev quadrature.
	//Input:
		// double (*f)(double*, double): function to integrate. It must be called 
		//   with first a string of arguments, and then the variable to integrate.
		// double* args: arguments to pass to the integrated function.
		// double a: "leftmost" (most negative) limit of the integration interval.
		// double b: "rightmost" (most positive) limit of the integration interval.
		// int n: Number of intervals to use for calculations.
	//Output:
		// Result of the numerical calculation of the integral of the desired 
		// function, with the established arguments, over the given interval, using 
		// 'n' intervals.
	
	if (n < 2)
	{
		fprintf(stderr, "ERROR: integrar_gauss_chebyshev necessita valor de 'n' major o igual a 2.\n\n");
		return -1;
   }
	
	//Let's create two (double) variables:
		// xi: will contain, on each iteration, the value of the corresponding chebyshev node.
		// sumatori: accumulator for the summation in the for loop below. Initialised at 0.
   double xi, sumatori=0;

	//Now, for each n
   for (int i=1; i<=n; i++) {
      xi = chebyshev_node(i,n);
      sumatori += sqrt((1-pow(xi,2)))*f(args,((b-a)*xi)/2+(a+b)/2);
   }
   return (M_PI*(b-a))/(2*(n+1)) * sumatori;
}

double chebyshev_node(int i,int n)
{
	return cos(((2*i+1)*M_PI)/(2*n));
}
