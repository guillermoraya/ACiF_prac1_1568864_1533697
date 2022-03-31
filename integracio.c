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


    double xn_1, xn, tol;
    xn = cos(M_PI * ((i-1/4)/(n+1/2)));
    tol = 1e-8;

    printf("n: %d, i: %d\n", n, i);
    printf("x0: %8.f\n", xn);
    do {
        xn_1 = xn;
        xn = xn_1 - P(xn_1,n)/dP(xn_1,n);
        printf("xn: %8.f\n", xn);
        printf("xn_1: %8.f\n", xn_1);
        
        printf("P: %8.f\n", P(xn_1,n));
        printf("dP: %8.f\n", dP(xn_1,n));
    } while(fabs(xn - xn_1) < tol);

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

double integrar_gauss_legendre(double (*f)(double*, double), double* args,double a, double b, int n)
{   
	if(n != 2 && n != 5 && n != 10)
	{
		fprintf(stderr, "The n value is incorrect, it only accepts 2, 5 or 10. \n");
		return -1;
	}

	double fx,mid_length, mid, xi[n], wi[n];
	fx = 0;
	mid_length = (b-a)/2;
	mid = (a+b)/2;
	
	if(n==2)
	{
		xi[0] = -1./sqrt(3.);
		xi[1] = 1./sqrt(3.);

		wi[0] = 1.;
		wi[1] = 1.;
	}
	else if (n==5)
	{
		xi[0] = 0.;
		wi[0] = 128./225.;

		xi[1] = -(1./3.)*sqrt(5.-2.*sqrt(10./7.));
		wi[1] = (322.+13.*sqrt(70.))/900.;

		xi[2] = (1./3.)*sqrt(5.-2.*sqrt(10./7.)); 
		wi[2] = wi[1];

		xi[3] = -(1./3.)*sqrt(5.+2.*sqrt(10./7.));
		wi[3] = (322.-13.*sqrt(70.))/900.;
		
		xi[4] = (1./3.)*sqrt(5.+2.*sqrt(10./7.));
		wi[4] = wi[3]; 
	}
	else if (n==10)
	{
		xi[0] = -0.973907;
		wi[0] = 0.0666701555807740;
		
		xi[1] = -0.865063;
		wi[1] = 0.149451725908325;

		xi[2] = -0.67941;
		wi[2] = 0.219086123824010;

		xi[3] = -0.433395;
		wi[3] = 0.269266832578992;

		xi[4] = -0.148874;
		wi[4] = 0.295524255222263;

		xi[5] = 0.148874;
		wi[5] = wi[4];

		xi[6] = 0.433395;
		wi[6] = wi[3];

		xi[7] = 0.67941;
		wi[7] = wi[2];

		xi[8] = 0.865063;
		wi[8] = wi[1];

		xi[9] = 0.973907;
		wi[9] = wi[0];
	}

	for (int i=0; i<n; i++)
	{
		fx += wi[i]*f(args, mid_length * xi[i] + mid);
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
