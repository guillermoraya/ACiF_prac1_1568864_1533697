#include "integracio.h"
#include "utils.h"
#include <string.h>
#include <stdio.h>
#include <math.h>

double integrar_trapezi_compost(double (*f)(double*,double),double* args,int numArgs, double a, double b, int n)
{
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

double integrar_simpson_compost(double (*f)(double*,double),double* args,int numArgs, double a, double b, int n)
{
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

double integrar_gauss_legendre(double (*f)(double), double a, double b, int n)
{
   double fx,xi,wi,mid_length, mid;
   fx = 0;
   mid_length = (b-a)/2;
   mid = (a+b)/2;

   for (int j=1; j<=n; j++){
      xi = newton(n,j);
      wi = weights(n, xi);

      fx += wi*f(mid_length * xi + mid);
   }

   return mid_length * fx;
}
