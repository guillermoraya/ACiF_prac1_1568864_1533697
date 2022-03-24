#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "integracio.h"

// N és la funció a integrar (funció de densitat de probabilitat de la 
// distribució Gaussiana).
double N(double* arguments,double x)

{
	double mu=arguments[0];
	double sigma=arguments[1];
	double exponent = -pow(x-mu,2)/(2*sigma*sigma);
	return exp(exponent)/(sigma*sqrt(2*M_PI));
}

double f_1(double x)
{
	return x*x/exp(x);
}

int main(int argc, char **argv)
{
	int n;
	if((argc!=4)&&(argc!=5))
	{
		fprintf(stderr, "ERROR a main: La funció necessita al menys tres arguments; sigma (std), mitjana (mu) i x (valor del qual calculem la probabilitat acumulada). Opcionalment pot rebre'n quatre, on el quart argument seria el nombre d'intervals per als mètodes composts de Simpson i del trapezi.\n\n");
		return -1;
	}
	else if(argc==5)
	{
		n = (int) atof(argv[4]);
	}
	else
	{
		n = 1000000;
	}
  
	double mu = atof(argv[1]);
	double sigma = atof(argv[2]);
	double x = atof(argv[3]);
	
	//We create the arguments array, called "args", in order to call our integration functions later.
	double args[4]={mu,sigma};
	
	printf("PARÀMETRES:\n");
	printf("	Mitjana:             %.8f\n",mu);
	printf("	Desviació estàndard: %.8f\n",sigma);
	printf("	Valor d'x:           %.8f\n\n",x);

	printf("Càlcul de la probabilitat acumulada:\n");
	printf("	-Mètode trapezi compost (n=%d): %.8f \n",n,integrar_trapezi_compost(N,args,2,mu,x,n));
	printf("	-Mètode Simpson compost (n=%d): %.8f \n",n,integrar_simpson_compost(N,args,2,mu,x,n));
	return 0;
}
