#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "integracio.h"

// TO-DO LIST:
	// TODO: Get output in table format (waaaay more user-friendly)
	// TODO: Implement timing of integration functions.
	// TODO: Calculate errors of integration functions(?).
	

// N és la funció a integrar (funció de densitat de probabilitat de la 
// distribució Gaussiana).
double N(double* arguments,double x){
	// N: Gaussian probability density function (pdf).
		// args: Array of double containing:
			// args[0]=mu: mean of the gaussian distribution.
			// args[1]=sigma: standard deviation of the gaussian distribution.
		// x: double variable of which we intend to calculate the probability.
		
	// We'll read the arguments we receive on 'args'
	double mu=arguments[0]; // <- Mean of the Gaussian
	double sigma=arguments[1];// <- Standard deviation of the Gaussian
	
	// And from then on, we'll just evaluate the pdf as defined by the 
	// arguments above on the variable 'x'.
	double exponent = -pow(x-mu,2)/(2*sigma*sigma);
	return exp(exponent)/(sigma*sqrt(2*M_PI));
}

double f_1(double x){
	return x*x/exp(x);
}

int main(int argc, char **argv){
	int n,gl_n;
	double mu,sigma,x;
	
	if((argc!=4)&&(argc!=5)&&(argc!=1))
	{
		fprintf(stderr, "ERROR a main: La funcio necessita al menys tres arguments; sigma (std), mitjana (mu) i x (valor del qual calculem la probabilitat acumulada). Opcionalment pot rebre'n quatre, on el quart argument seria el nombre d'intervals per als mètodes composts de Simpson i del trapezi.\n\n");
		return -1;
	}
	else
	{
		if(argc==5)
		{
			// If the user provides 4 arguments (argc==5), the last one will be their 
			// desired value for n.
			// We read command line values to double format using the 'atof' function.
			// Since 'n' is supposed to be an integer, we just convert it to integer
			// with the prefix '(int)'.
			n = (int) atof(argv[4]);
		}
		else
		{
			// In case the user doesn't provide us with 4 arguments, we've set this program
			// to use n=10^6 by default.
			n=1000000;
		}
		if(argc!=1)
		{
			// If the user provides us with either four or five arguments through command line,
			//we'll read the first as follows. 
		  	// Since the input is in the char** format, we'll use the 'atof' function (as 
		  	// previously explained), to transform it into doubles (lookup the "atof" function
		  	// for more information).
		  	
		  	// We shall remember that, as per this assignment's requests, the input of "gauss"'s
		  	// main will be comprised of:
		  		// mu   : mean of our N (gaussian  probability distribution).
				// sigma: standard deviation of our N (gaussian probability distribution).
				// x    : value for which we intend to calculate the cumulative probability on a 
				//        gaussian distribution as defined by the values above.
			mu = atof(argv[1]);
			sigma = atof(argv[2]);
			x = atof(argv[3]);
		}
		else
		{
			mu = 100;
			sigma = 15;
			x = 123;
		}
	}
	
	gl_n = n;
	while((gl_n!=2)&&(gl_n!=5)&&(gl_n!=10))
	{
		printf("El nombre d'intervals per a Gauss-Legendre ha de ser 2, 5 o 10.\n");
		printf("El nombre d'intervals per a Gauss-Legendre assignat ara mateix és gl_n=%d.\n",gl_n);
		printf("Si us plau, introduïu un valor acceptable: gl_n=");
		scanf("%d", &gl_n);
		if((gl_n!=2)&&(gl_n!=5)&&(gl_n!=10))
		{
			printf("Error: valor de 'gl_n' no vàlid.\n");
		}
		printf("\n");
	}
	
	// We'll create the arguments array, called "args", in order to call our integration functions later.
	// Said arguments will be, respectively, mu and sigma (as defined above).
	double args[2]={mu,sigma};
	
	// In order to get a more user-friendly interface, let's print out the values of the parameters 
	// (and their meanings, in catalan).
	printf("PARAMETRES:\n");
	printf("	Mitjana:                                      %.8f\n",mu);
	printf("	Desviacio estandard:                          %.8f\n",sigma);
	printf("	Valor d'x:                                    %.8f\n",x);
	printf("	Nº de intervals per al mètode Gauss-Legendre: %d\n",gl_n);
	printf("	Nº de intervals altres mètodes:               %d\n\n",n);

	// And now, for the main course:
	// We'll print the results of calling our integration functions.
	// The arguments passed to each integration function are:
		// N   : the gaussian probability density function, as implemented above our main in this file.
		// args: array of arguments we want to pass to "N" (read comments above for more details).
		// mu  : first (leftmost) point of the interval over which we want to integrate.
		// x   : the last (rightmost) point of the interval over which we want to integrate.
		// n   : the number of intervals to use for the integration method.
		
	// Note that we print 0.5+ the results of the integration.
	// That is because we are calculating the cumulative probability of x for x>mu.
	// Since the cumulative probability of the mean is 0.5, the cumulative probabilty of x>mu
	// Would be 0.5+Integral(probabilityDensityFunction(x)[from mu to x]).
	// (see this assignment's documentation for more details).
	printf("Calcul de la probabilitat acumulada:\n");
	printf("	-Metode trapezi compost      (n=%d): %.8f \n",n,0.5+integrar_trapezi_compost(N,args,mu,x,n));
	printf("	-Metode Simpson compost      (n=%d): %.8f \n",n,0.5+integrar_simpson_compost(N,args,mu,x,n));
	printf("	-Metode de Gauss-Txebixev    (n=%d): %.8f\n", n,0.5+integrar_gauss_chebyshev(N,args,mu,x,n));
	printf("	-Metode de Gauss-Legendre (gl=%d): %.8f \n\n\n\n", gl_n,0.5+integrar_gauss_legendre(N,args,mu,x,gl_n));
	
	printf("Prints per a debugar Gauss-Legendre:\n");
	int gl_n_values[3]={2,5,10};
	for(int i=0;i<3;i++)
	{
		printf("	-Metode de Gauss-Legendre (gl_n=%d): %.8f \n", gl_n_values[i],0.5+integrar_gauss_legendre(N,args,mu,x,gl_n_values[i]));
	}
	return 0;	
}
