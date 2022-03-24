update: clean all

all: gauss

gauss: gauss.o integracio.o
	gcc gauss.o integracio.o -o gauss -lm
gauss.o: gauss.c
	gcc -c gauss.c -Wall -o gauss.o -lm
integracio.o: integracio.c
	gcc -c integracio.c -Wall -o integracio.o -lm

clean: 
	rm integracio.o
	rm gauss.o
	rm gauss
