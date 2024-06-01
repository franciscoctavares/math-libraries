complex.o:
	gcc -c complex/complex.c
main: complex.o
	gcc -c main.c
	gcc -o main main.o complex.o
clean:
	rm -f *.o main