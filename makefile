matrix.o:
	g++ -c matrix/matrix.cpp
main: clean matrix.o
	g++ -c main.cpp
	g++ -o main main.o matrix.o
clean:
	rm -f *.o main