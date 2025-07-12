matrix:
	g++ -c src/matrix.cpp -o obj/matrix.o
lp:
	g++ -c src/lp.cpp -o obj/lp.o
ip:
	g++ -c src/ip.cpp -o obj/ip.o
file:
	g++ -c src/file.cpp -o obj/file.o
constraint:
	g++ -c src/constraint.cpp -o obj/constraint.o
main: clean matrix lp ip constraint
	g++ -c main.cpp -o obj/main.o
	g++ -o main obj/main.o obj/matrix.o obj/lp.o obj/ip.o obj/constraint.o
	clear
	./main
test: clean matrix file
	g++ -c test.cpp -o obj/test.o
	g++ -o test obj/test.o obj/matrix.o obj/file.o
	clear
	./test
clean:
	rm -f obj/* main