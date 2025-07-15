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
main: clean matrix lp constraint file
	g++ -c src/main.cpp -o obj/main.o
	g++ -o ./bin/main obj/main.o obj/matrix.o obj/lp.o obj/constraint.o obj/file.o
	clear
	./bin/main
clean:
	rm -f obj/* ./bin/main
	rm -rf obj/
	mkdir -p obj