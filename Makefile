all: test.o Info.o
	g++ -o test test.o Info.o
Info.o: Info.h Info.cpp
	g++ -c Info.cpp
test.o: test.cpp
	g++ -c test.cpp
clean:
	rm *.o
	rm test
