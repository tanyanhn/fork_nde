all: test.o Multigrid.h
	g++ -o test test.o Multigrid.h
test.o:
	g++ -c test.cpp
clean:
	rm *.o
	rm test
