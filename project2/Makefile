all: test.o Multigrid.h
	g++ -o test test.o Multigrid.h -lblas -std=c++11
test.o:
	g++ -c test.cpp -std=c++11
clean:
	rm *.o
	rm test
	rm *.m

