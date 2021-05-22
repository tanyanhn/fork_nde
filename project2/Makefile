all: main.o Multigrid.h
	g++ -o main main.o Multigrid.h -lblas -std=c++11
main.o:
	g++ -c main.cpp -std=c++11
clean:
	rm *.o
	rm main
	rm *.m
test:
	./main Inputfile/test
cleant:
	rm *.o
	rm test
	rm *.m
