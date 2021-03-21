all: test.o Info.o cal.h
	g++ -o test test.o Info.o cal.h -std=c++11
Info.o: Info.h Info.cpp
	g++ -c Info.cpp 
test.o: test.cpp
	g++ -c test.cpp -std=c++11
clean:
	rm *.o
	rm *.m
	rm test
