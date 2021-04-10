all: test2.o Info.o 
	g++ -o test2 test2.o Info.o -std=c++11
Info.o: Info.h Info.cpp
	g++ -c Info.cpp 
test2.o: test2.cpp 
	g++ -c test2.cpp  -std=c++11
Factory.o: Factory.h Factory.cpp
	g++ -c Factory.cpp -std=c++11
clean:
	rm *.o
	rm test2
cleanm:
	rm *.o
	rm *.m
	rm test2
