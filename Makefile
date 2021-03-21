all: test.o Info.o cal.h
	g++ -o test test.o Info.o cal.h
Info.o: Info.h Info.cpp
	g++ -c Info.cpp
test.o: test.cpp
	g++ -c test.cpp
clean:
	rm *.o
	rm test
