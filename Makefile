all: main.o Info.o 
	g++ -o main main.o Info.o -std=c++11
Info.o: Info.h Info.cpp
	g++ -c Info.cpp 
main.o: main.cpp 
	g++ -c main.cpp  -std=c++11
Factory.o: Factory.h Factory.cpp
	g++ -c Factory.cpp -std=c++11
clean:
	rm *.o
	rm main
cleanm:
	rm *.o
	rm *.m
	rm main
run:
	./main Inputdata/Input_test
story: doc.tex
	xelatex doc
	xelatex doc
cleanr:
	rm *.pdf
	rm *.aux
	rm *.log
