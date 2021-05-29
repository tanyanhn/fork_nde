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
test1:
	./main Inputfile/test1 > data
test22:
	./main Inputfile/test22
test231:
	./main Inputfile/test2311 > Outputfile/test231_128
	./main Inputfile/test2312 > Outputfile/test231_256
	./main Inputfile/test2313 > Outputfile/test231_512
	./main Inputfile/test2314 > Outputfile/test231_1024
test232:
	./main Inputfile/test2321 > Outputfile/test232_128
	./main Inputfile/test2322 > Outputfile/test232_256
	./main Inputfile/test2323 > Outputfile/test232_512
	./main Inputfile/test2324 > Outputfile/test232_1024
test233:
	./main Inputfile/test2331 > Outputfile/test233_128
	./main Inputfile/test2332 > Outputfile/test233_256
	./main Inputfile/test2333 > Outputfile/test233_512
	./main Inputfile/test2334 > Outputfile/test233_1024
cleant:
	rm *.o
	rm test
	rm *.m
