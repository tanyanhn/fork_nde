all: main.o Info.o 
	g++ -o main main.o Info.o -std=c++11
Info.o: Info.h Info.cpp
	g++ -c Info.cpp 
main.o: main.cpp 
	g++ -c main.cpp  -std=c++11
Factory.o: Factory.h 
	g++ -c Factory.h -std=c++11
clean:
	rm *.m
run:
	make test1
	./main Inputdata/test_all_1
	matlab -nodesktop -nosplash -r AB_2_Init1_analysis
	./main Inputdata/test_all_2
	matlab -nodesktop -nosplash -r AB_3_Init2_analysis
	./main Inputdata/test_all_3
	matlab -nodesktop -nosplash -r AM_3_Init2_analysis
	./main Inputdata/test_all_4
	matlab -nodesktop -nosplash -r BDF_3_Init2_analysis
	./main Inputdata/test_all_5
	matlab -nodesktop -nosplash -r RK_Init1_analysis
	./main Inputdata/test_all_6
	matlab -nodesktop -nosplash -r AB_4_21331
	./main Inputdata/test_all_7
	matlab -nodesktop -nosplash -r AB_2_3828
	./main Inputdata/test_all_8
	matlab -nodesktop -nosplash -r AM_4_8532
	./main Inputdata/test_all_9
	matlab -nodesktop -nosplash -r AM_2_1914
	./main Inputdata/test_all_10
	matlab -nodesktop -nosplash -r BDF_2_170652
	./main Inputdata/test_all_11
	matlab -nodesktop -nosplash -r BDF_4_1914
	./main Inputdata/test_all_12
	matlab -nodesktop -nosplash -r RK_11376
	./main Inputdata/test_all_13
	matlab -nodesktop -nosplash -r RK_479
	./main Inputdata/test_all_14
	echo All tests have done. Use "make clean" to clean all .m file.
test1:
	./main Inputdata/test1
	matlab -nodesktop -nosplash -r AB_1_24000
	matlab -nodesktop -nosplash -r RK_6000
	echo Use "make clean" to clean all .m file.
test2:
	make test21
	make test22
	make test23
	make test24
test21:
	./main Inputdata/test2_1
	matlab -nodesktop -nosplash -r AB_1_Init2_analysis
	matlab -nodesktop -nosplash -r AB_2_Init1_analysis
	matlab -nodesktop -nosplash -r AB_3_Init2_analysis
	matlab -nodesktop -nosplash -r AB_4_Init1_analysis
	echo Use "make clean" to clean all .m file.
test22:
	./main Inputdata/test2_2
	matlab -nodesktop -nosplash -r AM_3_Init2_analysis
	matlab -nodesktop -nosplash -r AM_5_Init1_analysis
	echo Use "make clean" to clean all .m file.
test23:
	./main Inputdata/test2_3
	matlab -nodesktop -nosplash -r BDF_1_Init1_analysis
	matlab -nodesktop -nosplash -r BDF_3_Init2_analysis
	echo Use "make clean" to clean all .m file.
test24:
	./main Inputdata/test2_4
	matlab -nodesktop -nosplash -r RK_Init1_analysis
	matlab -nodesktop -nosplash -r RK_Init2_analysis
	echo Use "make clean" to clean all .m file.
test31:
	./main Inputdata/test3_1
	./main Inputdata/test3_2
	./main Inputdata/test3_3
	./main Inputdata/test3_4
	echo Use "make clean" to clean all .m file.
test32:
	./main Inputdata/test3_5
	echo Use "make clean" to clean all .m file.
test4:
	./main Inputdata/test4
	echo Use "make clean" to clean all .m file.
story: doc.tex
	make
	xelatex doc
	xelatex doc
	echo You can read report by opening doc.pdf, and run tests by using command "make run".Exit matlab to get the next plot.
cleanc:
	rm *.o
	rm main
	rm *.m
cleanr:
	rm *.pdf
	rm *.aux
	rm *.log
