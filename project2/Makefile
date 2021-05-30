all: main.o Multigrid.h
	g++ -o main main.o Multigrid.h -lblas -std=c++11
main.o:
	g++ -c main.cpp -std=c++11
run:
	make
	echo Part 1, display of solution effect.
	echo For nonhomogeneous case,
	./main Inputfile/test01
	matlab -nodesktop -nosplash -r V_cycle_show_128
	echo For homogeneous case,
	./main Inputfile/test02
	matlab -nodesktop -nosplash -r fm_cycle_show_1024
	echo Part 2, analysis of V-cycle for nonhomogeneous case.
	echo For full-weighting and linear interpolation,
	make test231
	echo n = 128,
	matlab -nodesktop -nosplash -r V_cycle_128
	echo n = 256,
	matlab -nodesktop -nosplash -r V_cycle_256
	echo n = 512,
	matlab -nodesktop -nosplash -r V_cycle_512
	echo n = 1024,
	matlab -nodesktop -nosplash -r V_cycle_1024
	echo The maximum norm of the error vector can be found in corresponding output file.
	echo For full-weighting and quadratic interpolation,
	make test232
	echo n = 128,
	matlab -nodesktop -nosplash -r V_cycle_128
	echo n = 256,
	matlab -nodesktop -nosplash -r V_cycle_256
	echo n = 512,
	matlab -nodesktop -nosplash -r V_cycle_512
	echo n = 1024,
	matlab -nodesktop -nosplash -r V_cycle_1024
	echo The maximum norm of the error vector can be found in corresponding output file.
	echo For injection and linear interpolation,
	make test233
	echo n = 128,
	matlab -nodesktop -nosplash -r V_cycle_128
	echo n = 256,
	matlab -nodesktop -nosplash -r V_cycle_256
	echo n = 512,
	matlab -nodesktop -nosplash -r V_cycle_512
	echo n = 1024,
	matlab -nodesktop -nosplash -r V_cycle_1024
	echo The maximum norm of the error vector can be found in corresponding output file.
	echo For injection and quadratic interpolation,
	make test234
	echo n = 128,
	matlab -nodesktop -nosplash -r V_cycle_128
	echo n = 256,
	matlab -nodesktop -nosplash -r V_cycle_256
	echo n = 512,
	matlab -nodesktop -nosplash -r V_cycle_512
	echo n = 1024,
	matlab -nodesktop -nosplash -r V_cycle_1024
	echo The maximum norm of the error vector can be found in corresponding output file.
	echo Part 3, analysis of full multigrid cycle for nonhomogeneous case.
	echo For full-weighting and linear interpolation,
	make test241
	echo For full-weighting and quadratic interpolation,
	make test242
	echo For injection and linear interpolation,
	make test243
	echo For injection and quadratic interpolation,
	echo test244
	echo Part 4, find critical value of relative accuracy.
	make test3
	echo Part 5, analysis of V-cycle for homogeneous case. Only show one group.
	echo For full-weighting and linear interpolation,
	make test421
	echo n = 128,
	matlab -nodesktop -nosplash -r V_cycle_128
	echo n = 256,
	matlab -nodesktop -nosplash -r V_cycle_256
	echo n = 512,
	matlab -nodesktop -nosplash -r V_cycle_512
	echo n = 1024,
	matlab -nodesktop -nosplash -r V_cycle_1024
	echo The maximum norm of the error vector can be found in corresponding output file.
	echo Part 6, analysis of full multigrid cycle for homogeneous case. Only show one group.
	echo For injection and quadratic interpolation,
	make test434
	echo All tests have done. Use "make clean" to clean all .m file.
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
test234:
	./main Inputfile/test2341 > Outputfile/test234_128
	./main Inputfile/test2342 > Outputfile/test234_256
	./main Inputfile/test2343 > Outputfile/test234_512
	./main Inputfile/test2344 > Outputfile/test234_1024
test241:
	./main Inputfile/test241
test242:
	./main Inputfile/test242
test243:
	./main Inputfile/test243
test244:
	./main Inputfile/test244
test3:
	./main Inputfile/test31 > Outputfile/test31
	matlab -nodesktop -nosplash -r V_cycle_1024
	./main Inputfile/test32 > Outputfile/test32
	matlab -nodesktop -nosplash -r V_cycle_1024
	./main Inputfile/test33 > Outputfile/test33
	matlab -nodesktop -nosplash -r V_cycle_1024
	./main Inputfile/test34 > Outputfile/test34
	matlab -nodesktop -nosplash -r V_cycle_1024
test41:
	./main Inputfile/test41
test421:
	./main Inputfile/test4211 > Outputfile/test421_128
	./main Inputfile/test4212 > Outputfile/test421_256
	./main Inputfile/test4213 > Outputfile/test421_512
	./main Inputfile/test4214 > Outputfile/test421_1024
test422:
	./main Inputfile/test4221 > Outputfile/test422_128
	./main Inputfile/test4222 > Outputfile/test422_256
	./main Inputfile/test4223 > Outputfile/test422_512
	./main Inputfile/test4224 > Outputfile/test422_1024
test423:
	./main Inputfile/test4231 > Outputfile/test423_128
	./main Inputfile/test4232 > Outputfile/test423_256
	./main Inputfile/test4233 > Outputfile/test423_512
	./main Inputfile/test4234 > Outputfile/test423_1024
test424:
	./main Inputfile/test4241 > Outputfile/test424_128
	./main Inputfile/test4242 > Outputfile/test424_256
	./main Inputfile/test4243 > Outputfile/test424_512
	./main Inputfile/test4244 > Outputfile/test424_1024
test431:
	./main Inputfile/test431
test432:
	./main Inputfile/test432
test433:
	./main Inputfile/test433
test434:
	./main Inputfile/test434
answers: answers.tex
	xelatex answers
	xelatex answers
	echo You can get detailed answers to problems in Section 1 by opening answers.pdf, and run tests by using command "make run".Exit matlab to get the next plot.
math: math.tex
	xelatex math
	xelatex math
	echo You can get the math document by opening math.pdf, and run tests by using command "make run".Exit matlab to get the next plot.
design: design.tex
	xelatex design
	xelatex design
	echo You can get the degisn document by opening design.pdf, and run tests by using command "make run".Exit matlab to get the next plot.
cleanr:
	rm *.pdf
	rm *.aux
	rm *.log
clean:
	rm *.m
cleanm:
	rm *.o
	rm main
	rm *.m
