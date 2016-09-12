CC=g++
CFLAGS=-std=c++0x

help:
	@echo 'Makefile for footprint problem and visualizations                      '
	@echo '                                                                       '
	@echo 'Usage:                                                                 '
	@echo '   make all                         compile diskcover-all and diskover-plot-solution'
	@echo '   make sample                      compile diskcover-sample           '
	@echo '   make coverage                    compile diskcover-coverage         '
	@echo '   make solution                    compile diskcover-plot-solution	  '
	@echo '   make vis                         runs files and makes visualization '
	@echo '                                                                       '

all: 
	$(CC) $(CFLAGS) main-all.cpp diskcover.cpp Footprint.cpp IO.cpp Point.cpp -l stdc++ -o diskcover-all
	$(CC) $(CFLAGS) main-plot-solution.cpp diskcover.cpp Footprint.cpp IO.cpp Point.cpp -l stdc++ -o diskcover-plot-solution

sample:
	$(CC) $(CFLAGS) main-sample.cpp diskcover.cpp Footprint.cpp IO.cpp Point.cpp -l stdc++ -o diskcover-sample

coverage: 
	$(CC) $(CFLAGS) main-coverage.cpp diskcover.cpp Footprint.cpp IO.cpp Point.cpp -l stdc++ -o diskcover-coverage

solution:
	$(CC) $(CFLAGS) main-plot-solution.cpp diskcover.cpp Footprint.cpp IO.cpp Point.cpp -l stdc++ -o diskcover-plot-solution

vis:
	bash ./runfootprint.sh

