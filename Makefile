EIGEN_DIR = $(HOME)/lib/eigen3.1.3/
ODEINT_DIR = $(HOME)/lib/boost_1_53_0/
STLIB_DIR = $(HOME)/lib/stlib/src
CC = g++
CLIB = -L/usr/lib64
CINC = -I/usr/include -I. -I $(EIGEN_DIR) -I $(ODEINT_DIR) -I $(STLIB_DIR)
CFLAG = -O2 -msse2 -std=c++0x -fopenmp
MAIN = Test
EXE = run
ADDONS = Include.hh Halo.hh FiniteMethod.hh Advection.hh

$(EXE): Test.o 
	$(CC) $(CFLAG) $(CLIB) -lnetcdf_c++ -o $(EXE) $(<)
$(MAIN).o: $(MAIN).cc $(ADDONS)
	$(CC) $(CFLAG) $(CINC) -c $(<)
clean:
	@ rm -f run
	@ rm -f *.o
	@ rm -f start.in
.PHONY: clean
