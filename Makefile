EIGEN_DIR = $(HOME)/lib/eigen3.1.3/
ODEINT_DIR = $(HOME)/lib/boost_1_53_0/
STLIB_DIR = $(HOME)/lib/stlib/src
CC = g++
CLIB = -L/usr/lib64
CINC = -I/usr/include -I. -I $(EIGEN_DIR) -I $(ODEINT_DIR) -I $(STLIB_DIR)
CFLAG = -O2 -msse2 -std=c++0x -fopenmp
EXE = run
ADDONS = Include.hh Boundary.hh Tools.hh \
		 FiniteMethod.hh Advection.hh\
		 OdeSystem.hh ShallowWater.hh SlabHeating.hh ColumnHeating.hh\
		 MicroPhysics.hh

$(EXE): Main.o 
	$(CC) $(CFLAG) $(CLIB) -lnetcdf_c++ -o $(EXE) $(<)
Main.o: Main.cc $(ADDONS)
	$(CC) $(CFLAG) $(CINC) -c $(<)
clean:
	@ rm -f run
	@ rm -f *.o
	@ rm -f start.in
.PHONY: clean
