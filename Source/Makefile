#
# The optimizing Intel compiler options are included here.
#

CC = icpc
CFLAGS = -no-prec-div   -qopt-mem-layout-trans=4 -Wall -g -Ofast  \
  -qopt-subscript-in-range -axmic-avx512 -xHost -falign-functions \
  -qopt-mem-layout-trans=4 -ipo -inline-factor=2000 -inline-max-size=20000 \
  -inline-max-total-size=200000 -qopt-report=5 -qopt-report-phase=all  -qopt-zmm-usage=high 

OFLAGS = -no-prec-div   -qopt-mem-layout-trans=4 -Wall -g -Ofast  \
  -qopt-subscript-in-range -axmic-avx512 -xHost -falign-functions \
  -qopt-mem-layout-trans=4 -ipo -inline-factor=2000 -inline-max-size=20000  \
  -inline-max-total-size=200000 -qopt-report=5 -qopt-report-phase=all  -qopt-zmm-usage=high 

## MAKE SURE TO POINT THE ROOT DIRECTORY OF CLAPACK HERE!! 
CLAPACK_ROOT = ../CLAPACK-3.2.1
# specify directory to header files in CLAPACK
INCLUDE_DIRS = $(CLAPACK_ROOT)/INCLUDE 
# include f2clib dir here
LIB_DIRS = $(CLAPACK_ROOT)/F2CLIBS
# include lapack, blas and f2c libs here.
LIBS = $(CLAPACK_ROOT)/lapack_LINUX.a \
       $(CLAPACK_ROOT)/blas_LINUX.a \
       $(CLAPACK_ROOT)/F2CLIBS/libf2c.a

# leave OPENMP_FLAG empty if you want a serial compilation of REPEAT.  The compiler may complain 
#  about some unrecognized OPENMP directives, but should compile properly.
# For reference: the g++ compiler uses the -fopenmp flag, 
#		 the intel compiler recognizes the -qopenmp flag
 
OPENMP_FLAG = -qopenmp
OBJS = memory.o vars_arrays.o setup.o utilities.o io_processing.o ewald.o symmetry.o solver.o

io_processing.o : vars_arrays.h memory.h utilities.h constants.h symmetry.h io_processing.h io_processing.cpp
	$(CC) $(CFLAGS) -c io_processing.cpp

setup.o : constants.h vars_arrays.h memory.h utilities.h setup.h setup.cpp
	$(CC) $(CFLAGS) -c $(OPENMP_FLAG) setup.cpp

vars_arrays.o : constants.h vars_arrays.h vars_arrays.cpp
	$(CC) $(CFLAGS) -c vars_arrays.cpp

ewald.o : vars_arrays.h memory.h utilities.h ewald.h ewald.cpp
	$(CC) $(CFLAGS) -c ewald.cpp

memory.o : memory.h memory.cpp
	$(CC) $(CFLAGS) -c memory.cpp

symmetry.o : vars_arrays.h memory.h symmetry.h symmetry.cpp
	$(CC) $(CFLAGS) -c symmetry.cpp

solver.o : vars_arrays.h memory.h ewald.h utilities.h symmetry.h io_processing.h solver.h solver.cpp
	$(CC) $(CFLAGS) $(OPENMP_FLAG) -c solver.cpp -I$(INCLUDE_DIRS) 

utilities.o : vars_arrays.h memory.h utilities.h utilities.cpp
	$(CC) $(CFLAGS) -c utilities.cpp

all : repeat_main.cpp $(OBJS) 
	$(CC) $(OFLAGS) repeat_main.cpp $(OBJS) -o repeat.x -I$(INCLUDE_DIRS) -L$(LIB_DIRS) $(LIBS) $(OPENMP_FLAG)  -lstdc++ -lm -lpthread -lc -ldl

clean:
	rm -f *.o repeat.x
