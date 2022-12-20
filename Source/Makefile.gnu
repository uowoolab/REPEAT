CC = g++
CFLAGS = -Wall -g -O3
CFLAGS = -g -Ofast -finline-limit=200000 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
CFLAGS = -g -Ofast -finline-limit=200000 
CLAPACK_ROOT = ../CLAPACK-3.2.1
# specify directory to header files in CLAPACK
INCLUDE_DIRS = $(CLAPACK_ROOT)/INCLUDE 
# include f2clib dir here
LIB_DIRS = $(CLAPACK_ROOT)/F2CLIBS
# include lapack, blas and f2c libs here.
#LIBS = /usr/lib64/libasan.so.5.0.0 \
#
LIBS = $(CLAPACK_ROOT)/lapack_LINUX.a \
       $(CLAPACK_ROOT)/blas_LINUX.a \
       $(CLAPACK_ROOT)/F2CLIBS/libf2c.a

# leave OPENMP_FLAG empty if you want a serial compilation of REPEAT.  The compiler may complain 
#  about some unrecognized OPENMP directives, but should compile properly.
# For reference: the g++ compiler uses the -fopenmp flag, 
#		 the intel compiler recognizes the -openmp flag
OPENMP_FLAG = -fopenmp
OBJS = memory.o setup.o utilities.o io_processing.o ewald.o symmetry.o vars_arrays.o solver.o

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
	$(CC) $(CFLAGS) repeat_main.cpp $(OBJS) -o repeat.x -I$(INCLUDE_DIRS) -L$(LIB_DIRS) $(LIBS) $(OPENMP_FLAG) 

clean:
	rm -f *.o repeat.x
