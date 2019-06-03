# Makefile for compilation of 'fbm'
#
#  - Modify $INCLUDEPATHS and $LIBRARYPATHS as indicated
#  - Install by running $:'make'
#  - (Type make clean to delete object files and compilation target)
#
# COMMENTARY
# In a first step, fbm_main.c and fbm_functions.c are compiled using gcc
# Afterwards, they get linked using gfortran
#
# For the compilation step, it is necessary to adapt the $INCLUDEPATHS to include the relative paths to 'gsl.h', 'cblas.h', and 'lapacke.h' if they are not already contained in the standard directory /usr/include/ 
LAPACKEDIR=#../path to LAPACK installation/lapack-3.8.0
INCLUDEPATHS = -I${LAPACKEDIR}/CBLAS/include -I/${LAPACKEDIR}/LAPACKE/include

# For the linking step, it is necessary to provide the paths (by writing them into $LIBRARYPATHS) to following libraries
# fftw3, math, blas, gsl, gslcblas
# if they are not already stored in /usr/lib or any other standard directory (as they usually are after installing them with a package manager)
# Furthermore, one needs to link against the Lapack and Lapacke libraries. These are static by default (liblapack.a and liblapacke.a) and are to be found in the $LAPACK installation directory (see README for details). (Comment: If you are interested in linking with dynamical lapack(e)-libraries, you need to compile them first. See http://theoryno3.blogspot.com/2010/12/compiling-lapack-as-shared-library-in.html for more information)
LIBRARYPATHS = -L${LAPACKEDIR}



####### Below here, it is not necessary to manipulate the make-file

CC = gcc
FORTRAN = gfortran
OPTIM = -O3 
CFLAGS += -Wall 

OBJFILES = fbm_main.o fbm_functions.o

LDFLAGS = -lfftw3 -lm -llapacke -llapack -lblas -lgslcblas -lgsl

TARGET = fbm

$(TARGET): $(OBJFILES)  
	$(FORTRAN) $(OPTIM) $(LIBRARYPATHS) $(LDFLAGS) -o $@ $^ 

.c.o:
	$(CC) $(OPTIM) $(INCLUDEPATHS) $(CFLAGS)   -c -o $@ $^

clean:
	rm -f $(OBJFILES) $(TARGET) *~
