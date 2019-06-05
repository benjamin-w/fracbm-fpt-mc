# If you are experiencing any trouble compiling the code, you are welcome to write 
# b . w a l t e r 1 6 @ imperial.ac.uk
#
# Makefile for compilation of 'fbm'
#
#  - Modify $INCLUDEPATHS and $LIBRARYPATHS as indicated
#  - Install by running $:'make'
#  - (Type make clean to delete object files and compilation target)
#
# COMMENTARY
#
# In a first step, fbm_main.c and fbm_functions.c are compiled using gcc
# Afterwards, they get linked using gfortran
#
# Depending on how your computer is set up it might be necessary to provide the compiler with further informaton regarding header files (for compilation) and libraries (for linking).
#
# For the compilation step, it may be for example necessary to adapt the $INCLUDEPATHS to include the relative paths to 'gsl.h', 'cblas.h', and 'lapacke.h' if they are not already contained in the standard directory /usr/include/ 
# In that case uncomment the line below and modify as appropiate
# INCLUDEPATHS = -I/..path..to..cblas.h -I/..path..to..lapacke.h 
INCLUDEPATHS= 

# For the linking step, it may be necessary to provide the paths of following libraries
# fftw3, math, blas, gsl, gslcblas, lapack, lapacke
# if they are not already stored in /usr/lib or any other standard directory (as they usually are after installing them with a package manager)
#
# If the libraries are not in the default library directories of $FORTRAN (=gfortran), you need to provide them to it by uncommenting and modifying the line below
# LIBRARYPATHS = -L/..path..to..blas -L/..path.to..lapack(e)
LIBRARYPATHS=
#
# EXAMPLE CASES
# 
# Instruction for BREW package manager
#
# If you are using Brew as a package manager, install fftw, lapacke, cblas etc. if you don't have it already via
# brew install fftw
# brew install lapack
# brew install openblas
# and then set
# LAPACKEDIR=/usr/local/opt/lapack
# INCLUDEPATHS = -I/usr/local/opt/openblas/include -I${LAPACKEDIR}/include
# LIBRARYPATHS = -L/${LAPACKEDIR}/lib

####### Below here, it is not necessary to manipulate the make-file

CC = gcc
FORTRAN = gfortran
OPTIM = -O3 
CFLAGS += -Wall 

OBJFILES = fbm_main.o fbm_functions.o

LDFLAGS = -lfftw3 -lm -llapacke -llapack -lblas -lgslcblas -lgsl

TARGET = fbm

$(TARGET): $(OBJFILES)  
	$(FORTRAN) -o $@ $^ $(OPTIM) $(LIBRARYPATHS) $(LDFLAGS) 

.c.o:
	$(CC) $(OPTIM) $(INCLUDEPATHS) $(CFLAGS)   -c -o $@ $^

clean:
	rm -f $(OBJFILES) $(TARGET) *~
