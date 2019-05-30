# To use this makefile, uncomment commented lines and input your paths to the respective libraries
# Then run 'make' to compile
LDFLAGS = -fPIC -lfftw3 -lm 
OPTIM = -O3 
CFLAGS += -Wall 
OBJFILES = fbm_main.o fbm_functions.o
TARGET = fbm

#INCLUDE = -I/..your path../gsl-2.5 -I/..your path../lapack-3.8.0/LAPACKE/include

CC = gcc
FORTRAN = gfortran
#LIBRARIES = /..your path../lapack-3.8.0/liblapacke.a /..your path../lapack-3.8.0/liblapack.a /..your path../lapack-3.8.0/librefblas.a /..your path../gsl-2.5/cblas/.libs/libgslcblas.a /..your path/gsl-2.5/.libs/libgsl.a

$(TARGET): $(OBJFILES) $(LIBRARIES) 
	$(FORTRAN) $(OPTIM) $(LDFLAGS) -o $@ $^ 

.c.o:
	$(CC) $(OPTIM) $(INCLUDE) $(CFLAGS)   -c -o $@ $^

clean:
	rm -f $(OBJFILES) $(TARGET) *~
