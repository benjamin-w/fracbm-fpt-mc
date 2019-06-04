# To use this makefile, uncomment commented lines and input your paths to the respective libraries
# Then run 'make' to compile
LDFLAGS = -fPIC -lfftw3 -lm 
OPTIM = -O3 
CFLAGS += -Wall 
OBJFILES = fbm_main.o fbm_functions.o
TARGET = fbm

INCLUDE = -I/usr/local/Cellar/gsl/2.5/include/gsl -I/usr/local/opt/lapack/include

CC = gcc
FORTRAN = gfortran
LIBRARIES = /usr/local/opt/lapack/lib/liblapacke.a /usr/local/opt/lapack/lib/liblapack.a /usr/local/opt/lapack/lib/librefblas.a /usr/local/Cellar/gsl/2.5/cblas/.libs/libgslcblas.a /usr/local/Cellar/gsl/2.5/.libs/libgsl.a

$(TARGET): $(OBJFILES) $(LIBRARIES) 
	$(FORTRAN) $(OPTIM) $(LDFLAGS) -o $@ $^ 

.c.o:
	$(CC) $(OPTIM) $(INCLUDE) $(CFLAGS)   -c -o $@ $^

clean:
	rm -f $(OBJFILES) $(TARGET) *~
