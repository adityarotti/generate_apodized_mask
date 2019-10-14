# Choose which file you want to compile.
filename=pymain

# Path to the moduales
modules=./modules

# PATHS ON MY LAPTOP
include=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/include
healpixlib=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/lib
cfitsiolib=/Users/adityarotti/Documents/Work/Software/cfitsio/lib

LDFLAGS= -L$(healpixlib) -lhealpix  -lhpxgif -L$(cfitsiolib) -lcfitsio -lgomp
######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -march=native -cpp -ffixed-line-length-none #-Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN -fno-second-underscore #-fopenmp
############################################################

all:bips
bips:bipsobj
	f2py -c --fcompiler=gfortran $(LDFLAGS) *.o -m $(filename) $(filename).f90
	make clean

bipsobj:
	$(FC) -c -I$(include) $(F90flags) $(modules)/global_py.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/mask_operations_py.f90 

clean:
	$(RM) *.o *.mod *~
