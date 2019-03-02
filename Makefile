# Choose which file you want to compile.
filename=main

exexfilename = run-gen-apo-mask
# Path to the moduales
modules=./modules

# PATHS ON MY LAPTOP
include=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/include
healpixlib=/Users/adityarotti/Documents/Work/Software/Healpix_3.20/lib
cfitsiolib=/Users/adityarotti/Documents/Work/Software/cfitsio/lib

######################## GFORTRAN ##########################
FC=gfortran -fbounds-check -ffixed-line-length-none #-Wall -Wextra -Wconversion
F77flags=-c -fopenmp
F90flags= -DGFORTRAN -fno-second-underscore -fopenmp
############################################################

all:bips
bips:bipsobj
	$(FC) -O  -I$(include) $(F90flags) -o $(exexfilename) *.o -L$(healpixlib) -lhealpix  -lhpxgif -L$(cfitsiolib) -lcfitsio
#	$(FC) -O3 *o $(filename).f90 -o nrun

bipsobj:
	$(FC) -c -I$(include) $(F90flags) $(modules)/global.f90 
	$(FC) -c -I$(include) $(F90flags) $(modules)/mask_operations.f90 
	$(FC) -O3 -I$(include) $(F90flags) -c $(filename).f90 -o $(filename).o

clean:
	$(RM) *.o *.mod *~
