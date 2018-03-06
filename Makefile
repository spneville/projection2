
F90 = gfortran
F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O2
#-fbacktrace -fbounds-check
FLAGS=-lblas -llapack

PROJECTION2=parsemod.o \
	globalmod.o \
	projmod.o \
	projection2.o

DISPLACE2=parsemod.o \
	globalmod.o \
	dispmod.o \
	displace2.o

MKINFO2=parsemod.o \
	globalmod.o \
	infomod.o \
	mkinfo2.o

ALL=parsemod.o \
	globalmod.o \
	projmod.o \
	dispmod.o \
	infomod.o \
	projection2.o \
	displace2.o \
	mkinfo2.o

projection2: $(PROJECTION2)
	$(F90) $(F90OPTS) $(PROJECTION2) -o projection2 $(FLAGS)
	rm -f *.o *~ *.mod

displace2: $(DISPLACE2)
	$(F90) $(F90OPTS) $(DISPLACE2) -o displace2 $(FLAGS)
	rm -f *.o *~ *.mod

mkinfo2: $(MKINFO2)
	$(F90) $(F90OPTS) $(MKINFO2) -o mkinfo2 $(FLAGS)
	rm -f *.o *~ *.mod

all: $(ALL)
	$(F90) $(F90OPTS) $(PROJECTION2) -o projection2 $(FLAGS)
	$(F90) $(F90OPTS) $(DISPLACE2) -o displace2 $(FLAGS)
	$(F90) $(F90OPTS) $(MKINFO2) -o mkinfo2 $(FLAGS)
	rm -f *.o *~ *.mod

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.F90
	$(F90) -c $(F90OPTS) $<

