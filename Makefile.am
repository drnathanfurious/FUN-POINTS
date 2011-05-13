SUFFIXES = .f90
# this only works with the intel compiler... need to check somewhere for it
AM_FCFLAGS = -m64 -r8 -i8 -fast

bin_PROGRAMS = funpoints
funpoints_SOURCES = main.f90 optimization.f90 points.f90

.f90.o:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c $(FCFLAGS_f90) $< 
