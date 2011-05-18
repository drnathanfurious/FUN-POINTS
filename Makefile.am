SUFFIXES = .f90

bin_PROGRAMS = funpoints
funpoints_SOURCES = sort.f90 points.f90 stick.f90 optimization.f90 main.f90 

.f90.o:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c $(FCFLAGS_f90) $< 
