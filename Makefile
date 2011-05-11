# Makefile for FUN-POINTS

TARGET = out
F90 = ifort
FFLAGS = -m64 -r8 -i8 -fast -opt_report

INCLUDE_PATH = -I$(WOODYHOME)/include
LIBRARY_PATH = -L$(WOODYHOME)/lib

LIBS = -lm -lnlopt

OBJS = main.o

.SUFFIXES:
.SUFFIXES: .o .f90

all : $(TARGET)

$(TARGET) : $(OBJS)
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(LIBRARY_PATH) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) -c $*.f90 $(INCLUDE_PATH)

clean:
	rm *.o *.mod out
