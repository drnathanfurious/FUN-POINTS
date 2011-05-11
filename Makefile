# Makefile for FUN-POINTS

TARGET = out
F90 = ifort
FFLAGS = -m64 -r8 -i8 -fast -opt_report

ifeq ($(strip $(WOODYHOME)),)
CASA=$(HOME)
else
CASA=$(WOODYHOME)
endif

INCLUDE_PATH = -I$(CASA)/include
LIBRARY_PATH = -L$(CASA)/lib

LIBS = -lm -lnlopt

OBJS = data_types.o main.o 

.SUFFIXES:
.SUFFIXES: .o .f90

all : $(TARGET)

funpoints : $(OBJS)
	$(F90) $(FFLAGS) -o funpoints $(OBJS) $(LIBS)

$(TARGET) : $(OBJS)
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(LIBRARY_PATH) $(LIBS)

.f90.o:
	$(F90) $(FFLAGS) -c $*.f90 $(INCLUDE_PATH)

clean:
	rm *.o *.mod out

%.o: %.f90
	$(F90) $(FFLAGS) -c -o $@ $<
