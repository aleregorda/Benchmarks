.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.machine

default: code2D

OBJECTS3D = output_for_paraview3D.o simplefem3D.o

OBJECTS2D = output_for_paraview2D.o simplefem2D.o

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

code2D:	$(OBJECTS2D)
	$(F90) $(OPTIONS) $(OBJECTS2D) $(LIBS) -o simplefem

code3D:	$(OBJECTS3D)
	$(F90) $(OPTIONS) $(OBJECTS3D) $(LIBS) -o simplefem
