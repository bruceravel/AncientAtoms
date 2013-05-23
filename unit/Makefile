# compilation information for LINUX using g77

CC = gfortran -O2 -ffast-math -pedantic -Wall -g -std=legacy
LFLAGS =
CFLAGS =

# -O2:	optimization
# -ffast-math: another optimization flag
# -m486: use i486 primitives (also an optimization flag)
# -pedantic:  strictly follow the ansi fortran standard
# -Wall:  issue lots of warnings
# -g: debugging
# -pg:  profiling

.SUFFIXES:
.SUFFIXES: .o .f
.f.o:
	$(CC) -c $(CFLAGS) -I$(INC) $*.f

INC = ../headers

SRC =	 p1out.f realcl.f   unit.f unitcl.f writcl.f

OBJ =	 p1out.o realcl.o   unit.o unitcl.o writcl.o

unit.a:	$(OBJ)
	ar sr unit.a $(OBJ)
#=================================================================
# dependencies on header files
unit.o:	$(INC)/atparm.h $(INC)/crystl.h $(INC)/unit.h $(INC)/exafs.h \
$(INC)/version.h
#=================================================================
