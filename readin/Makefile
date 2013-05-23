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

SRC =	atchck.f atinit.f atinpt.f dopfix.f getatm.f readin.f

OBJ =	atchck.o atinit.o atinpt.o dopfix.o getatm.o readin.o


readin.a:	$(OBJ)
	ar sr readin.a $(OBJ)

#=================================================================
# dependencies on header files
#=================================================================
