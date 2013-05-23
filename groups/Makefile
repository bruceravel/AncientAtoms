# compilation information for LINUX using g77

CC = g77 -O2 -ffast-math -m486 -pedantic -Wall -g
# CC = g77 -O2 -ffast-math -m486 -pedantic -Wall -g -pg
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

SRC =	atspec.f groups.f origin.f rh2hex.f schfix.f settng.f spcchk.f \
	 systm.f

OBJ =	atspec.o groups.o origin.o rh2hex.o schfix.o settng.o spcchk.o \
	 systm.o

groups.a:	$(OBJ)
	ar sr groups.a $(OBJ)

groups.o:	$(INC)/atparm.h $(INC)/crystl.h

#=================================================================
# dependencies on header files
#=================================================================
