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

SRC =	abslen.f     i0.f    mcm.f mcmast.f slfabs.f mucerr.f

OBJ =	abslen.o     i0.o    mcm.o mcmast.o slfabs.o mucerr.o

mcm.a:	$(OBJ)
	ar sr mcm.a $(OBJ)

#=================================================================
# dependencies on header files
mcm.o:	$(INC)/atparm.h $(INC)/crystl.h $(INC)/exafs.h
#=================================================================
