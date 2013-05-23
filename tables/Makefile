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
	$(CC) -c $(CFLAGS) $(INC) $*.f

INC = -I../headers

SRC =	chantler_block.f        cromann.f           fcal.f \
	         mucal.f   sasaki_block.f

OBJ =	chantler_block.o        cromann.o           fcal.o \
	         mucal.o   sasaki_block.o

.PHONY: all
all:	$(OBJ)

#=================================================================
# dependencies on header files
#=================================================================
