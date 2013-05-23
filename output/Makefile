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

SRC =	  card.f feffpr.f  geout.f output.f

OBJ =	  card.o feffpr.o  geout.o output.o

output.a:	$(OBJ)
	ar sr output.a $(OBJ)

#=================================================================
# dependencies on header files
#=================================================================
