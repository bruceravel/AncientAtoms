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


SRC =	    ascat.f makea0.f  refls.f spcing.f  wrta0.f

OBJ =	    ascat.o makea0.o  refls.o spcing.o  wrta0.o

ascat.a:	$(OBJ)
	ar sr ascat.a $(OBJ)

#=================================================================
# dependencies on header files
#=================================================================
