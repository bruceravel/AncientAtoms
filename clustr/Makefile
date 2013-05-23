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


SRC =	atheap.f cellex.f clustr.f   dist.f ovrlap.f subshl.f tetrot.f \
	 trans.f

OBJ =	atheap.o cellex.o clustr.o   dist.o ovrlap.o subshl.o tetrot.o \
	 trans.o

clustr.a:	$(OBJ)
	ar sr clustr.a $(OBJ)

#=================================================================
# dependencies on header files
#=================================================================
