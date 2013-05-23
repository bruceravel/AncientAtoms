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

.SUFFIXES: .o .f
.f.o:
	$(CC) -c $(CFLAGS) -I$(INC) $*.f

INC = ../headers

SRC =	crystl.f basfil.f basort.f  bperm.f equipt.f  fperm.f genmul.f \
	ibravl.f metric.f multip.f syschk.f

OBJ =	crystl.o basfil.o basort.o  bperm.o equipt.o  fperm.o genmul.o \
	ibravl.o metric.o multip.o syschk.o

crystl.a:	$(OBJ)
	ar sr crystl.a $(OBJ)

#=================================================================
# dependencies on header files

crystl.o: $(INC)/atparm.h $(INC)/crystl.h
#=================================================================

MISC = $(INC)/case.f   $(INC)/dbglvl.f $(INC)/istrln.f $(INC)/messag.f \
       $(INC)/positn.f $(INC)/upper.f

.PHONY:	assemble
assemble:
	cat $(SRC) $(MISC) > crystl_full.f
#	uninclude crystl_full.f > foobar
