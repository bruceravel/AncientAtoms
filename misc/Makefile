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

SRC =	bwords.f   case.f determ.f fixsym.f getint.f getlgc.f getrea.f \
	gettit.f interp.f   is2z.f  isnum.f istrln.f  lower.f messag.f \
	  nofx.f nxtunt.f polyft.f    ref.f    s2e.f  triml.f  untab.f \
	 upper.f volume.f    z2s.f dbglvl.f positn.f

OBJ =	bwords.o   case.o determ.o fixsym.o getint.o getlgc.o getrea.o \
	gettit.o interp.o   is2z.o  isnum.o istrln.o  lower.o messag.o \
	  nofx.o nxtunt.o polyft.o    ref.o    s2e.o  triml.o  untab.o \
	 upper.o volume.o    z2s.o dbglvl.o positn.o

misc.a:	$(OBJ)
	ar sr misc.a $(OBJ)

#=================================================================
# dependencies on header files
dbglvl.o:	$(INC)/atparm.h
#=================================================================
