NAME = atoms
VNUM = 250
BINDIR = ./

# compilation information for LINUX using g77
F77 = g77
#CFLAGS = -O2 -ffast-math -m486 -Wall -pedantic -g
CFLAGS = -O2 -ffast-math -m486 -Wall -pedantic -g
LFLAGS =

# -O2:	optimization
# -ffast-math: another optimization flag
# -m486: use i486 primitives (also an optimization flag)
# -pedantic:  strictly follow (sorta) the ansi fortran standard
# -Wall:  issue lots of warnings
# -g: debugging
# -pg:  profiling

.SUFFIXES:
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(CFLAGS) -I$(INC) $*.f $(LFLAGS)


INC = ./headers

LIBS =	readin/readin.a crystl/crystl.a groups/groups.a \
	ascat/ascat.a   mcm/mcm.a       unit/unit.a \
	clustr/clustr.a output/output.a misc/misc.a

# TABLES = tables/mucal.o
# TABLES = tables/chantler_block.o tables/cromann.o tables/fcal.o \
#          tables/mucal.o tables/sasaki_block.o

atoms: atom.o $(LIBS)
	$(F77) $(CFLAGS) -I$(INC) atom.o $(LIBS) tables/mucal.o -o atoms $(LFLAGS)


atom.o:	$(INC)/version.h $(INC)/atparm.h $(INC)/crystl.h \
$(INC)/exafs.h $(INC)/unit.h

######################################################################
### Targets for making the libraries of Atoms, dependencies on source
### and headers

READIN_SRC = \
readin/readin.f readin/atchck.f readin/atinit.f \
readin/atinpt.f readin/dopfix.f readin/getatm.f
readin/readin.a:	$(READIN_SRC)
	$(MAKE) -C readin readin.a

CRYSTL_SRC = \
crystl/crystl.f crystl/basfil.f crystl/basort.f crystl/bperm.f  \
crystl/equipt.f crystl/fperm.f  crystl/genmul.f crystl/ibravl.f \
crystl/metric.f crystl/multip.f crystl/syschk.f
crystl/crystl.a:	$(CRYSTL_SRC) $(INC)/atparm.h $(INC)/crystl.h
	$(MAKE) -C crystl crystl.a

GROUPS_SRC = \
groups/groups.f groups/atspec.f groups/origin.f groups/rh2hex.f \
groups/schfix.f groups/settng.f groups/spcchk.f groups/systm.f
groups/groups.a:	$(GROUPS_SRC) $(INC)/atparm.h $(INC)/crystl.h
	$(MAKE) -C groups groups.a

ASCAT_SRC = \
ascat/ascat.f ascat/makea0.f ascat/refls.f ascat/spcing.f ascat/wrta0.f
ascat/ascat.a:	$(ASCAT_SRC)
	$(MAKE) -C ascat ascat.a

MCM_SRC = mcm/mcm.f mcm/abslen.f mcm/i0.f mcm/mcmast.f mcm/slfabs.f mcm/mucerr.f
mcm/mcm.a:	$(MCM_SRC) $(INC)/atparm.h $(INC)/crystl.h $(INC)/exafs.h
	$(MAKE) -C mcm mcm.a

UNIT_SRC = unit/unit.f unit/p1out.f unit/realcl.f unit/unitcl.f unit/writcl.f
unit/unit.a:	$(UNIT_SRC) $(INC)/atparm.h $(INC)/crystl.h \
$(INC)/unit.h $(INC)/exafs.h $(INC)/version.h
	$(MAKE) -C unit unit.a

CLUSTR_SRC = \
clustr/clustr.f clustr/atheap.f clustr/cellex.f clustr/dist.f \
clustr/ovrlap.f clustr/subshl.f clustr/tetrot.f clustr/trans.f
clustr/clustr.a:	$(CLUSTR_SRC)
	$(MAKE) -C clustr clustr.a

OUTPUT_SRC = output/output.f output/card.f output/feffpr.f output/geout.f
output/output.a:	$(OUTPUT_SRC)
	$(MAKE) -C output output.a

MISC_SRC = \
misc/bwords.f misc/case.f   misc/determ.f misc/fixsym.f misc/getint.f \
misc/getlgc.f misc/getrea.f misc/gettit.f misc/interp.f misc/is2z.f \
misc/isnum.f  misc/istrln.f misc/lower.f  misc/messag.f misc/nofx.f \
misc/nxtunt.f misc/polyft.f misc/ref.f    misc/s2e.f    misc/triml.f \
misc/untab.f  misc/upper.f  misc/volume.f misc/z2s.f    misc/dbglvl.f \
misc/positn.f
misc/misc.a:	$(MISC_SRC)
	$(MAKE) -C misc misc.a

TABLES_SRC = tables/mucal.f
#TABLES_SRC = tables/mucal.f tables/cromann.f tables/fcal.f \
#             tables/sasaki_block.f tables/chantler_block.f


######################################################################
### Making a monolithic distribution file

ALL_SRC = $(READIN_SRC) $(CRYSTL_SRC) $(GROUPS_SRC) $(MCM_SRC) \
          $(UNIT_SRC)   $(CLUSTR_SRC) $(OUTPUT_SRC) \
          $(TABLES_SRC) $(MISC_SRC)
#ALL_SRC = $(READIN_SRC) $(CRYSTL_SRC) $(GROUPS_SRC) $(MCM_SRC) \
#          $(ASCAT_SRC)  $(UNIT_SRC)   $(CLUSTR_SRC) $(OUTPUT_SRC) \
#          $(TABLES_SRC) $(MISC_SRC)

.PHONY:	monolithic ftnchek clean archive
monolithic:
	cat atom.f $(ALL_SRC) > headers/$(NAME)$(VNUM).f
	cd headers ; ./uninclude $(NAME)$(VNUM).f > ../$(NAME)$(VNUM).f ; \
        rm -f $(NAME)$(VNUM).f ; cd ..

monoexec:
	g77 $(CFLAGS) -fno-silent -o $(NAME)$(VNUM) $(NAME)$(VNUM).f

FTNCHEK = -nopretty -notruncation -nonovice -include=$(INC)
ftnchek:
	ftnchek  $(FTNCHEK) atom.f $(ALL_SRC) > from.ftnchek

clean:
	rm -f readin/*.[oa] crystl/*.[oa] groups/*.[oa] \
              ascat/*.[oa]  mcm/*.[oa]    unit/*.[oa]   \
              clustr/*.[oa] output/*.[oa] misc/*.[oa]   \
              atom.o atoms

HEADERS = headers/crystl.h  headers/exafs.h    headers/history \
          headers/uninclude headers/array.maps headers/dhead.h \
          headers/glossary  headers/runtime    headers/unit.h  \
          headers/atparm.h  headers/dversion.h headers/head.h  \
          headers/sample    headers/version.h

MAKEFILES = Makefile        readin/Makefile crystl/Makefile groups/Makefile \
            ascat/Makefile  mcm/Makefile    unit/Makefile   \
            clustr/Makefile output/Makefile misc/Makefile   \

DOC = doc/atoms.texi   doc/atoms.info   doc/atoms.info-1  \
      doc/atoms.info-2 doc/atoms.info-3 doc/atoms_toc.html \
      doc/atoms.html   doc/atoms.dvi    doc/atoms.ps \
      doc/Makefile


archive:
	rm -f atoms$(VNUM).tar.gz
	tar cvf - atom.f $(ALL_SRC) $(HEADERS) $(MAKEFILES) Wishlist \
                  copyright spinel.inp ybco.inp $(DOC)\
		  > atoms$(VNUM).tar
	gzip -fv atoms$(VNUM).tar

FTP = /home/ftp/pub/bruce/atoms-2.50

to_ftp:
	cp -v $(NAME)$(VNUM).f atoms$(VNUM).tar.gz $(DOC) \
              atoms.inp feff.inp $(FTP)
	rm -f $(FTP)/Makefile
