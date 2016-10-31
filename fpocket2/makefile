#------------------------------------------------------------
# Makefile.
#------------------------------------------------------------

PATH_GSL    = /home/user/gsl/
PATH_OBJ    = obj/
PATH_SRC    = src/
PATH_HEADER = headers/
PATH_BIN    = bin/
PATH_MAN    = man/
PATH_QHULL  = src/qhull/

BINDIR  = /usr/local/bin/
MANDIR  = /usr/local/man/man8/


FPOCKET     = fpocket
TPOCKET		= tpocket
DPOCKET		= dpocket
MDPOCKET	= mdpocket
CHECK		= pcheck
MYPROGS		= $(PATH_BIN)$(FPOCKET) $(PATH_BIN)$(TPOCKET) $(PATH_BIN)$(DPOCKET) $(PATH_BIN)$(MDPOCKET)

CC          = gcc
CCQHULL	    = gcc
LINKER		= gcc
LINKERQHULL = gcc

CGSL        = -DMD_NOT_USE_GSL -I$(PATH_GSL)include
COS         = -DM_OS_LINUX
CDEBUG		= -DMNO_MEM_DEBUG
CWARN       = -Wall -Wextra -Wwrite-strings -Wstrict-prototypes

CFLAGS      = $(CWARN) $(COS) $(CDEBUG) -O2 -g -pg -std=c99#$(CGSL)
QCFLAGS     = -O -g -pg -ansi

LGSL        = -L$(PATH_GSL)lib -lgsl -lgslcblas 
LFLAGS	    = -fno-underscoring -lm

#------------------------------------------------------------
# BINARIES OBJECTS 
#------------------------------------------------------------

QOBJS = $(PATH_QHULL)qvoronoi.o $(PATH_QHULL)user.o $(PATH_QHULL)global.o \
		$(PATH_QHULL)io.o $(PATH_QHULL)geom2.o $(PATH_QHULL)poly2.o \
		$(PATH_QHULL)merge.o $(PATH_QHULL)qhull.o $(PATH_QHULL)geom.o  \
		$(PATH_QHULL)poly.o $(PATH_QHULL)qset.o $(PATH_QHULL)mem.o \
		$(PATH_QHULL)stat.o

CHOBJ = $(PATH_OBJ)check.o $(PATH_OBJ)psorting.o $(PATH_OBJ)pscoring.o \
		$(PATH_OBJ)utils.o $(PATH_OBJ)pertable.o $(PATH_OBJ)memhandler.o \
		$(PATH_OBJ)voronoi.o $(PATH_OBJ)sort.o $(PATH_OBJ)calc.o \
		$(PATH_OBJ)writepdb.o $(PATH_OBJ)rpdb.o $(PATH_OBJ)tparams.o \
		$(PATH_OBJ)fparams.o $(PATH_OBJ)pocket.o $(PATH_OBJ)refine.o \
		$(PATH_OBJ)descriptors.o $(PATH_OBJ)cluster.o $(PATH_OBJ)aa.o \
		$(PATH_OBJ)fpocket.o $(PATH_OBJ)write_visu.o  $(PATH_OBJ)fpout.o \
		$(PATH_OBJ)atom.o $(PATH_OBJ)writepocket.o $(PATH_OBJ)voronoi_lst.o \
		$(PATH_OBJ)neighbor.o $(PATH_OBJ)asa.o \
		$(QOBJS)

FPOBJ = $(PATH_OBJ)fpmain.o $(PATH_OBJ)psorting.o $(PATH_OBJ)pscoring.o \
		$(PATH_OBJ)utils.o $(PATH_OBJ)pertable.o $(PATH_OBJ)memhandler.o \
		$(PATH_OBJ)voronoi.o $(PATH_OBJ)sort.o $(PATH_OBJ)calc.o \
		$(PATH_OBJ)writepdb.o $(PATH_OBJ)rpdb.o $(PATH_OBJ)tparams.o \
		$(PATH_OBJ)fparams.o $(PATH_OBJ)pocket.o $(PATH_OBJ)refine.o \
		$(PATH_OBJ)descriptors.o $(PATH_OBJ)cluster.o $(PATH_OBJ)aa.o \
		$(PATH_OBJ)fpocket.o $(PATH_OBJ)write_visu.o  $(PATH_OBJ)fpout.o \
		$(PATH_OBJ)atom.o $(PATH_OBJ)writepocket.o $(PATH_OBJ)voronoi_lst.o $(PATH_OBJ)asa.o \
		$(QOBJS)

TPOBJ = $(PATH_OBJ)tpmain.o $(PATH_OBJ)psorting.o $(PATH_OBJ)pscoring.o \
		$(PATH_OBJ)utils.o $(PATH_OBJ)pertable.o $(PATH_OBJ)memhandler.o \
		$(PATH_OBJ)voronoi.o $(PATH_OBJ)sort.o $(PATH_OBJ)calc.o \
		$(PATH_OBJ)writepdb.o $(PATH_OBJ)rpdb.o $(PATH_OBJ)tparams.o \
		$(PATH_OBJ)fparams.o $(PATH_OBJ)pocket.o $(PATH_OBJ)refine.o \
		$(PATH_OBJ)tpocket.o  $(PATH_OBJ)descriptors.o $(PATH_OBJ)cluster.o \
		$(PATH_OBJ)aa.o $(PATH_OBJ)fpocket.o $(PATH_OBJ)write_visu.o \
		$(PATH_OBJ)fpout.o $(PATH_OBJ)atom.o $(PATH_OBJ)writepocket.o \
		$(PATH_OBJ)voronoi_lst.o $(PATH_OBJ)neighbor.o $(PATH_OBJ)asa.o\
		$(QOBJS)

DPOBJ = $(PATH_OBJ)dpmain.o $(PATH_OBJ)psorting.o $(PATH_OBJ)pscoring.o \
		$(PATH_OBJ)dpocket.o $(PATH_OBJ)dparams.o  $(PATH_OBJ)voronoi.o \
		$(PATH_OBJ)sort.o  $(PATH_OBJ)rpdb.o $(PATH_OBJ)descriptors.o \
		$(PATH_OBJ)neighbor.o $(PATH_OBJ)atom.o $(PATH_OBJ)aa.o \
		$(PATH_OBJ)pertable.o $(PATH_OBJ)calc.o $(PATH_OBJ)utils.o \
		$(PATH_OBJ)writepdb.o $(PATH_OBJ)memhandler.o $(PATH_OBJ)pocket.o \
		$(PATH_OBJ)refine.o $(PATH_OBJ)cluster.o $(PATH_OBJ)fparams.o \
		$(PATH_OBJ)fpocket.o $(PATH_OBJ)fpout.o $(PATH_OBJ)writepocket.o \
		$(PATH_OBJ)write_visu.o $(PATH_OBJ)asa.o\
		$(PATH_OBJ)voronoi_lst.o $(QOBJS)

MDPOBJ = $(PATH_OBJ)mdpmain.o $(PATH_OBJ)mdpocket.o $(PATH_OBJ)mdpbase.o $(PATH_OBJ)mdpout.o $(PATH_OBJ)psorting.o $(PATH_OBJ)pscoring.o \
		$(PATH_OBJ)mdparams.o $(PATH_OBJ)voronoi.o \
		$(PATH_OBJ)sort.o  $(PATH_OBJ)rpdb.o $(PATH_OBJ)descriptors.o \
		$(PATH_OBJ)neighbor.o $(PATH_OBJ)atom.o $(PATH_OBJ)aa.o \
		$(PATH_OBJ)pertable.o $(PATH_OBJ)calc.o $(PATH_OBJ)utils.o \
		$(PATH_OBJ)writepdb.o $(PATH_OBJ)memhandler.o $(PATH_OBJ)pocket.o \
		$(PATH_OBJ)refine.o $(PATH_OBJ)cluster.o $(PATH_OBJ)fparams.o \
		$(PATH_OBJ)fpocket.o $(PATH_OBJ)fpout.o \
		$(PATH_OBJ)writepocket.o $(PATH_OBJ)write_visu.o $(PATH_OBJ)asa.o \
		$(PATH_OBJ)voronoi_lst.o $(QOBJS)

#------------------------------------------------------------
# GENERAL RULES FOR COMPILATION
#------------------------------------------------------------

$(PATH_QHULL)%.o: $(PATH_QHULL)%.c
	$(CCQHULL) $(QCFLAGS) -c $< -o $@

$(PATH_OBJ)%.o: $(PATH_SRC)%.c
	$(CC) $(CFLAGS) -c $< -o $@
	
$(PATH_OBJ)%.o: $(PATH_SRC)%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
		
#-----------------------------------------------------------
# RULES FOR EXECUTABLES
#-----------------------------------------------------------

all: $(MYPROGS) $(PATH_BIN)$(CHECK)
		
$(PATH_BIN)$(CHECK): $(CHOBJ) $(QOBJS)
	$(LINKER) $(LFLAGS) $^ -o $@

$(PATH_BIN)$(FPOCKET): $(FPOBJ) $(QOBJS)
	$(LINKER) $(LFLAGS) $^ -o $@

$(PATH_BIN)$(TPOCKET): $(TPOBJ) $(QOBJS)
	$(LINKER) $(LFLAGS) $^ -o $@

$(PATH_BIN)$(DPOCKET): $(DPOBJ) $(QOBJS)
	$(LINKER) $(LFLAGS) $^ -o $@

$(PATH_BIN)$(MDPOCKET): $(MDPOBJ) $(QOBJS)
	$(LINKER) $(LFLAGS) $^ -o $@

install:
	mkdir -p $(BINDIR)
	mkdir -p $(MANDIR)
	cp $(PATH_BIN)$(FPOCKET) $(BINDIR)
	cp $(PATH_BIN)$(TPOCKET) $(BINDIR)
	cp $(PATH_BIN)$(DPOCKET) $(BINDIR)
	cp $(PATH_BIN)$(MDPOCKET) $(BINDIR)
	cp $(PATH_MAN)* $(MANDIR)

check:
	./$(PATH_BIN)$(CHECK)
		
test:
	./$(PATH_BIN)$(CHECK)

clean:
	rm -f $(PATH_QHULL)*.o
	rm -f $(PATH_OBJ)*.o

uninstall:
	rm -f $(PATH_BIN)$(FPOCKET) $(BINDIR)$(FPOCKET)
	rm -f $(PATH_BIN)$(TPOCKET) $(BINDIR)$(TPOCKET)
	rm -f $(PATH_BIN)$(DPOCKET) $(BINDIR)$(DPOCKET)
	rm -f $(PATH_BIN)$(MDPOCKET) $(BINDIR)$(MDPOCKET)
	rm -f $(MANDIR)fpocket.8 $(MANDIR)tpocket.8 $(MANDIR)dpocket.8
	
