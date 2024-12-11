SHELL = /bin/sh

include ../Path.incl

LIBRARY=libgurobi.a

GUROBI=/home/jpc/gurobi1200/linux64
GUROBILIB=-lgurobi120

GUROBI_INC= -I $(GUROBI)/include 
GUROBI_LIB= -L$(GUROBI)/lib -Wl,-R$(GUROBI)/lib $(GUROBILIB)

OBJS= nspgurobi-IN.o nspgurobi.o 

include $(SCIDIR)/Makefile.incl

CFLAGS = $(CC_OPTIONS) $(GUROBI_INC)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(GUROBI_LIB)

include $(SCIDIR)/config/Makeso.incl

all:: Makefile.mak

distclean:: clean 

clean:: 
	@$(RM) *.o *.lo 
	@$(RM) -r */.libs 

