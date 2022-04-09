#-----------------------------
# generated from Makefile: DO NOT EDIT
# -----------------------------
SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libgurobi.lib

GUROBI=/home/jpc/gurobi951/linux64
GUROBILIB=-lgurobi95

GUROBI_INC= -I $(GUROBI)/include 
GUROBI_LIB= -L$(GUROBI)/lib -Wl,-R$(GUROBI)/lib $(GUROBILIB)

OBJS= nspgurobi-IN.obj nspgurobi.obj 

include $(SCIDIR)/Makefile.incl.mak

CFLAGS = $(CC_OPTIONS) $(GUROBI_INC)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(GUROBI_LIB)

include $(SCIDIR)/config/Makeso.incl



Makefile.mak	: Makefile
	$(SCIDIR)/scripts/Mak2VCMak Makefile

Makefile.libmk	: Makefile
	$(SCIDIR)/scripts/Mak2ABSMak Makefile

distclean:: clean 

clean:: 
	@$(RM) *.obj *.lo 
	@$(RM) -r */.libs 

