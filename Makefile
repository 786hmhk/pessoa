#---------------------------------------------------------------------------#
# Makefile for Pessoa							     #
# University of California, Los Angeles                                      #
# Version 1.0, June, 2009	                                     	     #
#----------------------------------------------------------------------------#

MEXSUFFIX  = mex
DIR        = ../..
MATLAB	   = /opt/netapps/MathWorks/amd64/matlab/mdw/R2010a/
 #/Applications/MATLAB74/
CUDDDIR	   = ./cudd-2.4.2/
SYSTEMEX   = glnxa64 #mexmaci #glnxa64 #mexmaci
LDFLAGS_LIB    = -arch x86_64 #-mtune=native -fPIC -arch i386 #-arch x86_64
MAKEARG    = 'MATLABHOME=$(MATLAB)' 'CUDDDIR=../../$(CUDDDIR)' 'SYSTEMEX=$(SYSTEMEX)' 'LDFLAGS_LIB=$(LDFLAGS_LIB)'

.PHONY: all clean utils cudd

all: cudd pessoa.o pessoa_abstract pessoa_design pessoa_controller utils

clean:
	$(MAKE) distclean -C $(CUDDDIR)
	rm -f ./library/src/pessoa.o ./utils/bin/*.$(MEXSUFFIX)* ./Pessoa_abstract/bin/*.$(MEXSUFFIX)* ./Pessoa_design/bin/*.$(MEXSUFFIX)* ./Pessoa_control/bin/*.$(MEXSUFFIX)*

cudd:
	$(MAKE) -C $(CUDDDIR) 

pessoa.o: cudd
	$(MAKE) $(MAKEARG) -C ./library/src   

pessoa_abstract: pessoa.o 
	$(MAKE) $(MAKEARG) -C ./Pessoa_abstract/src   

pessoa_design: pessoa.o 
	$(MAKE) $(MAKEARG) -C ./Pessoa_design/src   

pessoa_controller: pessoa.o
	$(MAKE) $(MAKEARG) -C ./Pessoa_control/src   

utils: pessoa.o
	$(MAKE) $(MAKEARG) -C ./utils/src
