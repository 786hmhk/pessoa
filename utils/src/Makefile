MEXSUFFIX  = mex
DIR        = ../..
CUDDHOME   = $(DIR)/cudd-2.4.2
MEX        = $(MATLABHOME)/bin/mex
CXX        = g++
OBJFILES   = $(CUDDHOME)/util/datalimit.o $(CUDDHOME)/util/cpu_time.o $(CUDDHOME)/util/safe_mem.o $(DIR)/library/src/pessoa.o
CFLAGS     = -L$(MATLABHOME)/bin/glnxa64 -L$(CUDDHOME)/dddmp -L$(CUDDHOME)/epd -L$(CUDDHOME)/util -L$(CUDDHOME)/st -L$(CUDDHOME)/mtr -L$(CUDDHOME)/cudd -L$(CUDDHOME)/include
LIBS       = -lmex -lmx -lmat -lcudd -lmtr -lst -lutil -lepd -ldddmp
INCLUDE    = -I$(MATLABHOME)/extern/include -I$(DIR)/cudd-2.4.2/include -I$(DIR)/library/src
MEXFLAGS   = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'

all: pessoa_combine.$(MEXSUFFIX) pessoa_set2bdd.$(MEXSUFFIX) pessoa_charf2bdd.$(MEXSUFFIX) pessoa_compare_sets.$(MEXSUFFIX)
 
pessoa_combine.$(MEXSUFFIX): pessoa_combine.cpp 
		  $(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_combine

pessoa_set2bdd.$(MEXSUFFIX): pessoa_set2bdd.cpp 
		  $(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_set2bdd

pessoa_charf2bdd.$(MEXSUFFIX): pessoa_charf2bdd.cpp 
		  $(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_charf2bdd

pessoa_compare_sets.$(MEXSUFFIX): pessoa_compare_sets.cpp 
		  $(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_compare_sets