  MEXSUFFIX  = mex
  DIR        = ../..
  CUDDHOME   = $(DIR)/cudd-2.5.0
  MEX        = $(MATLABHOME)/bin/mex
  CXX        = g++
  OBJFILES   = $(CUDDHOME)/util/datalimit.o $(CUDDHOME)/util/cpu_time.o $(CUDDHOME)/util/safe_mem.o $(DIR)/library/src/pessoa.o
  CFLAGS     = -L$(MATLABHOME)/bin/glnxa64 -L$(CUDDHOME)/dddmp -L$(CUDDHOME)/epd -L$(CUDDHOME)/util -L$(CUDDHOME)/st -L$(CUDDHOME)/mtr -L$(CUDDHOME)/cudd -L$(CUDDHOME)/include
  LIBS       = -lmex -lmx -lmat -lcudd -lmtr -lst -lutil -lepd -ldddmp
  INCLUDE    = -I$(MATLABHOME)/extern/include -I$(DIR)/cudd-2.5.0/include -I$(DIR)/library/src
  MEXFLAGS   = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'

  all: pessoa_design.$(MEXSUFFIX) 

        pessoa_design.$(MEXSUFFIX): pessoa_design.cpp 
		  $(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_design


