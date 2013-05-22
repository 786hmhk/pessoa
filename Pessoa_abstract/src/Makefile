  MEXSUFFIX  = mex
  DIR        = ../..
  CUDDHOME   = $(CUDDDIR)
  MEX        = $(MATLABHOME)/bin/mex
  CXX        = g++
  OBJFILES   = $(CUDDHOME)/util/datalimit.o $(CUDDHOME)/util/cpu_time.o $(CUDDHOME)/util/safe_mem.o $(DIR)/library/src/pessoa.o
  CFLAGS     = -L$(MATLABHOME)/bin/$(SYSTEMEX) -L$(CUDDHOME)/dddmp -L$(CUDDHOME)/epd -L$(CUDDHOME)/util -L$(CUDDHOME)/st -L$(CUDDHOME)/mtr -L$(CUDDHOME)/cudd -L$(CUDDHOME)/include
  XCFLAGS    = -DHAVE_IEEE_754 -DBSD -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8
  LIBS       = -lmex -lmx -lmat -lcudd -lmtr -lst -lutil -lepd -ldddmp
  INCLUDE    = -I$(MATLABHOME)/extern/include -I$(CUDDHOME)/include -I$(DIR)/library/src 
  MEXFLAGS   = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)'

  all: pessoa_abstract.$(MEXSUFFIX) pessoa_abstract_add.$(MEXSUFFIX) 

  pessoa_abstract.$(MEXSUFFIX): pessoa_abstract.cpp 
	  $(MEX) $(MEXFLAGS) $(CFLAGS) $(XCFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_abstract

  pessoa_abstract_add.$(MEXSUFFIX): pessoa_abstract_add.cpp 
	  $(MEX) $(MEXFLAGS) $(CFLAGS) $(XCFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_abstract_add

