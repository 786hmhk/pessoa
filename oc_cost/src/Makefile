MEXSUFFIX  = mex
DIR        = ../..

MATLABHOME = /home/tanasaki/Programs/MATLAB/
CUDD_DIR   = ../../libs/cudd-2.5.0
SYSTEMEX   = glnx86

SP_DIR     = ../../shortestPath

MEX        = $(MATLABHOME)/bin/mex
CXX        = g++
CXXFLAGS   = -O0 -g3 -Wall -c -mtune=native -malign-double -DHAVE_IEEE_754 -DBSD -fPIC 
OBJFILES   = $(CUDD_DIR)/util/datalimit.o $(CUDD_DIR)/util/cpu_time.o $(CUDD_DIR)/util/safe_mem.o
CFLAGS     = -L$(MATLABHOME)/bin/$(SYSTEMEX) -L$(CUDD_DIR)/dddmp -L$(CUDD_DIR)/epd -L$(CUDD_DIR)/util -L$(CUDD_DIR)/st -L$(CUDD_DIR)/mtr -L$(CUDD_DIR)/cudd
LIBS       = -lmex -lmx -lmat -lcudd -lmtr -lst -lutil -lepd -ldddmp -lm
INCLUDE    = -I$(MATLABHOME)/extern/include -I$(CUDD_DIR)/include -I$(SP_DIR)/src
MEXFLAGS   = -cxx CC='$(CXX)' CXX='$(CXX)' LD='$(CXX)' CXXFLAGS='$(CXXFLAGS)'

all: shortestPathMex.$(MEXSUFFIX) pessoa_cost2add.$(MEXSUFFIX) pessoa_controller.$(MEXSUFFIX) pessoa_dsp.$(MEXSUFFIX)

shortestPathMex.$(MEXSUFFIX): shortestPathMex.cpp $(SP_DIR)/src/cuddObj.cc $(SP_DIR)/src/ShortestPath.cpp
	$(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/shortestPathMex
	
pessoa_cost2add.$(MEXSUFFIX): pessoa_cost2add.cpp $(SP_DIR)/src/cuddObj.cc
	$(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_cost2add

pessoa_dsp.$(MEXSUFFIX): pessoa_dsp.cpp $(SP_DIR)/src/cuddObj.cc $(SP_DIR)/src/ShortestPath.cpp
	$(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_dsp
	
pessoa_controller.$(MEXSUFFIX): pessoa_controller.cpp $(SP_DIR)/src/cuddObj.cc $(SP_DIR)/src/ShortestPath.cpp
	$(MEX) $(MEXFLAGS) $(CFLAGS) $(INCLUDE) $(LIBS) $(OBJFILES) $^ -output ../bin/pessoa_controller
	
run: shortestPathMex.$(MEXSUFFIX)
	matlab -nosplash -nodesktop -r tb_oc_sp
