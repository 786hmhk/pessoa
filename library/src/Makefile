#----------------------------------------------------------------------------#
# Makefile for Pessoa							     #
# University of California, Los Angeles                                      #
# Version 1.0, June, 2009	                                     	     #
#----------------------------------------------------------------------------#

  MEXSUFFIX  = mex
  DIR        = ../..
  CUDDHOME   = $(DIR)/cudd-2.4.2
  CXX        = g++
  CFLAGS     = -fPIC -L$(MATLABHOME)/bin/glnxa64 -L$(CUDDHOME)/dddmp -L$(CUDDHOME)/epd -L$(CUDDHOME)/util -L$(CUDDHOME)/st -L$(CUDDHOME)/mtr -L$(CUDDHOME)/cudd -L$(CUDDHOME)/include
  INCLUDE    = -I$(MATLABHOME)/extern/include -I$(DIR)/cudd-2.5.0/include

  pessoa.o: pessoa.cpp pessoa.h 
	 $(CXX) $(CFLAGS) $(INCLUDE) -c pessoa.cpp
