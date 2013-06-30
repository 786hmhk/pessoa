# compiler
CC = gcc

# make
MAKE = make

# code root directory
ROOTDIR = /home/rungger/code/hscc13

# mac or *nix suffix 
SO_SUFFIX  = dylib #so


# dirs
DIRS = cfiles cudd-2.5.0

# matlab 
#MHOME = /usr/lib/matlab-7.13/
MHOME = /mnt/matlab
MHOME = /Applications/MATLAB_R2012a.app

MARCH = maci64 # glnxa64 glnx86


# mfiles
MDIR = $(ROOTDIR)/mfiles

# cudd library
CUDDHOME = $(ROOTDIR)/cudd-2.5.0

# header files
HDIR  = -I$(CUDDHOME)/include -I$(MHOME)/extern/include -I$(ROOTDIR)/hfiles


CFLAGS =  -ansi -ggdb -Wall -Wextra -fPIC -DPIC  -DHAVE_IEEE_754 -DBSD -DSIZEOF_VOID_P=8 -DSIZEOF_LONG=8  


# libraries
LIBDIR = -L$(CUDDHOME)/dddmp -L$(CUDDHOME)/epd -L$(CUDDHOME)/util -L$(CUDDHOME)/st -L$(CUDDHOME)/mtr -L$(CUDDHOME)/cudd -L$(MHOME)/bin/$(MARCH)

LIB = -lcudd -lmtr -lst -lutil -lepd -ldddmp -ldl -lm -lmex

LDFLAGS  = -shared -m64 


