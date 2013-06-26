
include include.mk

.PHONY: all cfiles

all: cfiles

cfiles:
	$(MAKE) -C $@ 

clean: 
	rm -rf cfiles/*.o 
	rm -rf mfiles/*.$(SO_SUFFIX)




