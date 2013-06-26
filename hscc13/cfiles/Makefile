
include ../include.mk

OBJ = mcis.o $(CUDDHOME)/util/datalimit.o $(CUDDHOME)/util/cpu_time.o $(CUDDHOME)/util/safe_mem.o

all: mcis.$(SO_SUFFIX)


%.o: %.c
	$(CC) $(CFLAGS) $(HDIR) -c $^ -o $@


%.$(SO_SUFFIX): $(OBJ)
	$(CC) $(LDFLAGS) $(LIBDIR) $(OBJ) -o $@ $(LIB)
	mv $@ $(MDIR)


clean:
	rm -rf *.o *.$(SO_SUFFIX)

