CC          := mpicc
LFLAGS      := -O3 -std=c99
ALL         := gpm-par.exe

all : $(ALL)


gpm-par.exe :
	$(CC) $(LFLAGS) -o $@ utils.c graph.c gpm-par.c

clean :
	rm -f *.o *.out *.err $(ALL)