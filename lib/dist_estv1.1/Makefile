.SUFFIXES: .c .f .o .a  .def .exp .dll .exe 

FFLAGS=-O3
CFLAGS=-O3
CC=gcc
F77=gfortran

SOURCES=$(wildcard -f *.c *.f)
OBJS=$(foreach i,$(SOURCES),$(basename $i).o)

dist_est: $(OBJS)
	$(F77) $(FFLAGS) -o dist_est $(OBJS) -lm

.c.o:
	$(CC)  $(CFLAGS) $($*-CFLAGS) -c $< -o $@ 

.f.o:
	$(F77) $(FFLAGS) $($*-FFLAGS) -c $< -o $@

clean :
	-rm *.o dist_est