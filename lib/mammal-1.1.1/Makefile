.SUFFIXES: .c .f .o .a  .def .exp .dll .exe 

FFLAGS=-O3
CFLAGS=-O3
CC=gcc
F77=gfortran

SOURCES=$(wildcard -f *.c *.f)
OBJS=$(foreach i,$(SOURCES),$(basename $i).o)

all: mammal-sigma mult-data mult-mix-lwt charfreq dgpe

mammal-sigma: mammal-sigma.o mammal-sigmaf.o
	$(CC) $(CFLAGS) -o mammal-sigma mammal-sigma.o mammal-sigmaf.o -lm

mult-data: mult-data.o mult-dataf.o
	$(CC) $(CFLAGS) -o mult-data mult-data.o mult-dataf.o -lm

mult-mix-lwt: mult-mix-lwt.o mult-mix-lwtf.o
	$(CC) $(CFLAGS) -o mult-mix-lwt mult-mix-lwt.o mult-mix-lwtf.o -lm

charfreq: charfreq.o charfreqf.o
	$(CC) $(CFLAGS) -o charfreq charfreq.o charfreqf.o -lm

dgpe: dgpe.o dgpef.o dcdflib.o ipmpar.o
	$(CC) $(CFLAGS) -o dgpe dgpe.o dgpef.o dcdflib.o ipmpar.o -lm

.c.o:
	$(CC)  $(CFLAGS) $($*-CFLAGS) -c $< -o $@ 

.f.o:
	$(F77) $(FFLAGS) $($*-FFLAGS) -c $< -o $@
