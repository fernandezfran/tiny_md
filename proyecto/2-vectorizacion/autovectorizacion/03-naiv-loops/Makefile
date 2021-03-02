include in.parameters

CC       = gcc
CFLAGS	 = -O1 -funroll-loops -ffast-math -march=native -ftree-vectorize -fopt-info-vec-optimized
WFLAGS	 = -std=gnu99 -Wall -Wextra -g
LDFLAGS	 = -lm -lgomp

PARAMETERS	= -DN=$(N) -DSEED=$(SEED) -DT0=$(T0) -DRhoi=$(Rhoi) -Drcut=$(rcut) \
			  -Ddt=$(dt) -Dteq=$(teq) -Dtrun=$(trun) -Dtmes=$(tmes) \

TARGETS		= tiny_md
SOURCES		= $(shell echo *.c)
OBJECTS     = 

all: $(TARGETS)

tiny_md: tiny_md.o
	$(CC) $(WFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tiny_md.o: tiny_md.c
	$(CC) $(WFLAGS) $(CFLAGS) $(PARAMETERS) -c $< 

clean:
	rm -f $(TARGETS) *.o *.xyz *.log

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all
