include in.parameters

CC       = gcc
CFLAGS	 = -O1 -funroll-loops -ffast-math -march=native -ftree-vectorize -fopenmp
WFLAGS	 = -std=gnu99 -Wall -Wextra -g
LDFLAGS	 = -lm -lgomp

PARAMETERS	= -DN=$(N) -DSEED=$(SEED) -DT0=$(T0) -DRhoi=$(Rhoi) -Drcut=$(rcut) \
			  -Ddt=$(dt) -Dteq=$(teq) -Dtrun=$(trun) -Dtmes=$(tmes) \

TARGETS		= tiny_md viz
SOURCES		= $(shell echo *.c)
OBJECTS     = init.o lj_force.o verlet.o pbc.o rescale.o md_step.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(WFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(WFLAGS) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(WFLAGS) $(CFLAGS) $(PARAMETERS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all
