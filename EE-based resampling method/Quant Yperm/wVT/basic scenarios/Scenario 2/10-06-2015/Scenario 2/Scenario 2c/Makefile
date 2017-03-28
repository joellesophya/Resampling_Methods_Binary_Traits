CC = gcc
CFLAGS = -O2
LDFLAGS =
LDLIBS = -lm -Lc:/lapack-3.5.0/build/lib/ -llapack 


OBJECTS = simulation.o nrutil.o rand1.o gadev.o cholesky.o svdcomp.o datacheck.o read.o sim_pheno_geno.o geno_sim.o trait_sim.o brent.o compute_residuals.o assotest.o permute_sim.o

all: simulation

.PHONY: CLAPACK_SVD

simulation: $(OBJECTS) 
	$(CC) $(LDFLAGS)  $(OBJECTS)  $(LDLIBS) -o Scen2VTc

debug: CFLAGS = -O0 -g 
debug: simulation

debug2: LDFLAGS = -pg 
debug2: CFLAGS = -O2 -pg 
debug2: simulation

clean:
	@rm -rf *.o gmon.out *.exe