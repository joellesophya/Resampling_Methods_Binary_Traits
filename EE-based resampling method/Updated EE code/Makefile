CC = R CMD SHLIB
CFLAGS = -O2
LDFLAGS =
LDLIBS = -lm -L/usr/local/lib -llapack 


OBJECTS = simulation.c nrutil.c rand1.c gadev.c cholesky.c svdcomp.c datacheck.c read.c permute_sim.c sim_pheno_geno.c geno_sim.c assotest.c

all: simulation

simulation: $(OBJECTS) 
	$(CC) $(LDFLAGS)  $(OBJECTS)  $(LDLIBS) -o EE_Ceramic

debug: CFLAGS = -O0 -g 
debug: simulation

debug2: LDFLAGS = -pg 
debug2: CFLAGS = -O2 -pg 
debug2: simulation

clean:
	@rm -rf *.o gmon.out *.exe