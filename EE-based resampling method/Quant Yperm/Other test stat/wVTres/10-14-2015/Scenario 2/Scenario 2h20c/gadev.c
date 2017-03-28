/* Function to draw an observation from the standard normal distribution.
   This is from the book Numerical recipes in C, and uses the Box-Muller
   algorithm (gadev in the book).
*/


#include <math.h>
#include "rand1.h"

double gadev(long *idum) {
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	
	if (*idum < 0) {iset = 0;}
	
	if (iset == 0) {
		do {
			v1 = 2.0*rand1(idum)-1.0;
			v2 = 2.0*rand1(idum)-1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		
		fac = sqrt(-2.0*(log(rsq)/rsq));
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	} else {
		iset = 0;
		return gset;
	}
}

