/* Generate a random number between 0 and 1, noninclusive
rand1 of numerical recipes*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMAX (1.0-EPS)

#include "rand1.h"

double rand1(idum)
long *idum;
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy){
		if (-(*idum) < 1) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--){
			k = (*idum) / IQ;
			*idum = IA*(*idum - k*IQ) - IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA*(*idum - k*IQ) - IR*k;
	if (*idum <0) *idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM*iy) > RNMAX) temp = RNMAX;
	return temp;
}



