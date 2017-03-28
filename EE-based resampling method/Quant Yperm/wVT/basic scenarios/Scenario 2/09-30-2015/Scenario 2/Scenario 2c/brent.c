#include "brent.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "compute_residuals.h"

/* brent's algorithm, modified from GSL root/brent.c
*/
typedef struct
{
	double a, b, c, d, e;
	double fa, fb, fc;
}
brent_state_t;

#define GSL_SUCCESS 1
#define GSL_CONTINUE -2
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define MAX_ITER 10000


static int
brent_iterate(brent_state_t * vstate, function_xi *f, double * root, double * x_lower, double * x_upper);
static int
brent_init(brent_state_t * vstate, function_xi *f, double * root, double x_lower, double x_upper);

int brent(function_xi *f, double *root) {
	int status;
	int iter = 0, max_iter = MAX_ITER;
	double r = 0;
	double x_lo = 0.0, x_hi = 1.0;
	double tol;
	double epsrel = GSL_SQRT_DBL_EPSILON;

	brent_state_t *current_state;
	current_state = malloc(sizeof(brent_state_t));
	brent_init(current_state, f, &r, x_lo, x_hi);
	

	/* printf ("\nusing %s method to estimate heritability:\n",
			 "brent");

			 printf ("%5s [%9s, %9s] %9s %9s\n",
			 "iter", "lower", "upper", "root","err(est)");*/

	do
	{
		iter++;
		status = brent_iterate(current_state, f, &r, &x_lo, &x_hi);
		tol = fmin(fabs(x_lo), fabs(x_hi))*epsrel;
		if (fabs(x_lo - x_hi) < tol) {
			// printf ("Converged:\n");
			status = GSL_SUCCESS;
			*root = r;
		}
		else status = GSL_CONTINUE;

		/*printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
				iter, x_lo, x_hi,
				r, x_hi - x_lo);*/
	} while (status == GSL_CONTINUE && iter < max_iter);
	//printf("\n");
	free(current_state);

	return status;
}

static int
brent_init(brent_state_t * vstate, function_xi *f, double * root, double x_lower, double x_upper)
{
	brent_state_t * state = vstate;

	double f_lower, f_upper;

	*root = 0.5 * (x_lower + x_upper);
	f_lower = f->func_of_xi(x_lower, f->parameters);
	f_upper = f->func_of_xi(x_upper, f->parameters);

	state->a = x_lower;
	state->fa = f_lower;

	state->b = x_upper;
	state->fb = f_upper;

	state->c = x_upper;
	state->fc = f_upper;

	state->d = x_upper - x_lower;
	state->e = x_upper - x_lower;

	if ((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
	{
		printf("endpoints do not straddle y=0");
	}

	return GSL_SUCCESS;

}


static int
brent_iterate(brent_state_t * vstate, function_xi *f, double * root, double * x_lower, double * x_upper)
{
	brent_state_t * state = vstate;

	double tol, m;

	int ac_equal = 0;

	double a = state->a, b = state->b, c = state->c;
	double fa = state->fa, fb = state->fb, fc = state->fc;
	double d = state->d, e = state->e;

	if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
	{
		ac_equal = 1;
		c = a;
		fc = fa;
		d = b - a;
		e = b - a;
	}

	if (fabs(fc) < fabs(fb))
	{
		ac_equal = 1;
		a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;
	}

	tol = 0.5 * GSL_DBL_EPSILON * fabs(b);
	m = 0.5 * (c - b);

	if (fb == 0)
	{
		*root = b;
		*x_lower = b;
		*x_upper = b;

		return GSL_SUCCESS;
	}

	if (fabs(m) <= tol)
	{
		*root = b;

		if (b < c)
		{
			*x_lower = b;
			*x_upper = c;
		}
		else
		{
			*x_lower = c;
			*x_upper = b;
		}

		return GSL_SUCCESS;
	}

	if (fabs(e) < tol || fabs(fa) <= fabs(fb))
	{
		d = m;		/* use bisection */
		e = m;
	}
	else
	{
		double p, q, r;	/* use inverse cubic interpolation */
		double s = fb / fa;

		if (ac_equal)
		{
			p = 2 * m * s;
			q = 1 - s;
		}
		else
		{
			q = fa / fc;
			r = fb / fc;
			p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
			q = (q - 1) * (r - 1) * (s - 1);
		}

		if (p > 0)
		{
			q = -q;
		}
		else
		{
			p = -p;
		}

		if (2 * p < fmin(3 * m * q - fabs(tol * q), fabs(e * q)))
		{
			e = d;
			d = p / q;
		}
		else
		{
			/* interpolation failed, fall back to bisection */

			d = m;
			e = m;
		}
	}

	a = b;
	fa = fb;

	if (fabs(d) > tol)
	{
		b += d;
	}
	else
	{
		b += (m > 0 ? +tol : -tol);
	}

	// SAFE_FUNC_CALL (f, b, &fb);
	fb = f->func_of_xi(b, f->parameters);

	state->a = a;
	state->b = b;
	state->c = c;
	state->d = d;
	state->e = e;
	state->fa = fa;
	state->fb = fb;
	state->fc = fc;

	/* Update the best estimate of the root and bounds on each
	   iteration */

	*root = b;

	if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
	{
		c = a;
	}

	if (b < c)
	{
		*x_lower = b;
		*x_upper = c;
	}
	else
	{
		*x_lower = c;
		*x_upper = b;
	}

	return GSL_SUCCESS;
}

