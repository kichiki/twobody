/* coef of 2-body res functions on Jeffrey & Onishi 1984 JFM vol.139 p.261
 * using GMP library
 * Copyright (C) 1999 Kengo ICHIKI (kengo@caltech.edu)
 * $Id: two-body-JO-gmp.c,v 3.1 1999/08/27 22:38:37 ichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h> /* malloc() free() */
#include "gmp.h"

#include "../bench.h" /* ptime_ms_d() */


/* function prototypes */
int
main (void);

void
XA_sum_l (int nmax, mpq_t lambda);
void
YA_sum_l (int nmax, mpq_t lambda);
void
YB_sum_l (int nmax, mpq_t lambda);
void
XC_sum_l (int nmax, mpq_t lambda);
void
YC_sum_l (int nmax, mpq_t lambda);
void
XG_sum_l (int nmax, mpq_t lambda);
void
YG_sum_l (int nmax, mpq_t lambda);
void
YH_sum_l (int nmax, mpq_t lambda);
void
XM_sum_l (int nmax, mpq_t lambda);
void
YM_sum_l (int nmax, mpq_t lambda);
void
ZM_sum_l (int nmax, mpq_t lambda);

void
XP_sum_l (int nmax, mpq_t lambda);
void
XQ_sum_l (int nmax, mpq_t lambda);


void
XC (int nmax);

void
X_p (int n, int p, int q, int inp, int inv, mpq_t coef_p);
void
X_v (int n, int p, int q, int inp, int inv, mpq_t coef_v);
void
X_q (int n, int p, int q, mpq_t coef_q);

void
Y_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_p);
void
Y_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_v);
void
Y_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);

void
YM_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);
void
YM_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);
void
YM_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);

void
ZM_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);
void
ZM_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);
void
ZM_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q);

int
comb (mpq_t comb, int n, int m);

void
print_mpq (mpq_t x);
void
fprint_mpq (FILE *out, mpq_t x);

int
search_results (char *file, int n, int p, int q, mpq_t coef);
void
append_result (char *file, int n, int p, int q, mpq_t coef);


int
main (void)
{
  int nmax;
  mpq_t lambda;
  double time;


  nmax = 100;

  mpq_init (lambda);
  mpq_set_ui (lambda, 1, 1);

  /*XC (nmax);*/

  ptime_ms_d();

  XA_sum_l (nmax, lambda);
  YA_sum_l (nmax, lambda);
  YB_sum_l (nmax, lambda);
  XC_sum_l (nmax, lambda);
  YC_sum_l (nmax, lambda);
  XG_sum_l (nmax, lambda);
  YG_sum_l (nmax, lambda);
  YH_sum_l (nmax, lambda);
  XM_sum_l (nmax, lambda);
  YM_sum_l (nmax, lambda);
  ZM_sum_l (nmax, lambda);
  XP_sum_l (nmax, lambda);
  XQ_sum_l (nmax, lambda);

  time = ptime_ms_d();

  fprintf (stdout, "# ptime = %f\n", time);

  mpq_clear (lambda);

  return 0;
}


void
XA_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  X_p (1, k - q, q, 1, 1, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	/* mul 2^k */
	mpq_mul (f, f, two_k);

      fprintf (stdout, "XAf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YA_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  Y_p (1, k - q, q, 0, 1, 1, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	/* mul 2^k */
	mpq_mul (f, f, two_k);

      fprintf (stdout, "YAf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YB_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 0, 1, 1, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	      mpq_add (f, f, coef_q);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	  /* mul 2 */
	  mpq_set_ui (coef_q, 2, 1);
	  mpq_mul (f, f, coef_q);
	}

      fprintf (stdout, "YBf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XC_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  X_q (1, k - q, q, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);
	      mpq_add (f, f, coef_q);
	    }
	}
      fprintf (stdout, "XCf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YC_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 1, 0, 0, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	      mpq_add (f, f, coef_q);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	/* mul 2^k */
	mpq_mul (f, f, two_k);

      fprintf (stdout, "YCf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XG_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 1, 1, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	  /* mul 3 */
	  mpq_set_ui (coef_p, 3, 1);
	  mpq_mul (f, f, coef_p);
	  /* mul 4 */
	  mpq_set_ui (coef_p, 4, 1);
	  mpq_div (f, f, coef_p);
	}

      fprintf (stdout, "XGf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YG_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 0, 1, 1, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	  /* mul 3 */
	  mpq_set_ui (coef_p, 3, 1);
	  mpq_mul (f, f, coef_p);
	  /* div 4 */
	  mpq_set_ui (coef_p, 4, 1);
	  mpq_div (f, f, coef_p);
	}

      fprintf (stdout, "YGf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YH_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 1, 0, 0, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	  /* mul -3 */
	  mpq_set_si (coef_p, -3, 1);
	  mpq_mul (f, f, coef_p);
	  /* div 8 */
	  mpq_set_ui (coef_p, 8, 1);
	  mpq_div (f, f, coef_p);
	}

      fprintf (stdout, "YHf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XM_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	}

      fprintf (stdout, "XMf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YM_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  YM_p (2, k - q, q, 0, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	}

      fprintf (stdout, "YMf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ZM_sum_l (int nmax, mpq_t lambda)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 0; q <= k; q ++)
	{
	  ZM_p (2, k - q, q, 0, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	      mpq_add (f, f, coef_p);
	}
      if (mpz_cmp_si (mpq_numref (f), 0))
	{
	  /* mul 2^k */
	  mpq_mul (f, f, two_k);
	}

      fprintf (stdout, "ZMf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XP_sum_l (int nmax, mpq_t lambda)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 1; q <= k - 1; q ++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      X_p (n, q - n, k - q - 1, 1, 1, b);

	      if (mpz_cmp_si (mpq_numref (b), 0))
		{
		  /* mul (two_k) */
		  mpq_mul (b, b, two_k);
		  /* mul 3 */
		  mpq_set_ui (tmp, 3, 1);
		  mpq_mul (b, b, tmp);
		  /* mul 2 */
		  mpq_set_ui (tmp, 2, 1);
		  mpq_div (b, b, tmp);

		  mpq_add (a, a, b);
		}
	    }
	  if (mpz_cmp_si (mpq_numref (a), 0))
	      mpq_add (f, f, a);
	}

      fprintf (stdout, "XPf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XQ_sum_l (int nmax, mpq_t lambda)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      mpq_set_ui (f, 0, 1);
      for (q = 1; q <= k - 1; q ++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      X_p (n, q - n, k - q - 1, 2, 2, b);

	      if (mpz_cmp_si (mpq_numref (b), 0))
		{
		  /* mul (two_k) */
		  mpq_mul (b, b, two_k);
		  /* mul 5 */
		  mpq_set_ui (tmp, 5, 1);
		  mpq_mul (b, b, tmp);
		  /* mul 2 */
		  mpq_set_ui (tmp, 2, 1);
		  mpq_div (b, b, tmp);

		  mpq_add (a, a, b);
		}
	    }
	  if (mpz_cmp_si (mpq_numref (a), 0))
	      mpq_add (f, f, a);
	}

      fprintf (stdout, "XPf [%d] = ", k);
      print_mpq (f);
      fprintf (stdout, "\n");
    }

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}


void
XC (int nmax)
{
  int k, q;
  mpq_t two;
  mpq_t two_k;
  mpq_t coef_q;


  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      fprintf (stdout, "XCf_%d = ", k);
      for (q = 0; q <= k; q ++)
	{
	  X_q (1, k - q, q, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);

	      print_mpq (coef_q);
	      fprintf (stdout, " l^%d ", q);
	    }
	}
      fprintf (stdout, "\n");
    }

  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}



/* calc coef Pnpq on X[AGP]
 * INPUT
 *   n, p, q : (int)
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
X_p (int n, int p, int q, int inp, int inv, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Pnpq */
      mpq_set_ui (coef_p, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inp != 0
	  && n == inp)
	mpq_set_ui (coef_p, 1, 1);
      else
	mpq_set_ui (coef_p, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.X_p.%d.%d.%d", 0, inp, inv);
  if (search_results (filename, n, p, q, coef_p) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      X_p (s, q - s, p - n + 1, inp, inv, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2n+1) */
	  mpq_set_si (tmp, 2 * n + 1, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2ns-n-s+2) */
	  mpq_set_si (tmp, 2 * n * s - n - s + 2, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2(2s-1) */
	  mpq_set_si (tmp, 2 * (2 * s - 1), 1);
	  mpq_div (b, b, tmp);
	  /* div (n+s) */
	  mpq_set_si (tmp, n + s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
      X_p (s, q - s, p - n - 1, inp, inv, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -n */
	  mpq_set_si (tmp, - n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2n-1) */
	  mpq_set_si (tmp, 2 * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
      X_v (s, q - s - 2, p - n + 1, inp, inv, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -n */
	  mpq_set_si (tmp, - n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (4n^2-1) */
	  mpq_set_si (tmp, 4 * n * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2(2s+1) */
	  mpq_set_si (tmp, 2 * (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n) */
	  comb (tmp, n + s, n);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_p, coef_p, a);
	}
    }

  append_result (filename, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Vnpq on X[AGP]
 * INPUT
 *   n, p, q : (int)
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
X_v (int n, int p, int q, int inp, int inv, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_v, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inv != 0
	  && n == inv)
	mpq_set_ui (coef_v, 1, 1);
      else
	mpq_set_ui (coef_v, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.X_v.%d.%d.%d", 0, inp, inv);
  if (search_results (filename, n, p, q, coef_v) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  X_p (n, p, q, inp, inv, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      X_p (s, q - s, p - n - 1, inp, inv, a);

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n) */
	  comb (tmp, n + s, n);
	  mpq_mul (a, a, tmp);
	  /* mul (-2n) */
	  mpq_set_si (tmp, - 2 * n, 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);
	  /* div (2n+3) */
	  mpq_set_si (tmp, 2 * n + 3, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_v, coef_v, a);
	}
    }

  append_result (filename, n, p, q, coef_v);

  mpq_clear (a);
  mpq_clear (tmp);
}

/* calc coef Qnpq on XC
 * INPUT
 *   n, p, q : (int)
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
X_q (int n, int p, int q, mpq_t coef_q)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Qnpq */
      mpq_set_ui (coef_q, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (n == 1)
	mpq_set_ui (coef_q, 1, 1);
      else
	mpq_set_ui (coef_q, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.X_q.%d.%d.%d", 1, 0, 0);
  if (search_results (filename, n, p, q, coef_q) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (tmp);

  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);
  for (s = 0; s <= q; s ++)
    {
      /* Qnpq : Qs(q-s-1)(p-n) */
      X_q (s, q - s - 1, p - n, a);

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n) */
	  comb (tmp, n + s, n);
	  mpq_mul (a, a, tmp);
	  /* mul (s) */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  append_result (filename, n, p, q, coef_q);

  mpq_clear (a);
  mpq_clear (tmp);
}


/* calc coef Pnpq on Y[ABCGH]
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
Y_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Pnpq */
      mpq_set_ui (coef_p, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inp != 0
	  && n == inp)
	mpq_set_ui (coef_p, 1, 1);
      else
	mpq_set_ui (coef_p, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.Y_p.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_p) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      Y_p (s, q - s, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (2n+1) */
	  mpq_set_si (tmp, 2 * n + 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  /* mul 3(n+s)-(ns+1)(2ns-s-n+2) */
	  /* calc -(ns+1)(2ns-s-n+2) */
	  mpq_set_si (tmp0, - (n * s + 1), 1);
	  mpq_set_si (tmp, 2 * n * s - s - n + 2, 1);
	  mpq_mul (tmp, tmp, tmp0);
	  /* add 3(n+s) */
	  mpq_set_si (tmp0, 3 * (n + s), 1);
	  mpq_add (tmp, tmp, tmp0);
	  mpq_mul (b, b, tmp);

	  /* div (s) */
	  mpq_set_si (tmp, s, 1);
	  mpq_div (b, b, tmp);
	  /* div (n+s) */
	  mpq_set_si (tmp, n + s, 1);
	  mpq_div (b, b, tmp);
	  /* div (2s-1) */
	  mpq_set_si (tmp, 2 * s - 1, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
      Y_p (s, q - s, p - n - 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2n-1) */
	  mpq_set_si (tmp, 2 * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
      Y_v (s, q - s - 2, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (4n^2-1) */
	  mpq_set_si (tmp, 4 * n * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2(2s+1) */
	  mpq_set_si (tmp, 2 * (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
      Y_q (s, q - s - 1, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -2 */
	  mpq_set_si (tmp, - 2, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (4n^2-1) */
	  mpq_set_si (tmp, 4 * n * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 3 */
	  mpq_set_si (tmp, 3, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_p, coef_p, a);
	}
    }

  append_result (filename, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
}

/* calc coef Vnpq on Y[ABCGH]
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
Y_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_v, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inv != 0
	  && n == inv)
	mpq_set_ui (coef_v, 1, 1);
      else
	mpq_set_ui (coef_v, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.Y_v.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_v) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  Y_p (n, p, q, inq, inp, inv, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n+1) */
      Y_p (s, q - s, p - n - 1, inq, inp, inv, a);
      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* mul (2n) */
	  mpq_set_si (tmp, 2 * n, 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);
	  /* div (2n+3) */
	  mpq_set_si (tmp, 2 * n + 3, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_v, coef_v, a);
	}
    }

  append_result (filename, n, p, q, coef_v);

  mpq_clear (a);
  mpq_clear (tmp);
}

/* calc coef Qnpq on Y[ABCGH]
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
Y_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_q, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inq != 0
	  && n == inq)
	mpq_set_ui (coef_q, 1, 1);
      else
	mpq_set_ui (coef_q, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.Y_q.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_q) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Qnpq */
      mpq_set_ui (a, 0, 1);
      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
      Y_q (s, q - s - 1, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul s */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-n) */
      Y_p (s, q - s, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -3 */
	  mpq_set_si (tmp, - 3, 1);
	  mpq_mul (b, b, tmp);

	  /* div (2ns) */
	  mpq_set_si (tmp, 2 * n * s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  append_result (filename, n, p, q, coef_q);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Pnpq on YM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
YM_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0, tmp1;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Pnpq */
      mpq_set_ui (coef_p, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inp != 0
	  && n == inp)
	mpq_set_ui (coef_p, 1, 1);
      else
	mpq_set_ui (coef_p, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.YM_p.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_p) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      YM_p (s, q - s, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (2n+1) */
	  mpq_set_si (tmp, 2 * n + 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  /* mul (ns+4)(n+s)-2(ns+1)^2 */
	  /* calc -2(ns+1)^2 */
	  mpq_set_si (tmp, (n * s + 1), 1);
	  mpq_set_si (tmp0, (n * s + 1), 1);
	  mpq_mul (tmp, tmp, tmp0);
	  mpq_set_si (tmp0, - 2, 1);
	  mpq_mul (tmp, tmp, tmp0);
	  /* add (ns+4)(n+s) */
	  mpq_set_si (tmp0, n + s, 1);
	  mpq_set_si (tmp1, n * s + 4, 1);
	  mpq_mul (tmp0, tmp0, tmp1);
	  mpq_add (tmp, tmp, tmp0);
	  /* mul it with b */
	  mpq_mul (b, b, tmp);

	  /* div (s) */
	  mpq_set_si (tmp, s, 1);
	  mpq_div (b, b, tmp);
	  /* div (2s-1) */
	  mpq_set_si (tmp, 2 * s - 1, 1);
	  mpq_div (b, b, tmp);
	  /* div (n+s) */
	  mpq_set_si (tmp, n + s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
      YM_p (s, q - s, p - n - 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2n-1) */
	  mpq_set_si (tmp, 2 * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
      YM_v (s, q - s - 2, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (4n^2-1) */
	  mpq_set_si (tmp, 4 * n * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2(2s+1) */
	  mpq_set_si (tmp, 2 * (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
      YM_q (s, q - s - 1, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -(4n^2-1) */
	  mpq_set_si (tmp, - (4 * n * n - 1), 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_p, coef_p, a);
	}
    }

  append_result (filename, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
}

/* calc coef Vnpq on YM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
YM_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_v, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inv != 0
	  && n == inv)
	mpq_set_ui (coef_v, 1, 1);
      else
	mpq_set_ui (coef_v, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.YM_v.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_v) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  YM_p (n, p, q, inq, inp, inv, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      YM_p (s, q - s, p - n - 1, inq, inp, inv, a);
      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* mul (2n) */
	  mpq_set_si (tmp, 2 * n, 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);
	  /* div (2n+3) */
	  mpq_set_si (tmp, 2 * n + 3, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_v, coef_v, a);
	}
    }

  append_result (filename, n, p, q, coef_v);

  mpq_clear (a);
  mpq_clear (tmp);
}

/* calc coef Qnpq on YM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
YM_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_q, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inq != 0
	  && n == inq)
	mpq_set_ui (coef_q, 1, 1);
      else
	mpq_set_ui (coef_q, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.YM_q.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_q) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Qnpq */
      mpq_set_ui (a, 0, 1);
      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
      YM_q (s, q - s - 1, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul s */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-n) */
      YM_p (s, q - s, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* div -(ns) */
	  mpq_set_si (tmp, - n * s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+1) */
	  comb (tmp, n + s, n + 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  append_result (filename, n, p, q, coef_q);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Pnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
ZM_p (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0, tmp1;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Pnpq */
      mpq_set_ui (coef_p, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inp != 0
	  && n == inp)
	mpq_set_ui (coef_p, 1, 1);
      else
	mpq_set_ui (coef_p, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.ZM_p.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_p) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      ZM_p (s, q - s, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  mpq_set_ui (tmp1, 1, 1);
	  /* calc (ns+16)(n+s) */
	  mpq_set_si (tmp, n * s + 16, 1);
	  mpq_mul (tmp1, tmp1, tmp);
	  mpq_set_si (tmp, n + s, 1);
	  mpq_mul (tmp1, tmp1, tmp);

	  /* calc -2(ns+1)(ns+4) */
	  mpq_set_si (tmp, -2, 1);
	  mpq_set_si (tmp0, n * s + 1, 1);
	  mpq_mul (tmp, tmp, tmp0);
	  mpq_set_si (tmp0, n * s + 4, 1);
	  mpq_mul (tmp, tmp, tmp0);
	  /* calc (ns+16)(n+s)-2(ns+1)(ns+4) */
	  mpq_add (tmp1, tmp1, tmp);

	  mpq_mul (b, b, tmp1);

	  /* mul (2n+1) */
	  mpq_set_si (tmp, 2 * n + 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  /* div (s) */
	  mpq_set_si (tmp, s, 1);
	  mpq_div (b, b, tmp);
	  /* div (2s-1) */
	  mpq_set_si (tmp, 2 * s - 1, 1);
	  mpq_div (b, b, tmp);
	  /* div (n+s) */
	  mpq_set_si (tmp, n + s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
      ZM_p (s, q - s, p - n - 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (2n-1) */
	  mpq_set_si (tmp, 2 * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2 */
	  mpq_set_si (tmp, 2, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
      ZM_v (s, q - s - 2, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul n */
	  mpq_set_si (tmp, n, 1);
	  mpq_mul (b, b, tmp);
	  /* mul (4n^2-1) */
	  mpq_set_si (tmp, 4 * n * n - 1, 1);
	  mpq_mul (b, b, tmp);
	  /* div 2(2s+1) */
	  mpq_set_si (tmp, 2 * (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
      ZM_q (s, q - s - 1, p - n + 1, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -2(4n^2-1) */
	  mpq_set_si (tmp, - 2 * (4 * n * n - 1), 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+2) */
	  comb (tmp, n + s, n + 2);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_p, coef_p, a);
	}
    }

  append_result (filename, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
}

/* calc coef Vnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
ZM_v (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_v, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inv != 0
	  && n == inv)
	mpq_set_ui (coef_v, 1, 1);
      else
	mpq_set_ui (coef_v, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.ZM_v.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_v) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  ZM_p (n, p, q, inq, inp, inv, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      ZM_p (s, q - s, p - n - 1, inq, inp, inv, a);
      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+2) */
	  comb (tmp, n + s, n + 2);
	  mpq_mul (a, a, tmp);
	  /* mul (2n) */
	  mpq_set_si (tmp, 2 * n, 1);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);
	  /* div (2n+3) */
	  mpq_set_si (tmp, 2 * n + 3, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_v, coef_v, a);
	}
    }

  append_result (filename, n, p, q, coef_v);

  mpq_clear (a);
  mpq_clear (tmp);
}

/* calc coef Qnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   inq : index of initial condition where Q(inp)00 = 1
 *         in=0 means that all Qn00 = 0
 *   inp : index of initial condition where P(inp)00 = 1
 *         in=0 means that all Pn00 = 0
 *   inv : index of initial condition where V(inv)00 = 1
 *         in=0 means that all Vn00 = 0
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
ZM_q (int n, int p, int q, int inq, int inp, int inv, mpq_t coef_q)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];


  if (n < 0
      || p < 0
      || q < 0)
    {
      /* clear Vnpq */
      mpq_set_ui (coef_q, 0, 1);
      return;
    }

  if (p == 0
      && q == 0)
    {
      /* initial condition */
      if (inq != 0
	  && n == inq)
	mpq_set_ui (coef_q, 1, 1);
      else
	mpq_set_ui (coef_q, 0, 1);

      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.ZM_q.%d.%d.%d", inq, inp, inv);
  if (search_results (filename, n, p, q, coef_q) == 0)
    return; /* found */


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Qnpq */
      mpq_set_ui (a, 0, 1);
      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
      ZM_q (s, q - s - 1, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul s */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-n) */
      ZM_p (s, q - s, p - n, inq, inp, inv, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* div -2 */
	  mpq_set_si (tmp, - 2, 1);
	  mpq_mul (b, b, tmp);
	  /* div (ns) */
	  mpq_set_si (tmp, n * s, 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (n+s,n+2) */
	  comb (tmp, n + s, n + 2);
	  mpq_mul (a, a, tmp);
	  /* div (n+1) */
	  mpq_set_si (tmp, n + 1, 1);
	  mpq_div (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  append_result (filename, n, p, q, coef_q);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}




/* calc comb
 * INPUT
 *   comb : pointer
 *   n, m :
 * OUTPUT
 *   return value : -1 error
 *   comb
 */
int
comb (mpq_t comb, int n, int m)
{
  int i;
  mpq_t n1;

  if(n < 0
     || m < 0
     || n < m)
    {
      mpq_set_ui (comb, 0, 1);      
      return (-1);
    }

  if ((n - m) < m)
    m = n - m;

  /* Create two mrefs */
  mpq_init (n1);

  /* comb = 1 */
  mpq_set_ui (comb, 1, 1);
  for (i=1; i <= m; i++)
    {
      /* comb *= (n + 1 - i) */
      mpq_set_si (n1, n + 1 - i, 1);
      mpq_mul (comb, comb, n1);

      /* comb /= (m + 1 - i) */
      mpq_set_si (n1, m + 1 - i, 1);
      mpq_div (comb, comb, n1);
    }
  mpq_clear (n1);

  return (0);
}

/* output rational number on stdout
 */
void
print_mpq (mpq_t x)
{
  mpz_out_str (stdout, 10, mpq_numref (x));

  if (mpz_cmp_ui (mpq_denref (x), 1))
    {
      printf ("/");
      mpz_out_str (stdout, 10, mpq_denref (x));
    }
}

/* output rational number on stdout
 */
void
fprint_mpq (FILE *out, mpq_t x)
{
  mpz_out_str (out, 10, mpq_numref (x));

  if (mpz_cmp_ui (mpq_denref (x), 1))
    {
      fprintf (out, "/");
      mpz_out_str (out, 10, mpq_denref (x));
    }
}


/* search coef in result file
 * INTPUT
 *   char *file : file name
 *   n, p, q :
 * OUTPUT
 *   mpq_t coef :
 *   (RETURN VALUE) : -1 no entry
 */
int
search_results (char *file, int n, int p, int q, mpq_t coef)
{
  FILE *res;
  mpz_t num, den;
  int nn, pp, qq;
  char string [80];
  int i;

  /* search buffer */
  res = fopen (file, "r");
  if (res != NULL)
    {
      for (;;)
	{
	  i = fscanf (res, "%d %d %d %s\n", &nn, &pp, &qq, string);
	  if (i == EOF)
	    break;

	  if (n == nn
	      && p == pp
	      && q == qq)
	    {
	      mpz_init (num);
	      mpz_init (den);

	      for (i=0;
		   string [i] != '/'
		     && string [i] != '\0';
		   i ++);

	      if (string [i] == '/')
		{
		  string [i] = '\0';
		  mpz_set_str (num, &string [0], 10);
		  mpz_set_str (den, &string [i + 1], 10);
		}
	      else
		{
		  mpz_set_str (num, &string [0], 10);
		  mpz_set_ui (den, 1);
		}

	      mpq_set_num (coef, num);
	      mpq_set_den (coef, den);

	      mpz_clear (num);
	      mpz_clear (den);

	      fclose (res);
	      return (0);
	    }
	}
      fclose (res);
    }
  /* no entry */
  return (-1);
}

/* append result (n, p, q, coef) in result file
 * INTPUT
 *   char *file : file name
 *   n, p, q :
 *   mpq_t coef :
 * OUTPUT
 */
void
append_result (char *file, int n, int p, int q, mpq_t coef)
{
  FILE *res;

  res = fopen (file, "a");
  if (res != NULL)
    {
      fprintf (res, "%d %d %d ", n, p, q);
      fprint_mpq (res, coef);
      fprintf (res, "\n");

      fclose (res);
    }
  else
    fprintf (stderr, "cannot open %s\n", file);
}
