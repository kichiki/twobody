/* coef of 2-body res functions on Jeffrey & Onishi 1984 JFM vol.139 p.261
 * using GMP library
 * Copyright (C) 1999 Kengo ICHIKI (kengo@caltech.edu)
 * $Id$
 */
#include <stdio.h>
#include <stdlib.h> /* malloc() free() */
#include "gmp.h"


/* function prototypes */
int
main (void);

void
make_coef_f (int nmax, mpq_t lambda);

void
print_coef_lambda (int nmax, mpq_t *f, mpq_t labmda, char *label);
void
print_coef (int nmax, mpq_t *f, char *label);

void
XAGP (int nmax, mpq_t *fxa, mpq_t *fxg, mpq_t *fxp);
void
YABG (int nmax, mpq_t *fya, mpq_t *fyb, mpq_t *fyg);
void
XC (int nmax, mpq_t *fxc);
void
YCH (int nmax, mpq_t *fyc, mpq_t *fyh);
void
XMQ (int nmax, mpq_t *fxm, mpq_t *fxq);
void
YM (int nmax, mpq_t *fym);
void
ZM (int nmax, mpq_t *fym);


void
X_problem (int nmax, mpq_t*coef_p, mpq_t*coef_v);
void
Y_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q);
void
XC_problem (int nmax, mpq_t *coef_q);
void
YM_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q);
void
Z_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q);

int
comb (mpq_t comb, int n, int m);
int
pointer_npq (int nmax, int n, int p, int q);
int
pointer_mq (int nmax, int m, int q);

void
print_mpq (mpq_t x);

int
main (void)
{
  int nmax;
  mpq_t lambda;


  nmax = 50;

  mpq_init (lambda);
  mpq_set_ui (lambda, 1, 1);

  make_coef_f (nmax, lambda);

  mpq_clear (lambda);

  return 0;
}

void
make_coef_f (int nmax, mpq_t lambda)
{
  mpq_t *fxa;
  mpq_t *fya;
  mpq_t *fyb;
  mpq_t *fxc;
  mpq_t *fyc;
  mpq_t *fxg;
  mpq_t *fyg;
  mpq_t *fyh;
  mpq_t *fxm;
  mpq_t *fym;
  mpq_t *fzm;
  mpq_t *fxp;
  mpq_t *fxq;

  mpq_t *cur;
  mpq_t a;
  mpq_t b;
  mpq_t lq;

  int i;
  int m, q;


  mpq_init (a);
  mpq_init (b);
  mpq_init (lq);


  fxa = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fya = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fyb = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fxc = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fyc = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fxg = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fyg = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fyh = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fxm = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fym = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fzm = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fxp = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));
  fxq = (mpq_t *) malloc (sizeof (mpq_t)
			   * (nmax + 1) * (nmax + 1));

  /* zero clear */
  for (m=0; m <= nmax; m++)
    {
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = fxa + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fya + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fyb + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fxc + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fyc + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fxg + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fyg + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fyh + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fxm + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fym + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fzm + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fxp + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	  cur = fxq + i;
	  mpq_init (*cur);
	  mpq_set_ui ((*cur), 0, 1);
	}
    }

  /* OK */
  XAGP (nmax, fxa, fxg, fxp);
  YABG (nmax, fya, fyb, fyg);
  XC (nmax, fxc);
  YCH (nmax, fyc, fyh);
  YM (nmax, fym);
  ZM (nmax, fzm);
  /* NG? */
  XMQ (nmax, fxm, fxq);


  /* print coef f_m
  print_coef (nmax, fxa, "XAf");
  print_coef (nmax, fya, "YAf");
  print_coef (nmax, fyb, "YBf");
  print_coef (nmax, fxc, "XCf");
  print_coef (nmax, fyc, "YCf");
  print_coef (nmax, fxg, "XGf");
  print_coef (nmax, fyg, "YGf");
  print_coef (nmax, fyh, "YHf");
  print_coef (nmax, fxm, "XMf");
  print_coef (nmax, fym, "YMf");
  print_coef (nmax, fzm, "ZMf");
  print_coef (nmax, fxp, "XPf");
  print_coef (nmax, fxq, "XQf"); */

  /* calc coef f_m and print */
  print_coef_lambda (nmax, fxa, lambda, "XAf");
  print_coef_lambda (nmax, fya, lambda, "YAf");
  print_coef_lambda (nmax, fyb, lambda, "YBf");
  print_coef_lambda (nmax, fxc, lambda, "XCf");
  print_coef_lambda (nmax, fyc, lambda, "YCf");
  print_coef_lambda (nmax, fxg, lambda, "XGf");
  print_coef_lambda (nmax, fyg, lambda, "YGf");
  print_coef_lambda (nmax, fyh, lambda, "YHf");
  print_coef_lambda (nmax, fxm, lambda, "XMf");
  print_coef_lambda (nmax, fym, lambda, "YMf");
  print_coef_lambda (nmax, fzm, lambda, "ZMf");
  print_coef_lambda (nmax, fxp, lambda, "XPf");
  print_coef_lambda (nmax, fxq, lambda, "XQf");


  /* free */
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (lq);

  for (m=0; m <= nmax; m++)
    {
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = fxa + i;
	  mpq_clear (*cur);
	  cur = fya + i;
	  mpq_clear (*cur);
	  cur = fyb + i;
	  mpq_clear (*cur);
	  cur = fxc + i;
	  mpq_clear (*cur);
	  cur = fyc + i;
	  mpq_clear (*cur);
	  cur = fxg + i;
	  mpq_clear (*cur);
	  cur = fyg + i;
	  mpq_clear (*cur);
	  cur = fyh + i;
	  mpq_clear (*cur);
	  cur = fxm + i;
	  mpq_clear (*cur);
	  cur = fym + i;
	  mpq_clear (*cur);
	  cur = fzm + i;
	  mpq_clear (*cur);
	  cur = fxp + i;
	  mpq_clear (*cur);
	  cur = fxq + i;
	  mpq_clear (*cur);
	}
    }
  free (fxa);
  free (fya);
  free (fyb);
  free (fxc);
  free (fyc);
  free (fxg);
  free (fyg);
  free (fyh);
  free (fxm);
  free (fym);
  free (fzm);
  free (fxp);
  free (fxq);
}

/* calc coef summed up with a lambda and print out
 * INPUT
 *   nmax :
 *   f [nmax][nmax] :
 *   labmda
 *   label
 * OUTPUT (stdout)
 */
void
print_coef_lambda (int nmax, mpq_t *f, mpq_t lambda, char *label)
{
  mpq_t *cur;
  mpq_t a;
  mpq_t b;
  mpq_t lq;

  int i;
  int m, q;


  mpq_init (a);
  mpq_init (b);
  mpq_init (lq);

  for (m=0; m <= nmax; m++)
    {
      mpq_set_ui (a, 0, 1);
      for (q=0, mpq_set_ui (lq, 1, 1);
	   q <= nmax;
	   q++, mpq_mul (lq, lq, lambda))
	{
	  i = pointer_mq (nmax, m, q);

	  cur = f + i;
	  mpq_set (b, (*cur));
	  mpq_mul (b, b, lq);
	  mpq_add (a, a, b);
	}
      fprintf (stdout, "%s [%d] = ",
	       label,
	       m);
      print_mpq (a);
      fprintf (stdout, "\n");
    }

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (lq);
}

/* print coef (without sum up)
 * INPUT
 *   nmax :
 *   f [nmax][nmax] :
 *   labmda
 *   label
 * OUTPUT (stdout)
 */
void
print_coef (int nmax, mpq_t *f, char *label)
{
  mpq_t *cur;

  int i;
  int m, q;

  for (m=0; m <= nmax; m++)
    {
      fprintf (stdout, "%s_%d = ", label, m);
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = f + i;
	  if (mpz_cmp_si (mpq_numref (*cur), 0))
          /*if (mref_notzero (*cur))*/
	    {
	      print_mpq (*cur);
	      fprintf (stdout, " l^%d ", q);
	    }
	}
      fprintf (stdout, "\n");
    }
}

/* calc coef f_m for XA, XG and XP
 * INPUT
 *   nmax
 * OUTPUT
 *   fxa [(nmax + 1) * (nmax + 1)]
 *   fxg [(nmax + 1) * (nmax + 1)]
 *   fxp [(nmax + 1) * (nmax + 1)]
 */
void
XAGP (int nmax, mpq_t *fxa, mpq_t *fxg, mpq_t *fxp)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  X_problem (nmax, coef_p, coef_v);


  mpq_set_ui (two, 2, 1);
  /* XA */
  n = 1;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      cur = fxa + pointer_mq (nmax, k, q);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* XG */
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);
	      /* mul 3 */
	      mpq_set_ui (tmp, 3, 1);
	      mpq_mul (a, a, tmp);
	      /* mul 4 */
	      mpq_set_ui (tmp, 4, 1);
	      mpq_div (a, a, tmp);

	      cur = fxg + pointer_mq (nmax, k, q);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* XP */
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=1; q <= k - 1; q++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  mpq_set (b, (*cur));
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
	  cur = fxp + pointer_mq (nmax, k, q);
	  mpq_set ((*cur), a);
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
}

/* calc coef f_m for YA, YB and YG
 * INPUT
 *   nmax
 * OUTPUT
 *   fya [(nmax + 1) * (nmax + 1)]
 *   fyb [(nmax + 1) * (nmax + 1)]
 *   fyg [(nmax + 1) * (nmax + 1)]
 */
void
YABG (int nmax, mpq_t *fya, mpq_t *fyb, mpq_t *fyg)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *coef_q;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL
      || coef_q == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_q + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  Y_problem (nmax, coef_p, coef_v, coef_q);


  mpq_set_ui (two, 2, 1);
  /* YA */
  n = 1;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      cur = fya + pointer_mq (nmax, k, q);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* YB */
  n = 1;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);
	      /* mul 2 */
	      mpq_set_ui (tmp, 2, 1);
	      mpq_mul (a, a, tmp);

	      cur = fyb + pointer_mq (nmax, k, q);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* YG */
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);
	      /* mul 3 */
	      mpq_set_ui (tmp, 3, 1);
	      mpq_mul (a, a, tmp);
	      /* div 4 */
	      mpq_set_ui (tmp, 4, 1);
	      mpq_div (a, a, tmp);

	      cur = fyg + pointer_mq (nmax, k, q);
	      mpq_set ((*cur), a);
	    }
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	      cur = coef_q + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
  free (coef_q);
}

/* calc coef f_m for XC
 * INPUT
 *   nmax
 * OUTPUT
 *   fxc [(nmax + 1) * (nmax + 1)]
 */
void
XC (int nmax, mpq_t *fxc)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_q;
  mpq_t *cur;
  mpq_t a;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (tmp);

  coef_q = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_q == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* init and zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_q + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }


  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  XC_problem (nmax, coef_q);


  mpq_set_ui (two, 2, 1);
  /* XC */
  n = 1;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fxc + pointer_mq (nmax, k, q);
	      else
		cur = fxc + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (tmp);
  /* free coef_q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_q + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_q);
}

/* calc coef f_m for YC and YH
 * INPUT
 *   nmax
 * OUTPUT
 *   fyc [(nmax + 1) * (nmax + 1)]
 *   fyh [(nmax + 1) * (nmax + 1)]
 */
void
YCH (int nmax, mpq_t *fyc, mpq_t *fyh)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *coef_q;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL
      || coef_q == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_q + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  Y_problem (nmax, coef_p, coef_v, coef_q);


  mpq_set_ui (two, 2, 1);
  /* YC */
  n = 1;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fyc + pointer_mq (nmax, k, q);
	      else
		cur = fyc + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* YH */
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);
	      /* mul -3 */
	      mpq_set_si (tmp, -3, 1);
	      mpq_mul (a, a, tmp);
	      /* div 8 */
	      mpq_set_ui (tmp, 8, 1);
	      mpq_div (a, a, tmp);

	      if (k % 2 == 0)
		cur = fyh + pointer_mq (nmax, k, q);
	      else
		cur = fyh + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	      cur = coef_q + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
  free (coef_q);
}

/* calc coef f_m for XM and XQ
 * INPUT
 *   nmax
 * OUTPUT
 *   fxm [(nmax + 1) * (nmax + 1)]
 *   fxq [(nmax + 1) * (nmax + 1)]
 */
void
XMQ (int nmax, mpq_t *fxm, mpq_t *fxq)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  X_problem (nmax, coef_p, coef_v);


  mpq_set_ui (two, 2, 1);
  /* XM */
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fxm + pointer_mq (nmax, k, q);
	      else
		cur = fxm + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  /* XQ */
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=1; q <= k - 1; q++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  mpq_set (b, (*cur));
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
	  if (k % 2 == 0)
	    cur = fxq + pointer_mq (nmax, k, q);
	  else
	    cur = fxq + pointer_mq (nmax, k, q + 1);
	  mpq_set ((*cur), a);
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
}

/* calc coef f_m for YM
 * INPUT
 *   nmax
 * OUTPUT
 *   fym [(nmax + 1) * (nmax + 1)]
 */
void
YM (int nmax, mpq_t *fym)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *coef_q;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL
      || coef_q == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_q + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  YM_problem (nmax, coef_p, coef_v, coef_q);


  mpq_set_ui (two, 2, 1);
  /* YM */
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fym + pointer_mq (nmax, k, q);
	      else
		cur = fym + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	      cur = coef_q + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
  free (coef_q);
}

/* calc coef f_m for ZM
 * INPUT
 *   nmax
 * OUTPUT
 *   fzm [(nmax + 1) * (nmax + 1)]
 */
void
ZM (int nmax, mpq_t *fzm)
{
  int i;
  int n, p, q;
  int k;

  mpq_t two;
  mpq_t two_k;
  mpq_t *coef_p;
  mpq_t *coef_v;
  mpq_t *coef_q;
  mpq_t *cur;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (two);
  mpq_init (two_k);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  coef_p = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mpq_t *) malloc (sizeof (mpq_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_p == NULL
      || coef_v == NULL
      || coef_q == NULL)
    {
      fprintf (stderr, "cannot alloacted coef_q [%d]\n",
	       (nmax + 1) * (nmax + 1) * (nmax + 1));
      exit (-1);
    }

  /* zero clear */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_v + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);

	      cur = coef_q + i;
	      mpq_init (*cur);
	      mpq_set_ui ((*cur), 0, 1);
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mpq_set_ui ((*cur), 1, 1);

  Z_problem (nmax, coef_p, coef_v, coef_q);


  mpq_set_ui (two, 2, 1);
  n = 2;
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mpq_set (a, (*cur));
	      /* mul (two_k) */
	      mpq_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fzm + pointer_mq (nmax, k, q);
	      else
		cur = fzm + pointer_mq (nmax, k, q + 1);
	      mpq_set ((*cur), a);
	    }
	}
    }

  mpq_clear (two);
  mpq_clear (two_k);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mpq_clear (*cur);
	      cur = coef_v + i;
	      mpq_clear (*cur);
	      cur = coef_q + i;
	      mpq_clear (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
  free (coef_q);
}



void
X_problem (int nmax, mpq_t*coef_p, mpq_t*coef_v)
{
  int i;
  int n, p, q, s;

  mpq_t *cur;
  mpq_t *cur0;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (i=1; i <= nmax; i++)
    {
      for(q = 0; q <= nmax; q++)
	{
	  p = i - q;
	  if (p >= 0)
	    {
	      for (n=1; n <= nmax; n++)
		{
		  /* clear Pnpq and Vnpq */
		  cur = coef_p + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (a, (*cur0));

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

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mpq_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (b, (*cur0));

			  /* mul -n */
			  mpq_set_si (tmp, - n, 1);
			  mpq_mul (b, b, tmp);
			  /* mul (2n-1) */
			  mpq_set_si (tmp, 2 * n - 1, 1);
			  mpq_mul (b, b, tmp);
			  /* div 2 */
			  mpq_set_ui (tmp, 2, 1);
			  mpq_div (b, b, tmp);

			  mpq_add (a, a, b);
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      /* mul comb (n+s,n) */
		      comb (tmp, n + s, n);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  mpq_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

void
Y_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q)
{
  int i;
  int n, p, q, s;

  mpq_t *cur;
  mpq_t *cur0;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t tmp0;


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);

  for (i=1; i <= nmax; i++)
    {
      for(q = 0; q <= nmax; q++)
	{
	  p = i - q;
	  if (p >= 0)
	    {
	      for (n=1; n <= nmax; n++)
		{
		  /* clear Pnpq, Vnpq and Qnpq */
		  cur = coef_p + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n+1) */
		      /*if (p - n + 1 >= 0)*/
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    /*+ pointer_npq (nmax, s, q - s, p - n + 1);*/
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (a, (*cur0));

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

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mpq_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mpq_set (b, (*cur0));

			  /* mul s */
			  mpq_set_si (tmp, s, 1);
			  mpq_mul (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mpq_set (b, (*cur0));

			  /* mul -3 */
			  mpq_set_si (tmp, - 3, 1);
			  mpq_mul (b, b, tmp);

			  /* div (2ns) */
			  mpq_set_si (tmp, 2 * n * s, 1);
			  mpq_div (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
}

void
XC_problem (int nmax, mpq_t *coef_q)
{
  int i;
  int n, p, q, s;

  mpq_t *cur;
  mpq_t *cur0;
  mpq_t a, b;
  mpq_t tmp;


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (i=1; i <= nmax; i++)
    {
      for(q = 0; q <= nmax; q++)
	{
	  p = i - q;
	  if (p >= 0)
	    {
	      for (n=1; n <= nmax; n++)
		{
		  /* clear Qnpq */
		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  for (s=0; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Qnpq : Qs(q-s-1)(p-n) */
		      if (p - n >= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mpq_set (a, (*cur0));

			  /* mul comb (n+s,n) */
			  comb (tmp, n + s, n);
			  mpq_mul (a, a, tmp);
			  /* mul (s) */
			  mpq_set_si (tmp, s, 1);
			  mpq_mul (a, a, tmp);
			  /* div (n+1) */
			  mpq_set_si (tmp, n + 1, 1);
			  mpq_div (a, a, tmp);

			  cur = coef_q + pointer_npq (nmax, n, p, q);
			  mpq_add ((*cur), (*cur), a);
			}
		    }
		}
	    }
	}
    }

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

void
YM_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q)
{
  int i;
  int n, p, q, s;

  mpq_t *cur;
  mpq_t *cur0;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t tmp0;
  mpq_t tmp1;


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  for (i=1; i <= nmax; i++)
    {
      for(q = 0; q <= nmax; q++)
	{
	  p = i - q;
	  if (p >= 0)
	    {
	      for (n=1; n <= nmax; n++)
		{
		  /* clear Pnpq, Vnpq and Qnpq */
		  cur = coef_p + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (a, (*cur0));

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

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mpq_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mpq_set (b, (*cur0));

			  /* mul -(4n^2-1) */
			  mpq_set_si (tmp, - (4 * n * n - 1), 1);
			  mpq_mul (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mpq_set (b, (*cur0));

			  /* mul s */
			  mpq_set_si (tmp, s, 1);
			  mpq_mul (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mpq_set (b, (*cur0));

			  /* div -(ns) */
			  mpq_set_si (tmp, - n * s, 1);
			  mpq_div (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
}

void
Z_problem (int nmax, mpq_t *coef_p, mpq_t *coef_v, mpq_t *coef_q)
{
  int i;
  int n, p, q, s;

  mpq_t *cur;
  mpq_t *cur0;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t tmp0;
  mpq_t tmp1;


  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  for (i=1; i <= nmax; i++)
    {
      for(q = 0; q <= nmax; q++)
	{
	  p = i - q;
	  if (p >= 0)
	    {
	      for (n=1; n <= nmax; n++)
		{
		  /* clear Pnpq, Vnpq and Qnpq */
		  cur = coef_p + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mpq_set_ui ((*cur), 0, 1);

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (a, (*cur0));

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

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mpq_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  mpq_set_ui (b, 1, 1);
			  /* calc (ns+16)(n+s) */
			  mpq_set_si (tmp, n * s + 16, 1);
			  mpq_mul (b, b, tmp);
			  mpq_set_si (tmp, n + s, 1);
			  mpq_mul (b, b, tmp);

			  /* calc -2(ns+1)(ns+4) */
			  mpq_set_si (tmp, -2, 1);
			  mpq_set_si (tmp0, n * s + 1, 1);
			  mpq_mul (tmp, tmp, tmp0);
			  mpq_set_si (tmp0, n * s + 4, 1);
			  mpq_mul (tmp, tmp, tmp0);

			  /* calc (ns+16)(n+s)-2(ns+1)(ns+4) */
			  mpq_add (b, b, tmp);

			  /* mul (2n+1) */
			  mpq_set_si (tmp, 2 * n + 1, 1);
			  mpq_mul (b, b, tmp);
			  /* div 2 */
			  mpq_set_ui (tmp, 2, 1);
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

			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);

			  mpq_mul (b, b, (*cur0));

			  mpq_add (a, a, b);
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mpq_set (b, (*cur0));

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
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mpq_set (b, (*cur0));

			  /* mul -2(4n^2-1) */
			  mpq_set_si (tmp, - 2 * (4 * n * n - 1), 1);
			  mpq_mul (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* mul comb (n+s,n+2) */
		      comb (tmp, n + s, n + 2);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mpq_set_ui (a, 0, 1);
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mpq_set (b, (*cur0));

			  /* mul s */
			  mpq_set_si (tmp, s, 1);
			  mpq_mul (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mpq_set (b, (*cur0));

			  /* div -2 */
			  mpq_set_si (tmp, - 2, 1);
			  mpq_mul (b, b, tmp);
			  /* div (ns) */
			  mpq_set_si (tmp, n * s, 1);
			  mpq_div (b, b, tmp);

			  mpq_add (a, a, b);
			}
		      /* mul comb (n+s,n+2) */
		      comb (tmp, n + s, n + 2);
		      mpq_mul (a, a, tmp);
		      /* div (n+1) */
		      mpq_set_si (tmp, n + 1, 1);
		      mpq_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mpq_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mpq_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
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

  if(n < 0 || m < 0 || n < m)
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

/* return pointer of [nmax][nmax][nmax]
 * INPUT
 *   nmax :
 *   n, p, q :
 * OUTPUT (return value)
 */
int
pointer_npq (int nmax, int n, int p, int q)
{
  return ((n * (nmax + 1) + p) * (nmax + 1) + q);
}

/* return pointer of [nmax][nmax]
 * INPUT
 *   nmax :
 *   m, q :
 * OUTPUT (return value)
 */
int
pointer_mq (int nmax, int m, int q)
{
  return (m * (nmax + 1) + q);
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
