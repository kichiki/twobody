/* coef of 2-body res functions on Jeffrey & Onishi 1984 JFM vol.139 p.261
 * Copyright (C) 1999 Kengo ICHIKI (kengo@caltech.edu)
 * $Id: two-body-JO.c,v 1.4 1999/08/19 20:40:23 ichiki Exp $
 */
#include <stdio.h>
#include <stdlib.h> /* malloc() free() */
#include "saml.h"


/* function prototypes */
int
main (void);

void
make_coef_f (int nmax, mref_t lambda);

void
print_coef_lambda (int nmax, mref_t *f, mref_t labmda, char *label);
void
print_coef (int nmax, mref_t *f, char *label);

void
XAGP (int nmax, mref_t *fxa, mref_t *fxg, mref_t *fxp);
void
YABG (int nmax, mref_t *fya, mref_t *fyb, mref_t *fyg);
void
XC (int nmax, mref_t *fxc);
void
YCH (int nmax, mref_t *fyc, mref_t *fyh);
void
XMQ (int nmax, mref_t *fxm, mref_t *fxq);
void
YM (int nmax, mref_t *fym);
void
ZM (int nmax, mref_t *fym);


void
X_problem (int nmax, mref_t*coef_p, mref_t*coef_v);
void
Y_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q);
void
XC_problem (int nmax, mref_t *coef_q);
void
YM_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q);
void
Z_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q);

int
comb (mref_t comb, int n, int m);
int
pointer_npq (int nmax, int n, int p, int q);
int
pointer_mq (int nmax, int m, int q);


int
main (void)
{
  int nmax;
  mref_t lambda;


  nmax = 11;
  lambda = mref_new ();
  mref_build (lambda, ST_RATIONAL, "1");

  make_coef_f (nmax, lambda);

  mref_free (lambda);

  return 0;
}

void
make_coef_f (int nmax, mref_t lambda)
{
  mref_t *fxa;
  mref_t *fya;
  mref_t *fyb;
  mref_t *fxc;
  mref_t *fyc;
  mref_t *fxg;
  mref_t *fyg;
  mref_t *fyh;
  mref_t *fxm;
  mref_t *fym;
  mref_t *fzm;
  mref_t *fxp;
  mref_t *fxq;

  mref_t *cur;
  mref_t a;
  mref_t b;
  mref_t lq;

  int i;
  int m, q;


  a = mref_new ();
  b = mref_new ();
  lq = mref_new ();


  fxa = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fya = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fyb = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fxc = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fyc = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fxg = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fyg = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fyh = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fxm = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fym = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fzm = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fxp = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));
  fxq = (mref_t *) malloc (sizeof (mref_t)
			   * (nmax + 1) * (nmax + 1));

  /* zero clear */
  for (m=0; m <= nmax; m++)
    {
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = fxa + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fya + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fyb + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fxc + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fyc + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fxg + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fyg + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fyh + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fxm + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fym + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fzm + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fxp + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	  cur = fxq + i;
	  (*cur) = mref_new ();
	  mref_build ((*cur), ST_RATIONAL, "0");
	}
    }

  /* OK */
  XAGP (nmax, fxa, fxg, fxp);
  XC (nmax, fxc);
  YCH (nmax, fyc, fyh);
  YM (nmax, fym);
  ZM (nmax, fzm);
  YABG (nmax, fya, fyb, fyg);
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
  mref_free (a);
  mref_free (b);
  mref_free (lq);

  for (m=0; m <= nmax; m++)
    {
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = fxa + i;
	  mref_free (*cur);
	  cur = fya + i;
	  mref_free (*cur);
	  cur = fyb + i;
	  mref_free (*cur);
	  cur = fxc + i;
	  mref_free (*cur);
	  cur = fyc + i;
	  mref_free (*cur);
	  cur = fxg + i;
	  mref_free (*cur);
	  cur = fyg + i;
	  mref_free (*cur);
	  cur = fyh + i;
	  mref_free (*cur);
	  cur = fxm + i;
	  mref_free (*cur);
	  cur = fym + i;
	  mref_free (*cur);
	  cur = fzm + i;
	  mref_free (*cur);
	  cur = fxp + i;
	  mref_free (*cur);
	  cur = fxq + i;
	  mref_free (*cur);
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
print_coef_lambda (int nmax, mref_t *f, mref_t lambda, char *label)
{
  mref_t *cur;
  mref_t a;
  mref_t b;
  mref_t lq;

  int i;
  int m, q;


  a = mref_new ();
  b = mref_new ();
  lq = mref_new ();

  for (m=0; m <= nmax; m++)
    {
      mref_build (a, ST_RATIONAL, "0");
      for (q=0, mref_build (lq, ST_RATIONAL, "1");
	   q <= nmax;
	   q++, mref_mul (lq, lq, lambda))
	{
	  i = pointer_mq (nmax, m, q);

	  cur = f + i;
	  mref_copy (b, (*cur));
	  mref_mul (b, b, lq);
	  mref_add (a, a, b);
	}
      fprintf (stdout, "%s [%d] = %s\n",
	       label,
	       m,
	       mref_string (a));
    }

  mref_free (a);
  mref_free (b);
  mref_free (lq);
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
print_coef (int nmax, mref_t *f, char *label)
{
  mref_t *cur;

  int i;
  int m, q;

  for (m=0; m <= nmax; m++)
    {
      fprintf (stdout, "%s_%d = ", label, m);
      for (q=0; q <= nmax; q++)
	{
	  i = pointer_mq (nmax, m, q);

	  cur = f + i;
	  if (mref_notzero (*cur))
	    {
	      fprintf (stdout, "%s", mref_string (*cur));
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
XAGP (int nmax, mref_t *fxa, mref_t *fxg, mref_t *fxp)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  X_problem (nmax, coef_p, coef_v);


  mref_build (two, ST_RATIONAL, "2");
  /* XA */
  n = 1;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      cur = fxa + pointer_mq (nmax, k, q);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* XG */
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);
	      /* mul 3 */
	      mref_build (tmp, ST_RATIONAL, "3");
	      mref_mul (a, a, tmp);
	      /* mul 4 */
	      mref_build (tmp, ST_RATIONAL, "4");
	      mref_div (a, a, tmp);

	      cur = fxg + pointer_mq (nmax, k, q);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* XP */
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=1; q <= k - 1; q++)
	{
	  mref_build (a, ST_RATIONAL, "0");
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  mref_copy (b, (*cur));
		  /* mul (two_k) */
		  mref_mul (b, b, two_k);
		  /* mul 3 */
		  mref_build (tmp, ST_RATIONAL, "3");
		  mref_mul (b, b, tmp);
		  /* mul 2 */
		  mref_build (tmp, ST_RATIONAL, "2");
		  mref_div (b, b, tmp);

		  mref_add (a, a, b);
		}
	    }
	  cur = fxp + pointer_mq (nmax, k, q);
	  mref_copy ((*cur), a);
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
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
YABG (int nmax, mref_t *fya, mref_t *fyb, mref_t *fyg)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *coef_q;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_q + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  Y_problem (nmax, coef_p, coef_v, coef_q);

  /*for (n=1; n<=nmax; n++)
    for (p=0; p<=nmax; p++)
      for (q=0; q<=nmax; q++)
	{
	  fprintf (stdout, "P_%d %d %d = ", n, p, q);
	  cur = coef_p + pointer_npq (nmax, n, p, q);
	  fprintf (stdout, "%s", mref_string (*cur));
	  fprintf (stdout, " \n");
	}
  exit (0);*/

  mref_build (two, ST_RATIONAL, "2");
  /* YA */
  n = 1;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      cur = fya + pointer_mq (nmax, k, q);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* YB */
  n = 1;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);
	      /* mul 2 */
	      mref_build (tmp, ST_RATIONAL, "2");
	      mref_mul (a, a, tmp);

	      cur = fyb + pointer_mq (nmax, k, q);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* YG */
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);
	      /* mul 3 */
	      mref_build (tmp, ST_RATIONAL, "3");
	      mref_mul (a, a, tmp);
	      /* div 4 */
	      mref_build (tmp, ST_RATIONAL, "4");
	      mref_div (a, a, tmp);

	      cur = fyg + pointer_mq (nmax, k, q);
	      mref_copy ((*cur), a);
	    }
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
	      cur = coef_q + i;
	      mref_free (*cur);
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
XC (int nmax, mref_t *fxc)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_q;
  mref_t *cur;
  mref_t a;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  tmp = mref_new ();

  coef_q = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }


  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  XC_problem (nmax, coef_q);


  mref_build (two, ST_RATIONAL, "2");
  /* XC */
  n = 1;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fxc + pointer_mq (nmax, k, q);
	      else
		cur = fxc + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (tmp);
  /* free coef_q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_q + i;
	      mref_free (*cur);
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
YCH (int nmax, mref_t *fyc, mref_t *fyh)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *coef_q;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_q + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  Y_problem (nmax, coef_p, coef_v, coef_q);


  mref_build (two, ST_RATIONAL, "2");
  /* YC */
  n = 1;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fyc + pointer_mq (nmax, k, q);
	      else
		cur = fyc + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* YH */
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);
	      /* mul -3 */
	      mref_build (tmp, ST_RATIONAL, "-3");
	      mref_mul (a, a, tmp);
	      /* div 8 */
	      mref_build (tmp, ST_RATIONAL, "8");
	      mref_div (a, a, tmp);

	      if (k % 2 == 0)
		cur = fyh + pointer_mq (nmax, k, q);
	      else
		cur = fyh + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
	      cur = coef_q + i;
	      mref_free (*cur);
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
XMQ (int nmax, mref_t *fxm, mref_t *fxq)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  X_problem (nmax, coef_p, coef_v);


  /*for (n=1; n<=nmax; n++)
    for (p=0; p<=nmax; p++)
      for (q=0; q<=nmax; q++)
	{
	  fprintf (stdout, "P_%d %d %d = ", n, p, q);
	  cur = coef_p + pointer_npq (nmax, n, p, q);
	  fprintf (stdout, "%s", mref_string (*cur));
	  fprintf (stdout, " \n");
	}
  exit (0);*/

  mref_build (two, ST_RATIONAL, "2");
  /* XM */
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fxm + pointer_mq (nmax, k, q);
	      else
		cur = fxm + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  /* XQ */
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=1; q <= k - 1; q++)
	{
	  mref_build (a, ST_RATIONAL, "0");
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  mref_copy (b, (*cur));
		  /* mul (two_k) */
		  mref_mul (b, b, two_k);
		  /* mul 5 */
		  mref_build (tmp, ST_RATIONAL, "5");
		  mref_mul (b, b, tmp);
		  /* mul 2 */
		  mref_build (tmp, ST_RATIONAL, "2");
		  mref_div (b, b, tmp);

		  mref_add (a, a, b);
		}
	    }
	  if (k % 2 == 0)
	    cur = fxq + pointer_mq (nmax, k, q);
	  else
	    cur = fxq + pointer_mq (nmax, k, q + 1);
	  mref_copy ((*cur), a);
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
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
YM (int nmax, mref_t *fym)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *coef_q;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_q + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  YM_problem (nmax, coef_p, coef_v, coef_q);


  mref_build (two, ST_RATIONAL, "2");
  /* YM */
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fym + pointer_mq (nmax, k, q);
	      else
		cur = fym + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
	      cur = coef_q + i;
	      mref_free (*cur);
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
ZM (int nmax, mref_t *fzm)
{
  int i;
  int n, p, q;
  int k;

  mref_t two;
  mref_t two_k;
  mref_t *coef_p;
  mref_t *coef_v;
  mref_t *coef_q;
  mref_t *cur;
  mref_t a, b;
  mref_t tmp;


  two = mref_new ();
  two_k = mref_new ();
  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

  coef_p = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (mref_t *) malloc (sizeof (mref_t)
			      * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (mref_t *) malloc (sizeof (mref_t)
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
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_v + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");

	      cur = coef_q + i;
	      (*cur) = mref_new ();
	      mref_build ((*cur), ST_RATIONAL, "0");
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  mref_build ((*cur), ST_RATIONAL, "1");

  Z_problem (nmax, coef_p, coef_v, coef_q);


  /*for (n=1; n<=nmax; n++)
    for (p=0; p<=nmax; p++)
      for (q=0; q<=nmax; q++)
	{
	  fprintf (stdout, "P_%d %d %d = ", n, p, q);
	  cur = coef_p + pointer_npq (nmax, n, p, q);
	  fprintf (stdout, "%s", mref_string (*cur));
	  fprintf (stdout, " \n");
	}
  exit (0);*/

  mref_build (two, ST_RATIONAL, "2");
  n = 2;
  for (k=0, mref_build (two_k, ST_RATIONAL, "1");
       k <= nmax;
       k++, mref_mul (two_k, two_k, two))
    {
      for(q=0; q </*=*/ nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      mref_copy (a, (*cur));
	      /* mul (two_k) */
	      mref_mul (a, a, two_k);

	      if (k % 2 == 0)
		cur = fzm + pointer_mq (nmax, k, q);
	      else
		cur = fzm + pointer_mq (nmax, k, q + 1);
	      mref_copy ((*cur), a);
	    }
	}
    }

  mref_free (two);
  mref_free (two_k);
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  /* free coef_p,v,q */
  for (n=0; n <= nmax; n++)
    {
      for (p=0; p <= nmax; p++)
	{
	  for (q=0; q <= nmax; q++)
	    {
	      i = pointer_npq (nmax, n, p, q);

	      cur = coef_p + i;
	      mref_free (*cur);
	      cur = coef_v + i;
	      mref_free (*cur);
	      cur = coef_q + i;
	      mref_free (*cur);
	    }
	}
    }
  free (coef_p);
  free (coef_v);
  free (coef_q);
}



void
X_problem (int nmax, mref_t*coef_p, mref_t*coef_v)
{
  char string [20];
  int i;
  int n, p, q, s;

  mref_t *cur;
  mref_t *cur0;
  mref_t a, b;
  mref_t tmp;


  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

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
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (a, (*cur0));

			  /* mul comb (n+s,n) */
			  comb (tmp, n + s, n);
			  mref_mul (a, a, tmp);
			  /* mul (-2n) */
			  sprintf (string, "%d", - 2 * n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (a, a, tmp);
			  /* div (n+1) */
			  sprintf (string, "%d", n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);
			  /* div (2n+3) */
			  sprintf (string, "%d", 2 * n + 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mref_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2n+1) */
			  sprintf (string, "%d", 2 * n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2ns-n-s+2) */
			  sprintf (string, "%d", 2 * n * s - n - s + 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2(2s-1) */
			  sprintf (string, "%d", 2 * (2 * s - 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (n+s) */
			  sprintf (string, "%d", n + s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (b, (*cur0));

			  /* mul -n */
			  sprintf (string, "%d", - n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2n-1) */
			  sprintf (string, "%d", 2 * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul -n */
			  sprintf (string, "%d", - n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (4n^2-1) */
			  sprintf (string, "%d", 4 * n * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2(2s+1) */
			  sprintf (string, "%d", 2 * (2 * s + 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n) */
		      comb (tmp, n + s, n);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  mref_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
}

void
Y_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q)
{
  char string [20];
  int i;
  int n, p, q, s;

  mref_t *cur;
  mref_t *cur0;
  mref_t a, b;
  mref_t tmp;
  mref_t tmp0;


  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();
  tmp0 = mref_new ();

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
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

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
			  mref_copy (a, (*cur0));

			  /* mul comb (n+s,n+1) */
			  comb (tmp, n + s, n + 1);
			  mref_mul (a, a, tmp);
			  /* mul (2n) */
			  sprintf (string, "%d", 2 * n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (a, a, tmp);
			  /* div (n+1) */
			  sprintf (string, "%d", n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);
			  /* div (2n+3) */
			  sprintf (string, "%d", 2 * n + 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mref_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul (2n+1) */
			  sprintf (string, "%d", 2 * n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  /* mul 3(n+s)-(ns+1)(2ns-s-n+2) */
			  /* calc -(ns+1)(2ns-s-n+2) */
			  sprintf (string, "%d", - (n * s + 1));
			  mref_build (tmp0, ST_RATIONAL, string);
			  sprintf (string, "%d", 2 * n * s - s - n + 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (tmp, tmp, tmp0);
			  /* add 3(n+s) */
			  sprintf (string, "%d", 3 * (n + s));
			  mref_build (tmp0, ST_RATIONAL, string);
			  mref_add (tmp, tmp, tmp0);
			  mref_mul (b, b, tmp);

			  /* div (s) */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (n+s) */
			  sprintf (string, "%d", n + s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (2s-1) */
			  sprintf (string, "%d", 2 * s - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2n-1) */
			  sprintf (string, "%d", 2 * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (4n^2-1) */
			  sprintf (string, "%d", 4 * n * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2(2s+1) */
			  sprintf (string, "%d", 2 * (2 * s + 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul -2 */
			  sprintf (string, "%d", - 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (4n^2-1) */
			  sprintf (string, "%d", 4 * n * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 3 */
			  sprintf (string, "%d", 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mref_copy (b, (*cur0));

			  /* mul s */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mref_copy (b, (*cur0));

			  /* mul -3 */
			  sprintf (string, "%d", - 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  /* div (2ns) */
			  sprintf (string, "%d", 2 * n * s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  mref_free (tmp0);
}

void
XC_problem (int nmax, mref_t *coef_q)
{
  char string [20];
  int i;
  int n, p, q, s;

  mref_t *cur;
  mref_t *cur0;
  mref_t a, b;
  mref_t tmp;


  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();

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
		  mref_build ((*cur), ST_RATIONAL, "0");

		  for (s=0; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Qnpq : Qs(q-s-1)(p-n) */
		      if (p - n >= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mref_copy (a, (*cur0));

			  /* mul comb (n+s,n) */
			  comb (tmp, n + s, n);
			  mref_mul (a, a, tmp);
			  /* mul (s) */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (a, a, tmp);
			  /* div (n+1) */
			  sprintf (string, "%d", n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);

			  cur = coef_q + pointer_npq (nmax, n, p, q);
			  mref_add ((*cur), (*cur), a);
			}
		    }
		}
	    }
	}
    }

  mref_free (a);
  mref_free (b);
  mref_free (tmp);
}

void
YM_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q)
{
  char string [20];
  int i;
  int n, p, q, s;

  mref_t *cur;
  mref_t *cur0;
  mref_t a, b;
  mref_t tmp;
  mref_t tmp0;
  mref_t tmp1;


  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();
  tmp0 = mref_new ();
  tmp1 = mref_new ();

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
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (a, (*cur0));

			  /* mul comb (n+s,n+1) */
			  comb (tmp, n + s, n + 1);
			  mref_mul (a, a, tmp);
			  /* mul (2n) */
			  sprintf (string, "%d", 2 * n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (a, a, tmp);
			  /* div (n+1) */
			  sprintf (string, "%d", n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);
			  /* div (2n+3) */
			  sprintf (string, "%d", 2 * n + 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mref_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul (2n+1) */
			  sprintf (string, "%d", 2 * n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  /* mul (ns+4)(n+s)-2(ns+1)^2 */
			  /* calc -2(ns+1)^2 */
			  sprintf (string, "%d", (n * s + 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_build (tmp0, ST_RATIONAL, string);
			  mref_mul (tmp, tmp, tmp0);
			  sprintf (string, "%d", - 2);
			  mref_build (tmp0, ST_RATIONAL, string);
			  mref_mul (tmp, tmp, tmp0);
			  /* add (ns+4)(n+s) */
			  sprintf (string, "%d", n + s);
			  mref_build (tmp0, ST_RATIONAL, string);
			  sprintf (string, "%d", n * s + 4);
			  mref_build (tmp1, ST_RATIONAL, string);
			  mref_mul (tmp0, tmp0, tmp1);
			  mref_add (tmp, tmp, tmp0);
			  /* mul it with b */
			  mref_mul (b, b, tmp);

			  /* div (s) */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (2s-1) */
			  sprintf (string, "%d", 2 * s - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (n+s) */
			  sprintf (string, "%d", n + s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2n-1) */
			  sprintf (string, "%d", 2 * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (4n^2-1) */
			  sprintf (string, "%d", 4 * n * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2(2s+1) */
			  sprintf (string, "%d", 2 * (2 * s + 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul -(4n^2-1) */
			  sprintf (string, "%d", - (4 * n * n - 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mref_copy (b, (*cur0));

			  /* mul s */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mref_copy (b, (*cur0));

			  /* div -(ns) */
			  sprintf (string, "%d", - n * s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+1) */
		      comb (tmp, n + s, n + 1);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  mref_free (tmp0);
  mref_free (tmp1);
}

void
Z_problem (int nmax, mref_t *coef_p, mref_t *coef_v, mref_t *coef_q)
{
  char string [20];
  int i;
  int n, p, q, s;

  mref_t *cur;
  mref_t *cur0;
  mref_t a, b;
  mref_t tmp;
  mref_t tmp0;
  mref_t tmp1;


  a = mref_new ();
  b = mref_new ();
  tmp = mref_new ();
  tmp0 = mref_new ();
  tmp1 = mref_new ();

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
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  mref_build ((*cur), ST_RATIONAL, "0");

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (a, (*cur0));

			  /* mul comb (n+s,n+2) */
			  comb (tmp, n + s, n + 2);
			  mref_mul (a, a, tmp);
			  /* mul (2n) */
			  sprintf (string, "%d", 2 * n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (a, a, tmp);
			  /* div (n+1) */
			  sprintf (string, "%d", n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);
			  /* div (2n+3) */
			  sprintf (string, "%d", 2 * n + 3);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (a, a, tmp);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  mref_add ((*cur), (*cur), a);
			}

		      /* Pnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  mref_build (b, ST_RATIONAL, "1");
			  /* calc (ns+16)(n+s) */
			  sprintf (string, "%d", n * s + 16);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  sprintf (string, "%d", n + s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  /* calc -2(ns+1)(ns+4) */
			  mref_build (tmp, ST_RATIONAL, "-2");
			  sprintf (string, "%d", n * s + 1);
			  mref_build (tmp0, ST_RATIONAL, string);
			  mref_mul (tmp, tmp, tmp0);
			  sprintf (string, "%d", n * s + 4);
			  mref_build (tmp0, ST_RATIONAL, string);
			  mref_mul (tmp, tmp, tmp0);

			  /* calc (ns+16)(n+s)-2(ns+1)(ns+4) */
			  mref_add (b, b, tmp);

			  /* mul (2n+1) */
			  sprintf (string, "%d", 2 * n + 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  mref_build (tmp, ST_RATIONAL, "2");
			  mref_div (b, b, tmp);

			  /* div (s) */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (2s-1) */
			  sprintf (string, "%d", 2 * s - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);
			  /* div (n+s) */
			  sprintf (string, "%d", n + s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);

			  mref_mul (b, b, (*cur0));

			  mref_add (a, a, b);
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (2n-1) */
			  sprintf (string, "%d", 2 * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2 */
			  sprintf (string, "%d", 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul n */
			  sprintf (string, "%d", n);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* mul (4n^2-1) */
			  sprintf (string, "%d", 4 * n * n - 1);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div 2(2s+1) */
			  sprintf (string, "%d", 2 * (2 * s + 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  mref_copy (b, (*cur0));

			  /* mul -2(4n^2-1) */
			  sprintf (string, "%d", - 2 * (4 * n * n - 1));
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+2) */
		      comb (tmp, n + s, n + 2);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);

		      /* Qnpq */
		      mref_build (a, ST_RATIONAL, "0");
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  mref_copy (b, (*cur0));

			  /* mul s */
			  sprintf (string, "%d", s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  mref_copy (b, (*cur0));

			  /* div -2 */
			  sprintf (string, "%d", - 2);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_mul (b, b, tmp);
			  /* div (ns) */
			  sprintf (string, "%d", n * s);
			  mref_build (tmp, ST_RATIONAL, string);
			  mref_div (b, b, tmp);

			  mref_add (a, a, b);
			}
		      /* mul comb (n+s,n+2) */
		      comb (tmp, n + s, n + 2);
		      mref_mul (a, a, tmp);
		      /* div (n+1) */
		      sprintf (string, "%d", n + 1);
		      mref_build (tmp, ST_RATIONAL, string);
		      mref_div (a, a, tmp);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      mref_add ((*cur), (*cur), a);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p + pointer_npq (nmax, n, p, q);
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  mref_add ((*cur), (*cur), (*cur0));
		}
	    }
	}
    }
  mref_free (a);
  mref_free (b);
  mref_free (tmp);
  mref_free (tmp0);
  mref_free (tmp1);
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
comb (mref_t comb, int n, int m)
{
  int i;
  mref_t n1;
  char string [20];

  if(n < 0 || m < 0 || n < m)
    {
      mref_build (comb, ST_RATIONAL, "0");      
      return (-1);
    }

  if ((n - m) < m)
    m = n - m;

  /* Create two mrefs */
  n1 = mref_new ();

  /* comb = 1 */
  mref_build (comb, ST_RATIONAL, "1");
  for (i=1; i <= m; i++)
    {
      /* comb *= (n + 1 - i) */
      sprintf (string, "%d", n + 1 - i);
      mref_build (n1, ST_RATIONAL, string);
      mref_mul (comb, comb, n1);

      /* comb /= (m + 1 - i) */
      sprintf (string, "%d", m + 1 - i);
      mref_build (n1, ST_RATIONAL, string);
      mref_div (comb, comb, n1);
    }
  mref_free (n1);

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
