/* coef of 2-body res functions on Jeffrey & Onishi 1984 JFM vol.139 p.261
 * Copyright (C) 1999 Kengo ICHIKI (kengo@caltech.edu)
 * $Id: two-body-JO.c,v 1.3 1999/08/19 01:22:44 ichiki Exp $
 */
#include <math.h>
#include <stdio.h> /* printf() fprintf() */
#include <stdlib.h> /* exit() */

#include "ratio.h"

/* global variable */

/* function prototype */
int
main (int argc, char** argv);

void
XAGP (int nmax);
void
YABG (int nmax);
void
YCH (int nmax);
void
XC (int nmax);
void
XMQ (int nmax);
void
YM (int nmax);
void
ZM (int nmax);

void
X_problem (int nmax, ratio *coef_p, ratio *coef_v);
void
Y_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q);
void
XC_problem (int nmax, ratio *coef_q);
void
YM_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q);
void
Z_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q);


double
combi (int n, int m);
int
comb (int n, int m);

int
pointer_npq (int nmax, int n, int p, int q);


/* main program */
int
main (int argc, char** argv)
{
  char flg;
  int i;
  int nmax;


  /* option analysis */
  flg = 0;
  for (i=1; i<argc; i++)
    {
      if (strcmp (argv [i], "-n") == 0)
	{
	  if (i+1 < argc)
	    {
	      nmax = atoi (argv [++i]);
	      flg ++;
	    }
	}
      else
	{
	  fprintf (stderr, "USAGE\n");
	  fprintf (stderr, "XC [options]\n");
	  fprintf (stderr, "OPTIONS\n");
	  fprintf (stderr, "  -n number of coefficients\n");
	  exit (1);
	}
    }
  if (flg != 1)
    {
      fprintf (stderr, "You Must Believe In Spring.\n");
      exit (1);
    }

  XAGP (nmax);
  XMQ (nmax);
  /*XC (nmax);
  YCH (nmax);
  YM (nmax);
  ZM (nmax);*/
  /*YABG (nmax);*/

  return (0);
}


void
XAGP (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *cur;
  ratio a, b, c;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  X_problem (nmax, coef_p, coef_v);


  fprintf (stdout, "for XA\n");
  n = 1;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for XG\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k * 3;
	      a.den = (*cur).den * 4;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for XP\n");
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=1; q <= k - 1; q++)
	{
	  a.sgn = 1;
	  a.num = 0;
	  a.den = 1;
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num * two_k * 3;
		  b.den = (*cur).den * 2;
		  add_ratio (a, b, &c);

		  a.sgn = c.sgn;
		  a.num = c.num;
		  a.den = c.den;
		}
	    }
	  if (a.num != 0)
	    {
	      if (a.sgn > 0)
		fprintf (stdout, "+");
	      print_ratio (a);
	      fprintf (stdout, "l^%d ", q);
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
}

void
YABG (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *coef_q;
  ratio *cur;
  ratio a;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_q + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  cur = coef_v + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  Y_problem (nmax, coef_p, coef_v, coef_q);


  fprintf (stdout, "for YA\n");
  n = 1;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for YB\n");
  n = 1;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k * 2;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for YG\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k * 3;
	      a.den = (*cur).den * 4;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
  free (coef_q);
}

void
YCH (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *coef_q;
  ratio *cur;
  ratio a;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_q + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  Y_problem (nmax, coef_p, coef_v, coef_q);


  fprintf (stdout, "for YC\n");
  n = 1;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for YH\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = - (*cur).sgn;
	      a.num = (*cur).num * two_k * 3;
	      a.den = (*cur).den * 8;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
  free (coef_q);
}


void
XC (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_q;
  ratio *cur;
  ratio a;


  coef_q = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  if (coef_q == NULL)
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

	      cur = coef_q + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_q + pointer_npq (nmax, 1, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  XC_problem (nmax, coef_q);


  fprintf (stdout, "for XC\n");
  n = 1;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_q + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  fprintf (stdout, "l^%d ", q);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_q);
}

void
XMQ (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *cur;
  ratio a, b, c;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  X_problem (nmax, coef_p, coef_v);


  fprintf (stdout, "for XM\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  if (k % 2 == 0)
		    fprintf (stdout, "l^%d ", q);
		  else
		    fprintf (stdout, "l^%d ", q + 1);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  fprintf (stdout, "for XQ\n");
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=1; q <= k - 1; q++)
	{
	  a.sgn = 1;
	  a.num = 0;
	  a.den = 1;
	  for (n=1; n <= (q + 1) / 2; n++)
	    {
	      if (q - n >= 0
		  && k - q - 1 >= 0)
		{
		  cur = coef_p + pointer_npq (nmax, n, q - n, k - q - 1);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num * two_k * 5;
		  b.den = (*cur).den * 2;
		  add_ratio (a, b, &c);

		  a.sgn = c.sgn;
		  a.num = c.num;
		  a.den = c.den;
		}
	    }
	  if (a.num != 0)
	    {
	      if (a.sgn > 0)
		fprintf (stdout, "+");
	      print_ratio (a);
	      fprintf (stdout, "l^%d ", q);
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
}

void
YM (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *coef_q;
  ratio *cur;
  ratio a;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_q + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  YM_problem (nmax, coef_p, coef_v, coef_q);


  fprintf (stdout, "for YM\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  if (k % 2 == 0)
		    fprintf (stdout, "l^%d ", q);
		  else
		    fprintf (stdout, "l^%d ", q + 1);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
  free (coef_q);
}

void
ZM (int nmax)
{
  int i;
  int n, p, q;
  int k;
  int two_k;

  ratio *coef_p;
  ratio *coef_v;
  ratio *coef_q;
  ratio *cur;
  ratio a;


  coef_p = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_v = (ratio *) malloc (sizeof (ratio)
			     * (nmax + 1) * (nmax + 1) * (nmax + 1));
  coef_q = (ratio *) malloc (sizeof (ratio)
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
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_v + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;

	      cur = coef_q + i;
	      (*cur).sgn = 1;
              (*cur).num = 0;
	      (*cur).den = 1;
	    }
	}
    }
  /* initial condition */
  cur = coef_p + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  cur = coef_v + pointer_npq (nmax, 2, 0, 0);
  (*cur).sgn = 1;
  (*cur).num = 1;
  (*cur).den = 1;

  Z_problem (nmax, coef_p, coef_v, coef_q);


  fprintf (stdout, "for ZM\n");
  n = 2;
  for (k=0, two_k = 1;
       k <= nmax;
       k++, two_k *= 2)
    {
      fprintf (stdout, "f_%d = ", k);
      for(q=0; q <= nmax; q++)
	{
	  p = k - q;
	  if (p >= 0)
	    {
	      cur = coef_p + pointer_npq (nmax, n, p, q);
	      a.sgn = (*cur).sgn;
	      a.num = (*cur).num * two_k;
	      a.den = (*cur).den;
	      reduce_ratio (&a);
	      if (a.num != 0)
		{
		  if (a.sgn > 0)
		    fprintf (stdout, "+");
		  print_ratio (a);
		  if (k % 2 == 0)
		    fprintf (stdout, "l^%d ", q);
		  else
		    fprintf (stdout, "l^%d ", q + 1);
		}
	    }
	}
      fprintf (stdout, "\n");
    }

  free (coef_p);
  free (coef_v);
  free (coef_q);
}



void
X_problem (int nmax, ratio *coef_p, ratio *coef_v)
{
  int i;
  int n, p, q, s;

  ratio *cur;
  ratio *cur0;
  ratio a, b, c;


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
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;
		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  a.sgn = - (*cur0).sgn;
			  a.num = (*cur0).num
			    * comb (n + s, n)
			    * 2 * n;
			  a.den = (*cur0).den
			    * (n + 1) * (2 * n + 3);
			  reduce_ratio (&a);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  b.sgn = (*cur).sgn;
			  b.num = (*cur).num;
			  b.den = (*cur).den;

			  add_ratio (a, b, cur);
			}

		      /* Pnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (2 * n + 1) * (2 * n * s - n - s + 2);
			  b.den = (*cur0).den
			    * 2 * (2 * s - 1) * (n + s);
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (2 * n - 1);
			  b.den = (*cur0).den
			    * 2;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (4 * n * n - 1);
			  b.den = (*cur0).den
			    * 2 * (2 * s + 1);
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      a.num *= comb (n + s, n);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p
		    + pointer_npq (nmax, n, p, q);
		  a.sgn = (*cur0).sgn;
		  a.num = (*cur0).num;
		  a.den = (*cur0).den;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num;
		  b.den = (*cur).den;

		  add_ratio (a, b, cur);
		}
	    }
	}
    }
}

void
Y_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q)
{
  int i;
  int n, p, q, s;

  ratio *cur;
  ratio *cur0;
  ratio a, b, c;

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
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

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
			  a.sgn = (*cur0).sgn;
			  a.num = (*cur0).num
			    * comb (n + s, n + 1)
			    * 2 * n;
			  a.den = (*cur0).den
			    * (n + 1) * (2 * n + 3);
			  reduce_ratio (&a);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  b.sgn = (*cur).sgn;
			  b.num = (*cur).num;
			  b.den = (*cur).den;

			  add_ratio (a, b, cur);
			}

		      /* Pnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  b.sgn = 1;
			  b.num = 2 * n + 1;
			  b.den = 2;
			  reduce_ratio (&b);

			  b.num *=
			    3 * (n + s)
			    - (n * s  + 1) * (2 * n * s - s - n + 2);
			  b.den *=
			    s * (n + s) * (2 * s - 1);
			  reduce_ratio (&b);

			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  b.sgn *= (*cur0).sgn;
			  b.num *= (*cur0).num;
			  b.den *= (*cur0).den;
			  reduce_ratio (&b);

			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (2 * n - 1);
			  b.den = (*cur0).den
			    * 2;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (4 * n * n - 1);
			  b.den = (*cur0).den
			    * 2 * (2 * s + 1);
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * 2 * (4 * n * n - 1);
			  b.den = (*cur0).den
			    * 3;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 1);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);

		      /* Qnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * s;
			  b.den = (*cur0).den;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * 3;
			  b.den = (*cur0).den
			    * 2 * n * s;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 1);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p
		    + pointer_npq (nmax, n, p, q);
		  a.sgn = (*cur0).sgn;
		  a.num = (*cur0).num;
		  a.den = (*cur0).den;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num;
		  b.den = (*cur).den;

		  add_ratio (a, b, cur);
		}
	    }
	}
    }
}

void
XC_problem (int nmax, ratio *coef_q)
{
  int i;
  int n, p, q, s;

  ratio *cur;
  ratio *cur0;
  ratio a, b;

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
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  for (s=0; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Qnpq : Qs(q-s-1)(p-n) */
		      if (p - n >= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  a.sgn = (*cur0).sgn;
			  a.num = (*cur0).num
			    * comb (n + s, n)
			    * s;
			  a.den = (*cur0).den
			    * (n + 1);
			  reduce_ratio (&a);

			  cur = coef_q + pointer_npq (nmax, n, p, q);
			  b.sgn = (*cur).sgn;
			  b.num = (*cur).num;
			  b.den = (*cur).den;

			  add_ratio (a, b, cur);
			}
		    }
		}
	    }
	}
    }
}

void
YM_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q)
{
  int i;
  int n, p, q, s;

  ratio *cur;
  ratio *cur0;
  ratio a, b, c;

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
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  a.sgn = (*cur0).sgn;
			  a.num = (*cur0).num
			    * comb (n + s, n + 1)
			    * 2 * n;
			  a.den = (*cur0).den
			    * (n + 1) * (2 * n + 3);
			  reduce_ratio (&a);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  b.sgn = (*cur).sgn;
			  b.num = (*cur).num;
			  b.den = (*cur).den;

			  add_ratio (a, b, cur);
			}

		      /* Pnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  b.sgn = 1;
			  b.num = 2 * n + 1;
			  b.den = 2;
			  reduce_ratio (&b);

			  b.num *=
			    (n * s + 4) * (n + s)
			    - 2 * (n * s  + 1) * (n * s  + 1);
			  b.den *=
			    s * (2 * s - 1) * (n + s);
			  reduce_ratio (&b);

			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  b.sgn *= (*cur0).sgn;
			  b.num *= (*cur0).num;
			  b.den *= (*cur0).den;
			  reduce_ratio (&b);

			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (2 * n - 1);
			  b.den = (*cur0).den
			    * 2;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (4 * n * n - 1);
			  b.den = (*cur0).den
			    * 2 * (2 * s + 1);
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * (4 * n * n - 1);
			  b.den = (*cur0).den;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 1);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);

		      /* Qnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * s;
			  b.den = (*cur0).den;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num;
			  b.den = (*cur0).den
			    * n * s;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 1);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p
		    + pointer_npq (nmax, n, p, q);
		  a.sgn = (*cur0).sgn;
		  a.num = (*cur0).num;
		  a.den = (*cur0).den;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num;
		  b.den = (*cur).den;

		  add_ratio (a, b, cur);
		}
	    }
	}
    }
}

void
Z_problem (int nmax, ratio *coef_p, ratio *coef_v, ratio *coef_q)
{
  int i;
  int n, p, q, s;

  ratio *cur;
  ratio *cur0;
  ratio a, b, c;

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
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  cur = coef_q + pointer_npq (nmax, n, p, q);
		  (*cur).sgn = 1;
		  (*cur).num = 0;
		  (*cur).den = 1;

		  for (s=1; s <= q; s++)
		    {
		      /* q - s >= 0 */
		      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  a.sgn = (*cur0).sgn;
			  a.num = (*cur0).num
			    * comb (n + s, n + 2)
			    * 2 * n;
			  a.den = (*cur0).den
			    * (n + 1) * (2 * n + 3);
			  reduce_ratio (&a);

			  cur = coef_v + pointer_npq (nmax, n, p, q);
			  b.sgn = (*cur).sgn;
			  b.num = (*cur).num;
			  b.den = (*cur).den;

			  add_ratio (a, b, cur);
			}

		      /* Pnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
		      if (p - n + 1 >= 0)
			{
			  b.sgn = 1;
			  b.num = 2 * n + 1;
			  b.den = 2;
			  reduce_ratio (&b);

			  b.num *=
			    (n * s + 16) * (n + s)
			    - 2 * (n * s  + 4) * (n * s  + 1);
			  b.den *=
			    s * (2 * s - 1) * (n + s);
			  reduce_ratio (&b);

			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n + 1);
			  b.sgn *= (*cur0).sgn;
			  b.num *= (*cur0).num;
			  b.den *= (*cur0).den;
			  reduce_ratio (&b);

			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Pnpq : 2nd term Ps(q-s)(p-n-1) */
		      if (p - n - 1 >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n - 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (2 * n - 1);
			  b.den = (*cur0).den
			    * 2;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 3rd term Vs(q-s-2)(p-n+1) */
		      if (q - s - 2 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_v
			    + pointer_npq (nmax, s, q - s - 2, p - n + 1);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * n * (4 * n * n - 1);
			  b.den = (*cur0).den
			    * 2 * (2 * s + 1);
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}

		      /* Pnpq : 4th term Qs(q-s-1)(p-n+1) */
		      if (q - s - 1 >= 0
			  && p - n + 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n + 1);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * 2 * (4 * n * n - 1);
			  b.den = (*cur0).den;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 2);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_p + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);

		      /* Qnpq */
		      a.sgn = 1;
		      a.num = 0;
		      a.den = 1;
		      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
		      if (p - n>= 0
			  && q - s - 1 >= 0)
			{
			  cur0 = coef_q
			    + pointer_npq (nmax, s, q - s - 1, p - n);
			  b.sgn = (*cur0).sgn;
			  b.num = (*cur0).num
			    * s;
			  b.den = (*cur0).den;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      /* Qnpq : 2nd term Ps(q-s)(p-n) */
		      if (p - n >= 0)
			{
			  cur0 = coef_p
			    + pointer_npq (nmax, s, q - s, p - n);
			  b.sgn = - (*cur0).sgn;
			  b.num = (*cur0).num
			    * 2;
			  b.den = (*cur0).den
			    * n * s;
			  add_ratio (a, b, &c);

			  a.sgn = c.sgn;
			  a.num = c.num;
			  a.den = c.den;
			}
		      a.num *= comb (n + s, n + 2);
		      a.den *= n + 1;
		      reduce_ratio (&a);

		      cur = coef_q + pointer_npq (nmax, n, p, q);
		      b.sgn = (*cur).sgn;
		      b.num = (*cur).num;
		      b.den = (*cur).den;

		      add_ratio (a, b, cur);
		    }
		  /* Vnpq : 1st term Pnpq */
		  cur0 = coef_p
		    + pointer_npq (nmax, n, p, q);
		  a.sgn = (*cur0).sgn;
		  a.num = (*cur0).num;
		  a.den = (*cur0).den;

		  cur = coef_v + pointer_npq (nmax, n, p, q);
		  b.sgn = (*cur).sgn;
		  b.num = (*cur).num;
		  b.den = (*cur).den;

		  add_ratio (a, b, cur);
		}
	    }
	}
    }
}


double
combi (int n, int m)
{
  int i;
  double x, y, z;

  if(n < 0 || m < 0 || n < m)
    return (0.0);

  x = 1.0;
  for (i=1; i <= n; i++)
    {
      x *= (double)i;
    }
  y=1.0;
  for (i=1; i <= m; i++)
    {
      y *= (double)i;
    }
  z=1.0;
  for (i=1; i <= (n - m); i++)
    {
      z *= (double)i;
    }
  return (x / y / z);
}

int
comb (int n, int m)
{
  int i;
  ratio a;

  if(n < 0 || m < 0 || n < m)
    return (0);

  if ((n - m) < m)
    m = n - m;

  a.sgn = 1;
  a.num = 1;
  a.den = 1;
  for (i=1; i <= m; i++)
    {
      a.num *= (n + 1 - i);
      a.den *= (m + 1 - i);
      reduce_ratio (&a);
    }

  if (a.sgn != 1
      || a.den != 1)
    fprintf (stderr, "something is wrong (%d %d)\n", n, m);

  return (a.num);
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
