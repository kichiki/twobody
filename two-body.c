/* exact solution solver for 2 particles in Stokes flows using GMP library
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: two-body-JO-gmp.c,v 4.2 2005/01/22 18:42:56 ichiki Exp $
 *
 * References:
 * [JO-1984] D.J.Jeffrey and Onishi, J. Fluid Mech. 139 (1984) pp.261-290.
 *   XA, YA, YB, XC, YC, xa, ya, yb, xc, yc are given.
 * [J-1992] D.J.Jeffrey, Phys. Fluids A 4 (1992) pp.16-29.
 *   XG, YG, YH, XM, YM, ZM are given.
 * [JMB-1993] D.J.Jeffrey, J.F.Morris and J.F.Brady,
 *   Phys. Fluids A 5 (1993) pp.2317-2325.
 *   XP, XQ are given.
 * [KSB-2006] A.S.Khair, M.Swaroop and J.F.Brady,
 *   Phys. Fluids 18, 043102 (2006).
 *   TQ is given.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include <stdio.h>
#include <stdlib.h> /* malloc() free() */
#include <gmp.h>


/* function prototypes */
int
main (void);

void
XA (int nmax, int flag);
void
YA (int nmax, int flag);
void
YB (int nmax, int flag);
void
XC (int nmax, int flag);
void
YC (int nmax, int flag);
void
XG (int nmax, int flag);
void
YG (int nmax, int flag);
void
YH (int nmax, int flag);
void
XM (int nmax, int flag);
void
YM (int nmax, int flag);
void
ZM (int nmax, int flag);

void
XP (int nmax, int flag);
void
XQ (int nmax, int flag);

void
TQ (int nmax, int flag);


/** mobility functions **/
void
xa (int nmax, int flag);
void
ya (int nmax, int flag);
void
yb (int nmax, int flag);
void
xc (int nmax, int flag);
void
yc (int nmax, int flag);


/** recurrence relation solvers **/
void
X_p (int n, int p, int q, int mode, mpq_t coef_p);
void
X_v (int n, int p, int q, int mode, mpq_t coef_v);
void
X_q (int n, int p, int q, int mode, mpq_t coef_q);

void
Y_p (int n, int p, int q, int mode, mpq_t coef_p);
void
Y_v (int n, int p, int q, int mode, mpq_t coef_v);
void
Y_q (int n, int p, int q, int mode, mpq_t coef_q);

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

void
TQ_p (int n, int p, int q, mpq_t coef_p);
void
TQ_v (int n, int p, int q, mpq_t coef_v);

/** for mobility functions **/
void
x_u (int p, int q, int mode, mpq_t coef_u);
void
y_u (int p, int q, int mode, mpq_t coef_u);
void
y_q (int p, int q, int mode, mpq_t coef_q);
void
x_q (int p, int q, mpq_t coef_q);


/** utility functions **/
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

  nmax = 100;

  /* mobility functions */
  //xa (nmax, 0);
  //ya (nmax, 0);
  /*
  yb (nmax, 0);
  xc (nmax, 0);
  yc (nmax, 0);
  exit (0);
  */

  //xa (nmax, 1);
  //ya (nmax, 1);
  //yb (nmax, 1);
  //xc (nmax, 1);
  //yc (nmax, 1);

  //XA (nmax, 1);
  //YA (nmax, 1);
  //YB (nmax, 1);
  //XC (nmax, 1);
  //YC (nmax, 1);
  //XG (nmax, 1);
  //YG (nmax, 1);
  //YH (nmax, 1);
  //XM (nmax, 1);
  //YM (nmax, 1);
  //ZM (nmax, 1);
  //XP (nmax, 1);
  //XQ (nmax, 1);
  TQ (nmax, 1);

  return 0;
}


/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XA (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "XAf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  X_p (1, k - q, q, 10/* for XA */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XAf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YA (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "YAf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  Y_p (1, k - q, q, 10/* for YA (0, 1, 1) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YAf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YB (int nmax, int flag)
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
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k=0, mpq_set_ui (two_k, 2, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "YBf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 10/* for YB (0, 1, 1) */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      mpq_mul (coef_q, coef_q, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YBf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XC (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "XCf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  X_q (1, k - q, q, 1/* resistance */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  if (k%2 == 1)
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		  else
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XCf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YC (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "YCf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 11/* for YC (1, 0, 0) */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YCf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XG (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t three4;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three4);
  mpq_set_ui (three4, 3, 4);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "XGf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 10/* for X[AGP] */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (3/4) */
	      mpq_mul (coef_p, coef_p, three4);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XGf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YG (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t three4;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three4);
  mpq_set_ui (three4, 3, 4);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "YGf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 10/* for YG (0, 1, 1) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (3/4) */
	      mpq_mul (coef_p, coef_p, three4);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YGf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YH (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t mthree8;

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (mthree8);
  mpq_set_si (mthree8, -3, 8);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "YHf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 11/* for YH (1, 0, 0) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (-3/8) */
	      mpq_mul (coef_p, coef_p, mthree8);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YHf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (mthree8);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XM (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "XMf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 11/* for XM */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XMf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
YM (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "YMf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 0; q <= k; q ++)
	{
	  YM_p (2, k - q, q, 0, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "YMf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
ZM (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "ZMf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 0; q <= k; q ++)
	{
	  ZM_p (2, k - q, q, 0, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);

	      if (flag != 0)
		{
		  print_mpq (coef_p);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_p);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "ZMf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XP (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "XPf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 1; q <= k - 1; q ++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n = 1; n <= (q + 1) / 2; n++)
	    {
	      X_p (n, q - n, k - q - 1, 10/* for X[AGP] (1, 1, 1) */, b);

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
	    {
	      if (flag != 0)
		{
		  print_mpq (a);
		  fprintf (stdout, " l^%d ", q);
		}
	      else
		{
		  mpq_add (f, f, a);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XPf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

/* -- still something is wrong on f_8 (only) from JMB 1993
 * and question about the upper limit of "n"...
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
XQ (int nmax, int flag)
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
      if (flag != 0)
	{
	  fprintf (stdout, "XQf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 1; q <= k - 1; q ++)
	{
	  mpq_set_ui (a, 0, 1);
	  //for (n=1; n <= (q + 1) / 2; n++)
	  for (n = 1; n <= q; n++)
	    {
	      X_p (n, q - n, k - q - 1, 11/* for X[QM] (2, 2) */, b);

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
	    {
	      if (flag != 0)
		{
		  print_mpq (a);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, a);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "XQf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
TQ (int nmax, int flag)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "TQf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}

      for (q = 1; q <= k - 1; q ++)
	{
	  mpq_set_ui (a, 0, 1);
	  for (n = 1; n <= (q + 1) / 2; n++)
	    {
	      TQ_p (n, q - n, k - q - 1, b);

	      if (mpz_cmp_si (mpq_numref (b), 0))
		{
		  /* mul (two_k) */
		  mpq_mul (b, b, two_k);

		  mpq_add (a, a, b);
		}
	    }
	  if (mpz_cmp_si (mpq_numref (a), 0))
	    {
	      if (flag != 0)
		{
		  print_mpq (a);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q+1);
		    }
		}
	      else
		{
		  mpq_add (f, f, a);
		}
	    }
	}

      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "TQf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (two);
  mpq_clear (two_k);
}


/** mobility functions **/
/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
xa (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_u;
  mpq_t two;
  mpq_t two_k;

  mpq_init (f);
  mpq_init (coef_u);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "xaf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  x_u (k - q, q, 0/* mob for xa */, coef_u);

	  if (mpz_cmp_si (mpq_numref (coef_u), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_u, coef_u, two_k);
	      if (flag != 0)
		{
		  print_mpq (coef_u);
		  if (k%2 == 0) // k == even
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else // k == odd
		    {
		      fprintf (stdout, " l^%d ", q-1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_u);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "xaf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
ya (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_u;
  mpq_t two;
  mpq_t two_k;

  mpq_t minus1;

  mpq_init (f);
  mpq_init (coef_u);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (minus1);

  mpq_set_ui (two, 2, 1);
  mpq_set_si (minus1, -1, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "yaf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  y_u (k - q, q, 0/* for ya and yb*/, coef_u);

	  if (mpz_cmp_si (mpq_numref (coef_u), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_u, coef_u, two_k);
	      // another ad-hoc adjustment...
	      if (k != 0 && k%2 == 0)
		{
		  mpq_mul (coef_u, coef_u, minus1);
		}

	      if (flag != 0)
		{
		  print_mpq (coef_u);
		  if (k%2 == 0)
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		  else
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q-1);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_u);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "yaf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (minus1);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
yb (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  mpq_t three;
  mpq_t minus1;

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three);
  mpq_init (minus1);

  mpq_set_ui (two, 2, 1);
  mpq_set_ui (three, 3, 1); // for adjusting the result...
  mpq_set_si (minus1, -1, 1);
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k=0, mpq_set_ui (two_k, 2, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "ybf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  y_q (k - q, q, 0/* for ya and yb*/, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      mpq_mul (coef_q, coef_q, two_k);
	      // for adjusting the result by factor '1/3'...
	      mpq_div (coef_q, coef_q, three);
	      // another ad-hoc adjustment...
	      if (k%2 == 0)
		{
		  mpq_mul (coef_q, coef_q, minus1);
		}

	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  if (k%2 == 0)
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q-1);
		    }
		  else
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "ybf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three);
  mpq_clear (minus1);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
xc (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  //mpq_t two;
  //mpq_t two_k;

  mpq_t minus1;

  mpq_init (f);
  mpq_init (coef_q);
  //mpq_init (two);
  //mpq_init (two_k);

  mpq_init (minus1);

  //mpq_set_ui (two, 2, 1);
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  mpq_set_si (minus1, -1, 1);
  for (k=0/*, mpq_set_ui (two_k, 2, 1)*/;
       k <= nmax;
       k++/*, mpq_mul (two_k, two_k, two)*/)
    {
      if (flag != 0)
	{
	  fprintf (stdout, "xcf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  x_q (k - q, q, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      //mpq_mul (coef_q, coef_q, two_k);
	      if (k%2 == 1)
		{
		  // is this right?
		  mpq_mul (coef_q, coef_q, minus1);
		}

	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  if (k%2 == 1)
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q-2);
		    }
		  else
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "xcf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  //mpq_clear (two);
  //mpq_clear (two_k);

  mpq_clear (minus1);
}

/*
 * INPUT
 *  flag : 0 -- print only f_k (lambda=1)
 *         1 -- print each coef of l^q
 */
void
yc (int nmax, int flag)
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
  //mpq_set_ui (three, 3, 1); // for adjusting the result...
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k=0, mpq_set_ui (two_k, 2, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      if (flag != 0)
	{
	  fprintf (stdout, "ycf_%d = ", k);
	}
      else
	{
	  mpq_set_ui (f, 0, 1);
	}
      for (q = 0; q <= k; q ++)
	{
	  y_q (k - q, q, 1/* for yc */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      mpq_mul (coef_q, coef_q, two_k);
	      // for adjusting the result by factor '1/2'...
	      mpq_div (coef_q, coef_q, two);
	      if (flag != 0)
		{
		  print_mpq (coef_q);
		  if (k%2 == 1)
		    {
		      // is this right?
		      fprintf (stdout, " l^%d ", q-2);
		    }
		  else
		    {
		      fprintf (stdout, " l^%d ", q);
		    }
		}
	      else
		{
		  mpq_add (f, f, coef_q);
		}
	    }
	}
      if (flag != 0)
	{
	  fprintf (stdout, "\n");
	}
      else
	{
	  fprintf (stdout, "ycf [%d] = ", k);
	  print_mpq (f);
	  fprintf (stdout, "\n");
	}
    }

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}



/** recurrence relation solvers **/
/* calc coef Pnpq for X[AGP,MQ] and xa by Eqs.(JO-3.9) == (JO-8.17)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility function xa
 *             where P1pq = delta_0p delta_0q
 *   mode == 10 for resistance function X[AGP]
 *             where Pn00 = Vn00 = delta_1n
 *   mode == 11 for resistance function X[MQ]
 *             where Pn00 = Vn00 = delta_2n
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
X_p (int n, int p, int q, int mode, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];
  int inp = 0;
  int inv = 0;


  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // mobility for xa

      /* initial condition */
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.x_p.%d", mode);
      if (search_results (filename, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 10: // resistance function X[AGP]
      inp = 1;
      inv = 1;

      break;

    case 11: // resistance function X[MQ]
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  // for all resistance functions
  if (mode == 10
      || mode == 11)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inp != 0
	      && n == inp)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.X_p.%d.%d.%d", 0, inp, inv);
      if (search_results (filename, n, p, q, coef_p) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      X_p (s, q - s, p - n + 1, mode, b);

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
      X_p (s, q - s, p - n - 1, mode, b);

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
      X_v (s, q - s - 2, p - n + 1, mode, b);

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
/* calc coef Vnpq for X[AGP,MQ] and xa by Eqs.(JO-3.8) == (JO-8.15)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility function xa
 *             where P1pq = delta_0p delta_0q
 *   mode == 10 for resistance function X[AGP]
 *             where Pn00 = Vn00 = delta_1n
 *   mode == 11 for resistance function X[MQ]
 *             where Pn00 = Vn00 = delta_2n
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
X_v (int n, int p, int q, int mode, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];
  int inp = 0;
  int inv = 0;


  /* clear Vnpq */
  mpq_set_ui (coef_v, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // mobility for xa

      /* Vnpq for mobility has no initial condition */

      /* search buffer */
      sprintf (filename, "two-body.x_v.%d", mode);
      if (search_results (filename, n, p, q, coef_v) == 0)
	return; /* found */

      break;

    case 10: // resistance function X[AGP]
      inp = 1;
      inv = 1;

      break;

    case 11: // resistance function X[MQ]
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  // for all resistance functions
  if (mode == 10
      || mode == 11)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inv != 0
	      && n == inv)
	    {
	      mpq_set_ui (coef_v, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.X_v.%d.%d.%d", 0, inp, inv);
      if (search_results (filename, n, p, q, coef_v) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  X_p (n, p, q, mode, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      X_p (s, q - s, p - n - 1, mode, a);

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

/* calc coef Qnpq for XC and xc by Eqs.(JO-6.5b) and (JO-12.2)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility function xc
 *             where Q1pq = delta_0p delta_0q
 *   mode == 1 for resistance function XC
 *             where Qn00 = delta_1n
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
X_q (int n, int p, int q, int mode, mpq_t coef_q)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // mobility problem
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_q, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.x_q");
      if (search_results (filename, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    case 1: // resistance problem
      if (p == 0
	  && q == 0)
	{
	  /* initial condition */
	  if (n == 1)
	    {
	      mpq_set_ui (coef_q, 1, 1);
	    }
	  return;
	}
      /* necessary ?
      if (n > p + 1
	  || q == 0)
	{
	  mpq_set_ui (coef_q, 0, 1);
	  return;
	}
      */

      /* search buffer */
      sprintf (filename, "two-body.X_q.%d.%d.%d", 1, 0, 0);
      if (search_results (filename, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  mpq_init (a);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Qnpq : Qs(q-s-1)(p-n) */
      X_q (s, q - s - 1, p - n, mode, a);

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


/* calc coef Pnpq on y[agcgh] by Eq. (JO-9.14) == (JO-4.10)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *        == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
Y_p (int n, int p, int q, int mode, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // P1pq = delta_0p delta_0q and Q1pq = 0

      /* initial condition */
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  else
	    {
	      mpq_set_ui (coef_p, 0, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 1: // P1pq = 0 and Q1pq = delta_0p delta_0q

      /* initial condition */
      if (n == 1)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 10: // resistance functions Y[ABG]
      inq = 0;
      inp = 1;
      inv = 1;

      break;

    case 11: // resistance functions Y[CH]
      inq = 1;
      inp = 0;
      inv = 0;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inp != 0
	      && n == inp)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}
      if (n > p + 1
	  || q == 0)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.Y_p.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, n, p, q, coef_p) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);

  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      Y_p (s, q - s, p - n + 1, mode, b);
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
      Y_p (s, q - s, p - n - 1, mode, b);
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
      Y_v (s, q - s - 2, p - n + 1, mode, b);
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
      Y_q (s, q - s - 1, p - n + 1, mode, b);
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

/* calc coef Vnpq on y[abcgh] by Eq. (JO-9.13) == (JO-4.9)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *        == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
Y_v (int n, int p, int q, int mode, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Vnpq */
  mpq_set_ui (coef_v, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // P1pq = delta_0p delta_0q and Q1pq = 0
    case 1: // P1pq = 0 and Q1pq = delta_0p delta_0q

      /* Vnpq for mobility has no initial condition */

      /* search buffer */
      sprintf (filename, "two-body.y_v.%d", mode);
      if (search_results (filename, n, p, q, coef_v) == 0)
	return; /* found */

      break;

    case 10: // resistance functions Y[A]
      inq = 0;
      inp = 1;
      inv = 1;

      break;

    case 11: // resistance functions Y[]
      inq = 1;
      inp = 0;
      inv = 0;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inv != 0
	      && n == inv)
	    {
	      mpq_set_ui (coef_v, 1, 1);
	    }
	  return;
	}
      if (n > p + 1
	  || q == 0)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.Y_v.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, n, p, q, coef_v) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  Y_p (n, p, q, mode, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n+1) */
      Y_p (s, q - s, p - n - 1, mode, a);
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

/* calc coef Qnpq on y[agcgh] by Eq. (JO-9.15) == (JO-4.11)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *        == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
Y_q (int n, int p, int q, int mode, mpq_t coef_q)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);

  if (n < 1 // exclude n == 0, too.
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // P1pq = delta_0p delta_0q and Q1pq = 0

      /* initial condition */
      if (n == 1)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_q.%d", mode);
      if (search_results (filename, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    case 1: // P1pq = 0 and Q1pq = delta_0p delta_0q

      /* initial condition */
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_q, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_q.%d", mode);
      if (search_results (filename, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    case 10: // resistance functions Y[ABG]
      inq = 0;
      inp = 1;
      inv = 1;

      break;

    case 11: // resistance functions Y[CH]
      inq = 1;
      inp = 0;
      inv = 0;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inq != 0
	      && n == inq)
	    {
	      mpq_set_ui (coef_q, 1, 1);
	    }
	  return;
	}
      if (n > p + 1
	  || q == 0)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.Y_q.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, n, p, q, coef_q) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Qnpq */
      mpq_set_ui (a, 0, 1);
      /* Qnpq : 1st term Qs(q-s-1)(p-n) */
      Y_q (s, q - s - 1, p - n, mode, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul s */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-n) */
      Y_p (s, q - s, p - n, mode, b);
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
  if (n > p + 1
      || q == 0)
    {
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
  if (n > p + 1
      || q == 0)
    {
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
  if (n > p + 1
      || q == 0)
    {
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
  if (n > p + 1
      || q == 0)
    {
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
  if (n > p + 1
      || q == 0)
    {
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
  if (n > p + 1
      || q == 0)
    {
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


/* calc coef Pnpq for TQ
 * INPUT
 *   n, p, q : (int)
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
TQ_p (int n, int p, int q, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];


  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (n < 0 // include n == 0 here...
      || p < 0
      || q < 0)
    {
      return;
    }

  /* initial condition */
  if (n == 0)
    {
      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.TQ_p");
  if (search_results (filename, n, p, q, coef_p) == 0)
    return; /* found */

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 0; s <= q; s ++) // start from 0!!
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      TQ_p (s, q - s, p - n + 1, b);

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
      TQ_p (s, q - s, p - n - 1, b);

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
      TQ_v (s, q - s - 2, p - n + 1, b);

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
/* calc coef Vnpq for TQ
 * INPUT
 *   n, p, q : (int)
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
TQ_v (int n, int p, int q, mpq_t coef_v)
{
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];


  /* clear Vnpq */
  mpq_set_ui (coef_v, 0, 1);

  if (n < 0 // include n == 0 here...
      || p < 0
      || q < 0)
    {
      return;
    }

  /* initial condition */
  if (n == 0)
    {
      if (p == 0 && q == 0)
	{
	  mpq_set_si (coef_v, -1, 1);
	}
      return;
    }

  /* search buffer */
  sprintf (filename, "two-body.TQ_v");
  if (search_results (filename, n, p, q, coef_v) == 0)
    return; /* found */

  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  TQ_p (n, p, q, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      TQ_p (s, q - s, p - n - 1, a);

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


/** for mobility functions **/
/* calc coef Upq on xa by Eq. (JO-8.16)
 * INPUT
 *   p, q : (int)
 *   mode == 0 for mobility function xa
 *             where P1pq = delta_0p delta_0q
 *   mode == 10 for resistance function X[AGP]
 *             where Pn00 = Vn00 = delta_1n
 *   mode == 11 for resistance function X[MQ]
 *             where Pn00 = Vn00 = delta_2n
 * OUTPUT
 *   coef_u : (mpq_t)
 */
void
x_u (int p, int q, int mode, mpq_t coef_u)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;


  /* clear Upq */
  mpq_set_ui (coef_u, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  /* ad-hoc for U_{00} */
  if (p == 0 && q == 0)
    {
      mpq_set_ui (coef_u, 1, 1);
      return;
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Upq */
      mpq_set_ui (a, 0, 1);
      /* Upq : 1st term Ps(q-s)p */
      X_p (s, q - s, p, 0/* mob for xa */, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (3/4) */
	  mpq_set_si (tmp, 3, 4);
	  mpq_mul (b, b, tmp);
	  /* div (2s-1) */
	  mpq_set_si (tmp, (2 * s - 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Upq : 2nd term Ps(q-s)(p-2) */
      X_p (s, q - s, p - 2, 0/* mob for xa */, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -(1/4) */
	  mpq_set_si (tmp, - 1, 4);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      /* Upq : 3rd term Vs(q-s-2)p */
      X_v (s, q - s - 2, p, 0/* mob for xa */, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -(3/4) */
	  mpq_set_si (tmp, - 3, 4);
	  mpq_mul (b, b, tmp);
	  /* div (2s+1) */
	  mpq_set_si (tmp, (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul -(s+1) */
	  mpq_set_si (tmp, -(s + 1), 1);
	  mpq_mul (a, a, tmp);

	  mpq_add (coef_u, coef_u, a);
	}
    }

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Upq for ya by Eq. (JO-9.16)
 * INPUT
 *   p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *        == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[A]
 *             where Pn00 = Vn00 = delta_1n
 *   mode == 11 for resistance functions Y[]
 *             where Pn00 = Vn00 = delta_2n
 * OUTPUT
 *   coef_u : (mpq_t)
 */
void
y_u (int p, int q, int mode, mpq_t coef_u)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;


  /* clear Upq */
  mpq_set_ui (coef_u, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  /* ad-hoc for U_{00} */
  if (p == 0 && q == 0)
    {
      if (mode == 0)
	{
	  mpq_set_ui (coef_u, 1, 1);
	}
      else
	{
	  mpq_set_ui (coef_u, 0, 1);
	}
      return;
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Upq */
      mpq_set_ui (a, 0, 1);
      /* Upq : 1st term Ps(q-s)p */
      Y_p (s, q - s, p, mode, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (3/4s) */
	  mpq_set_si (tmp, 3, (4 * s));
	  mpq_mul (b, b, tmp);
	  /* mlu (2-s) */
	  mpq_set_si (tmp, (2 - s), 1);
	  mpq_mul (b, b, tmp);
	  /* div (2s-1) */
	  mpq_set_si (tmp, (2 * s - 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Upq : 2nd term Ps(q-s)(p-2) */
      Y_p (s, q - s, p - 2, mode, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (1/4) */
	  mpq_set_si (tmp, 1, 4);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      /* Upq : 3rd term Vs(q-s-2)p */
      Y_v (s, q - s - 2, p, mode, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (3/4) */
	  mpq_set_si (tmp, 3, 4);
	  mpq_mul (b, b, tmp);
	  /* div (2s+1) */
	  mpq_set_si (tmp, (2 * s + 1), 1);
	  mpq_div (b, b, tmp);

	  mpq_add (a, a, b);
	}

      /* Upq : 4th term Qs(q-s-1)p */
      Y_q (s, q - s - 1, p, mode, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -1 */
	  mpq_set_si (tmp, -1, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (s+1, 2) */
	  comb (tmp, s + 1, 2);
	  mpq_mul (a, a, tmp);

	  /* mul -1 */
	  mpq_set_si (tmp, -1, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (coef_u, coef_u, a);
	}
    }

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Qpq for ya by Eq. (JO-9.17)
 * INPUT
 *   p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *        == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[A]
 *             where Pn00 = Vn00 = delta_1n
 *   mode == 11 for resistance functions Y[]
 *             where Pn00 = Vn00 = delta_2n
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
y_q (int p, int q, int mode, mpq_t coef_q)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;


  /* clear Qpq */
  mpq_set_ui (coef_q, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  /* ad-hoc for Q_{00} */
  if (p == 0 && q == 0)
    {
      if (mode == 0)
	{
	  mpq_set_ui (coef_q, 0, 1);
	}
      else
	{
	  mpq_set_ui (coef_q, 1, 1);
	}
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);

  for (s = 1; s <= q; s ++)
    {
      /* Qnpq */
      mpq_set_ui (a, 0, 1);
      /* Qnpq : 1st term Qs(q-s-1)(p-1) */
      Y_q (s, q - s - 1, p - 1, mode, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul (s / 2) */
	  mpq_set_si (tmp, s, 2);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-1) */
      Y_p (s, q - s, p - 1, mode, b);
      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul -(3 / 4s) */
	  mpq_set_si (tmp, - 3, 4 * s);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (s+1,2) */
	  comb (tmp, s + 1, 2);
	  mpq_mul (a, a, tmp);

	  /* mul -1 */
	  mpq_set_si (tmp, -1, 1);
	  mpq_mul (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}


/* calc coef Qpq for xc by Eq. (JO-12.3)
 * INPUT
 *   p, q : (int)
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
x_q (int p, int q, mpq_t coef_q)
{
  int s;

  mpq_t a;
  mpq_t tmp;


  /* clear Qpq */
  mpq_set_ui (coef_q, 0, 1);

  if (p < 0
      || q < 0)
    {
      /* clear Qnpq */
      mpq_set_ui (coef_q, 0, 1);
      return;
    }

  /* ad-hoc for Q_{00} */
  if (p == 0 && q == 0)
    {
      mpq_set_ui (coef_q, 1, 1);
      return;
    }

  mpq_init (a);
  mpq_init (tmp);

  /* clear Qnpq */
  mpq_set_ui (coef_q, 0, 1);
  for (s = 1; s <= q; s ++)
    {
      /* Qnpq : Qs(q-s-1)(p-1) */
      X_q (s, q - s - 1, p - 1, 0/* mobility */, a);

      if (mpz_cmp_si (mpq_numref (a), 0))
	{
	  /* mul comb (s+1) */
	  mpq_set_si (tmp, (s + 1), 1);
	  mpq_mul (a, a, tmp);
	  /* mul -(s / 2) */
	  mpq_set_si (tmp, -s, 2);
	  mpq_mul (a, a, tmp);

	  mpq_add (coef_q, coef_q, a);
	}
    }

  mpq_clear (a);
  mpq_clear (tmp);

  mpq_canonicalize (coef_q);
}


/** utility functions **/
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
  char string [1024];
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
