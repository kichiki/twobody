/* exact solution solver for 2 particles in Stokes flows using GMP library
 * Copyright (C) 1999-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: two-body.c,v 5.3 2006/11/11 05:14:59 ichiki Exp $
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
#include <string.h>
#include <gmp.h>


/* function prototypes */

/** coefficients for functions **/
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
fprint_mpq (FILE *out, mpq_t x);
int
search_results (char *file, int n, int p, int q, mpq_t coef);
void
append_result (char *file, int n, int p, int q, mpq_t coef);


int
main (int argc, char** argv)
{
  int i;
  int nmax;
  int mode;
  int func;


  /* default values */
  nmax = 11;
  mode = 0;
  func = 0;

  /* option analysis */
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-n") == 0 ||
	  strcmp (argv [i], "--nmax") == 0)
	{
	  nmax = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-m") == 0 ||
	  strcmp (argv [i], "--mode") == 0)
	{
	  mode = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-f") == 0 ||
	       strcmp (argv [i], "--func") == 0)
	{
	  if (strcmp (argv [++i], "XA") == 0)
	    {
	      func =  1; // XA
	    }
	  else if (strcmp (argv [i], "YA") == 0)
	    {
	      func =  2; // YA
	    }
	  else if (strcmp (argv [i], "YB") == 0)
	    {
	      func =  3; // YB
	    }
	  else if (strcmp (argv [i], "XC") == 0)
	    {
	      func =  4; // XC
	    }
	  else if (strcmp (argv [i], "YC") == 0)
	    {
	      func =  5; // YC
	    }
	  else if (strcmp (argv [i], "XG") == 0)
	    {
	      func =  6; // XG
	    }
	  else if (strcmp (argv [i], "YG") == 0)
	    {
	      func =  7; // YG
	    }
	  else if (strcmp (argv [i], "YH") == 0)
	    {
	      func =  8; // YH
	    }
	  else if (strcmp (argv [i], "XM") == 0)
	    {
	      func =  9; // XM
	    }
	  else if (strcmp (argv [i], "YM") == 0)
	    {
	      func = 10; // YM
	    }
	  else if (strcmp (argv [i], "ZM") == 0)
	    {
	      func = 11; // ZM
	    }
	  else if (strcmp (argv [i], "XP") == 0)
	    {
	      func = 12; // XP
	    }
	  else if (strcmp (argv [i], "XQ") == 0)
	    {
	      func = 13; // XQ
	    }
	  else if (strcmp (argv [i], "TQ") == 0)
	    {
	      func = 14; // TQ
	    }
	  else if (strcmp (argv [i], "xa") == 0)
	    {
	      func = 15; // xa
	    }
	  else if (strcmp (argv [i], "ya") == 0)
	    {
	      func = 16; // ya
	    }
	  else if (strcmp (argv [i], "yb") == 0)
	    {
	      func = 17; // yb
	    }
	  else if (strcmp (argv [i], "xc") == 0)
	    {
	      func = 18; // xc
	    }
	  else if (strcmp (argv [i], "yc") == 0)
	    {
	      func = 19; // yc
	    }
	}
      else
	{
	  fprintf (stderr, "$Id$\n");
	  fprintf (stderr, "USAGE\n");
	  fprintf (stderr, "%s [OPTIONS]\n", argv [0]);
	  fprintf (stderr, "\t-h or --help : show this message.\n");
          fprintf (stderr, "\t-n or --nmax : max order\n");
          fprintf (stderr, "\t-m or --mode : output format\n");
	  fprintf (stderr, "\t\t0 plain txt for f_k with lambda=1\n");
          fprintf (stderr, "\t\t1 plain txt each coef of l^q\n");
          fprintf (stderr, "\t\t2 C source\n");
          fprintf (stderr, "\t\t3 latex source for f_k with lambda=1\n");
          fprintf (stderr, "\t\t4 latex source for each coef of l^q\n");
          fprintf (stderr, "\t-f or --func : function name to calculate;\n");
	  fprintf (stderr, "\t\tXA, YA, YB, XC, YC,\n");
	  fprintf (stderr, "\t\tXG, YG, YH, XM, YM, ZM,\n");
	  fprintf (stderr, "\t\tXP, XQ, TQ,\n");
	  fprintf (stderr, "\t\txa, ya, yb, xc, yc.\n");
	  exit (1);
	}
    }
  
  switch (func)
    {
    case 1:
      XA (nmax, mode);
      break;
    case 2:
      YA (nmax, mode);
      break;
    case 3:
      YB (nmax, mode);
      break;
    case 4:
      XC (nmax, mode);
      break;
    case 5:
      YC (nmax, mode);
      break;
    case 6:
      XG (nmax, mode);
      break;
    case 7:
      YG (nmax, mode);
      break;
    case 8:
      YH (nmax, mode);
      break;
    case 9:
      XM (nmax, mode);
      break;
    case 10:
      YM (nmax, mode);
      break;
    case 11:
      ZM (nmax, mode);
      break;
    case 12:
      XP (nmax, mode);
      break;
    case 13:
      XQ (nmax, mode);
      break;
    case 14:
      TQ (nmax, mode);
      break;
    case 15:
      xa (nmax, mode);
      break;
    case 16:
      ya (nmax, mode);
      break;
    case 17:
      yb (nmax, mode);
      break;
    case 18:
      xc (nmax, mode);
      break;
    case 19:
      yc (nmax, mode);
      break;
    default:
      fprintf (stderr, "invalid function\n");
      exit (1);
      break;
    }

  return 0;
}


/** I/O utility routines **/

/* for C source codes */
void
print_mpq_double (mpq_t x)
{
  double d;

  d = mpq_get_d (x);
  fprintf (stdout, "%.15e", d);
}

void
print_C_header (const char * label)
{
  fprintf (stdout, "void twobody_%s (int n, double l, double * f)\n",
	   label);
  fprintf (stdout, "{\n");
}

/*
 * INPUT
 *  q :
 *  q0: order of q in the last term (0 for the first term)
 */
void
print_C_coef_q (mpq_t coef, int q, int q0)
{
  fprintf (stdout, "    + ");

  int i;
  for (i = q0; i < q; i ++)
    {
      fprintf (stdout, "l * ");
    }

  fprintf (stdout, "(");

  print_mpq_double (coef);

  fprintf (stdout, " // q = %d\n", q);
}

void
print_C_close_q (int nq, int k)
{
  fprintf (stdout, "    ");
  int i;
  for (i = 0; i < nq; i ++)
    {
      fprintf (stdout, ")");
    }
  if (nq == 0)
    {
      fprintf (stdout, "0.0;\n");
    }
  else
    {
      fprintf (stdout, ";\n");
    }
  fprintf (stdout, "  if (n == %d) return;\n", k);
}

void
print_C_footer (void)
{
  fprintf (stdout, "}\n\n");
}


/* for text */
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


void
print_text_coef_q (mpq_t x, int q, int nq)
{
  if (nq > 0 && mpz_sgn (mpq_numref (x)) == 1)
    {
      fprintf (stdout, " +");
    }
  else
    {
      fprintf (stdout, " ");
    }
  print_mpq (x);

  if (q == 1)
    {
      fprintf (stdout, " l");
    }
  else
    {
      fprintf (stdout, " l^%d", q);
    }
}

void
print_text_fk (mpq_t f, int k, const char * label)
{
  fprintf (stdout, "%sf [%d] = ", label, k);
  print_mpq (f);
  fprintf (stdout, "\n");
}


/* for latex source */
void
print_mpz_formed (mpz_t x)
{
  char buf  [1024];
  char buf2 [1024];

  int len;
  int i, j;

  gmp_snprintf (buf, 1024, "%Zd", x);

  len = strlen (buf);
  j = len + 2 * ((len-1)/3);
  buf2 [j] = '\0';
  j --;
  for (i = 0; i < (len - 3); i += 3, j -= 5)
    {
      buf2 [j  ] = buf [len-1-i  ];
      buf2 [j-1] = buf [len-1-i-1];
      buf2 [j-2] = buf [len-1-i-2];
      buf2 [j-3] = ' ';
      buf2 [j-4] = '\\';
    }
  for (; i < len; i ++, j --)
    {
      buf2 [j] = buf [len-1-i];
    }
  if (j != -1)
    {
      fprintf (stderr, "something is wrong...\n");
      fprintf (stderr, "check %s, %s\n", buf, buf2);
      exit (1);
    }

  fprintf (stdout, "%s", buf2);
}

/*
 * INPUT
 *  flag_plus : 0, don't print "+" for positive x.
 *              1, print "+" for positive x.
 */
void
print_mpq_formed (mpq_t x, int flag_plus)
{
  mpz_t tmp;
  mpz_init (tmp);

  if (flag_plus == 1
      && mpz_sgn (mpq_numref (x)) == 1)
    {
      fprintf (stdout, "+");
    }
  mpz_set (tmp, mpq_numref (x));
  if (mpz_sgn (mpq_numref (x)) == -1) // negative
    {
      mpz_neg (tmp, tmp);
      fprintf (stdout, "-");
    }
  if (mpz_cmp_ui (mpq_denref (x), 1))
    {
      /* fraction */
      fprintf (stdout, "\\frac{");
      //print_mpz_formed (mpq_numref (x));
      print_mpz_formed (tmp);
      fprintf (stdout, "}{");
      print_mpz_formed (mpq_denref (x));
      fprintf (stdout, "}");
    }
  else
    {
      //print_mpz_formed (mpq_numref (x));
      print_mpz_formed (tmp);
    }
  mpz_clear (tmp);
}

void
print_latex_coef_q (mpq_t x, int q, int nq)
{
  if (nq != 0 && nq % 4 == 0)
    {
      fprintf (stdout, "  \\nonumber\\\\\n  &&\n");
    }
  fprintf (stdout, "  ");
  if (nq == 0) print_mpq_formed (x, 0); // no "+"
  else         print_mpq_formed (x, 1); // with "+"

  if (q == 0)
    {
      fprintf (stdout, "\n");
    }
  else if (q == 1)
    {
      fprintf (stdout, " \\lambda\n");
    }
  else
    {
      fprintf (stdout, " \\lambda^{%d}\n", q);
    }
}

void
print_latex_fk (mpq_t f, int k, const char * label)
{
  fprintf (stdout, "\\begin{equation}\n  f^{\\rm %s}_{%d} = ", label, k);
  print_mpq_formed (f, 0);
  fprintf (stdout, "\n\\end{equation}\n");
}

/** top-level I/O routines for C source codes and latex source;
 *  flag : 0 -- print f_k with lambda=1
 *         1 -- print each coef of l^q
 *         2 -- print C source
 *         3 -- print latex source for f_k with lambda=1
 *         4 -- print latex source for each coef of l^q
 **/
/* there are 5 stages
 * 1) output overall header for the function
 * 2) some header statement for k
 * 3) output the coeficieint for q
 * 4) some footer statement at the end of k
 * 5) overall footer for the function
 */

/* 1 : output overall header for the function */
void
print_fk_header (int flag, const char * label)
{
  if (flag == 2) // C source
    {
      print_C_header (label);
    }
}

/* 2 : some header statement for k */
void
print_fk_header_k (int flag, int k, const char * label)
{
  if (flag == 1) // text for each coef of l^q
    {
      fprintf (stdout, "%sf [%d] = ",
	       label, k);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      fprintf (stdout, "\\begin{eqnarray}\n  f^{\\rm %s}_{%d} &=&\n",
	       label, k);
    }
  else if (flag == 2) // C source
    {
      fprintf (stdout, "  f [%d] =\n", k);
    }
}

/* 3 : output the coeficieint for q */
void
print_fk_q (int flag, mpq_t coef, mpq_t f, int q, int q0, int nq)
{
  if (flag == 1) // text for each coef of l^q
    {
      print_text_coef_q (coef, q, nq);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      print_latex_coef_q (coef, q, nq);
    }
  else if (flag == 2) // C source
    {
      print_C_coef_q (coef, q, q0);
    }
}

/* 4 : some footer statement at the end of k */
void
print_fk_footer_k (int flag, int k, mpq_t f, int nq, const char * label)
{
  if (flag == 0) // text for f_k with lambda=1
    {
      print_text_fk (f, k, label);
    }
  else if (flag == 1) // text for each coef of l^q
    {
      fprintf (stdout, "\n");
    }
  else if (flag == 3) // latex for f_k with lambda=1
    {
      print_latex_fk (f, k, label);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      if (nq == 0)
	{
	  if (nq == 0)
	    {
	      fprintf (stdout, "  0\n");
	    }
	}
      fprintf (stdout, "\\end{eqnarray}\n");
    }
  else if (flag == 2) // C source
    {
      print_C_close_q (nq, k);
    }
}

/* 5 : overall footer for the function */
void
print_fk_footer (int flag)
{
  if (flag == 2) // C source
    {
      print_C_footer ();
    }
}


/** coefficients for functions **/
/*  flag : 0 -- print f_k with lambda=1
 *         1 -- print each coef of l^q
 *         2 -- print C source
 *         3 -- print latex source for f_k with lambda=1
 *         4 -- print latex source for each coef of l^q
 */
void
XA (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XA");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XA");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  X_p (1, k - q, q, 10/* for XA */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XA");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YA (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YA");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YA");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_p (1, k - q, q, 10/* for YA (0, 1, 1) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YA");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YB (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YB");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YB");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 10/* for YB (0, 1, 1) */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      mpq_mul (coef_q, coef_q, two_k);
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, coef_q, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YB");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XC (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XC");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XC");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  X_q (1, k - q, q, 1/* resistance */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 1)
		{
		  // is this right?
		  print_fk_q (flag, coef_q, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      else
		{
		  print_fk_q (flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XC");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YC (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YC");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YC");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_q (1, k - q, q, 11/* for YC (1, 0, 0) */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_q, coef_q, two_k);
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_q, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YC");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XG (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t three4;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XG");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XG");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 10/* for X[AGP] */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (3/4) */
	      mpq_mul (coef_p, coef_p, three4);

	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XG");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

void
YG (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t three4;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YG");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YG");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 10/* for YG (0, 1, 1) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (3/4) */
	      mpq_mul (coef_p, coef_p, three4);

	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YG");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

void
YH (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  mpq_t mthree8;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YH");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YH");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 11/* for YH (1, 0, 0) */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      /* mul (-3/8) */
	      mpq_mul (coef_p, coef_p, mthree8);

	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YH");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (mthree8);
}

void
XM (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XM");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  X_p (2, k - q, q, 11/* for XM */, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YM (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "YM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "YM");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Y_p (2, k - q, q, 12/* for YM (0, 2, 2)*/, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "YM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ZM (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "ZM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "ZM");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  ZM_p (2, k - q, q, 0, 2, 2, coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "ZM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XP (int nmax, int flag)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t tmp;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XP");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XP");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
		  /* mul (3 / 2) */
		  mpq_set_ui (tmp, 3, 2);
		  mpq_mul (b, b, tmp);

		  mpq_add (a, a, b);
		}
	    }
	  if (mpz_cmp_si (mpq_numref (a), 0))
	    {
	      mpq_add (f, f, a);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (flag, a, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XP");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

/* -- still something is wrong on f_8 (only) from JMB 1993
 * and question about the upper limit of "n"...
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

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "XQ");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "XQ");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
		  /* mul (5 / 2) */
		  mpq_set_ui (tmp, 5, 2);
		  mpq_mul (b, b, tmp);

		  mpq_add (a, a, b);
		}
	    }
	  if (mpz_cmp_si (mpq_numref (a), 0))
	    {
	      mpq_add (f, f, a);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, a, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, a, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "XQ");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
TQ (int nmax, int flag)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "TQ");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "TQ");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
	      mpq_add (f, f, a);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, a, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, a, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "TQ");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (two);
  mpq_clear (two_k);
}


/** mobility functions **/
void
xa (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_u;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "xa");

  mpq_init (f);
  mpq_init (coef_u);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  for (k=0, mpq_set_ui (two_k, 1, 1);
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "xa");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  x_u (k - q, q, 0/* mob for xa */, coef_u);

	  if (mpz_cmp_si (mpq_numref (coef_u), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_u, coef_u, two_k);
	      mpq_add (f, f, coef_u);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (flag, coef_u, f, q, q0, nq);
		  q0 = q;
		}
	      else // k == odd
		{
		  print_fk_q (flag, coef_u, f, q-1, q0, nq);
		  q0 = q - 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "xa");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ya (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_u;
  mpq_t two;
  mpq_t two_k;

  mpq_t minus1;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "ya");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "ya");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
	      mpq_add (f, f, coef_u);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0)
		{
		  print_fk_q (flag, coef_u, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (flag, coef_u, f, q-1, q0, nq);
		  q0 = q-1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "ya");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (minus1);
}

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

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "yb");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "yb");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0)
		{
		  // is this right?
		  print_fk_q (flag, coef_q, f, q-1, q0, nq);
		  q0 = q-1;
		}
	      else
		{
		  print_fk_q (flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "yb");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three);
  mpq_clear (minus1);
}

void
xc (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  //mpq_t two;
  //mpq_t two_k;

  mpq_t minus1;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "xc");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "xc");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
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
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 1)
		{
		  // is this right?
		  print_fk_q (flag, coef_q, f, q-2, q0, nq);
		  q0 = q-2;
		}
	      else
		{
		  print_fk_q (flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "xc");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  //mpq_clear (two);
  //mpq_clear (two_k);

  mpq_clear (minus1);
}

void
yc (int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (flag, "yc");

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
      /* 2 : some header statement for k */
      print_fk_header_k (flag, k, "yc");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  y_q (k - q, q, 1/* for yc */, coef_q);

	  if (mpz_cmp_si (mpq_numref (coef_q), 0))
	    {
	      /* mul 2^k and another 2 */
	      mpq_mul (coef_q, coef_q, two_k);
	      // for adjusting the result by factor '1/2'...
	      mpq_div (coef_q, coef_q, two);
	      mpq_add (f, f, coef_q);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 1)
		{
		  // is this right?
		  print_fk_q (flag, coef_q, f, q-2, q0, nq);
		  q0 = q-2;
		}
	      else
		{
		  print_fk_q (flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (flag, k, f, nq, "yc");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (flag);

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


/* calc coef Pnpq on Y[ABCGHM] and y[abc] by Eqs. (JO-4.10) == (JO-9.14)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *   mode == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 *   mode == 12 for resistance functions YM
 *             where Pn00 = Vn00 = delta_2n, Qn00 = 0
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

    case 12: // resistance functions YM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11
      || mode == 12)
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

/* calc coef Vnpq on Y[ABCGHM] and y[abc] by Eqs. (JO-4.9) == (JO-9.13)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *   mode == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 *   mode == 12 for resistance functions YM
 *             where Pn00 = Vn00 = delta_2n, Qn00 = 0
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

    case 12: // resistance functions YM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11
      || mode == 12)
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

/* calc coef Qnpq on Y[ABCGHM] and y[abc] by Eqs. (JO-4.11) == (JO-9.15)
 * note that Qnpq in J-1992 is (2/3)Qnpq in JO-1984.
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions ya and yb
 *             where P1pq = delta_0p delta_0q
 *             and   Q1pq = 0
 *   mode == 1 for mobility functions yc
 *             where P1pq = 0
 *             and   Q1pq = delta_0p delta_0q
 *   mode == 10 for resistance functions Y[ABG]
 *             where Pn00 = Vn00 = delta_1n, Qn00 = 0
 *   mode == 11 for resistance functions Y[CH]
 *             where Pn00 = Vn00 = 0, Qn00 = delta_1n
 *   mode == 12 for resistance functions YM
 *             where Pn00 = Vn00 = delta_2n, Qn00 = 0
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

    case 12: // resistance functions YM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }


  // for all resistance functions
  if (mode == 10
      || mode == 11
      || mode == 12)
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
