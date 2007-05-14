/* exact solution solver for 2 particles in Stokes flows using GMP library
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: two-body.c,v 5.6 2007/04/12 20:19:59 ichiki Exp $
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


struct table_cache {
  int number;
  int *n;
  int *p;
  int *q;
  mpq_t *coef;
};

struct table_cache *cache_p = NULL;
struct table_cache *cache_v = NULL;
struct table_cache *cache_q = NULL;
int cache_max; // max size of cache
int cache_i; // current index


/* function prototypes */

void
print_mpq (FILE *out, mpq_t x);


/** coefficients for functions **/
void XA (FILE *out, int n0, int nmax, int flag);
void YA (FILE *out, int n0, int nmax, int flag);
void YB (FILE *out, int n0, int nmax, int flag);
void XC (FILE *out, int n0, int nmax, int flag);
void YC (FILE *out, int n0, int nmax, int flag);
void XG (FILE *out, int n0, int nmax, int flag);
void YG (FILE *out, int n0, int nmax, int flag);
void YH (FILE *out, int n0, int nmax, int flag);
void XM (FILE *out, int n0, int nmax, int flag);
void YM (FILE *out, int n0, int nmax, int flag);
void ZM (FILE *out, int n0, int nmax, int flag);

void XP (FILE *out, int n0, int nmax, int flag);
void XQ (FILE *out, int n0, int nmax, int flag);

void TQ (FILE *out, int n0, int nmax, int flag);


/** mobility functions **/
void xa (FILE *out, int n0, int nmax, int flag);
void ya (FILE *out, int n0, int nmax, int flag);
void yb (FILE *out, int n0, int nmax, int flag);
void xc (FILE *out, int n0, int nmax, int flag);
void yc (FILE *out, int n0, int nmax, int flag);

void xg (FILE *out, int n0, int nmax, int flag);
void yg (FILE *out, int n0, int nmax, int flag);
void yh (FILE *out, int n0, int nmax, int flag);
void xm (FILE *out, int n0, int nmax, int flag);
void ym (FILE *out, int n0, int nmax, int flag);
void zm (FILE *out, int n0, int nmax, int flag);


/** recurrence relation solvers **/
void X_p (int n, int p, int q, int mode, mpq_t coef_p);
void X_v (int n, int p, int q, int mode, mpq_t coef_v);
void X_q (int n, int p, int q, int mode, mpq_t coef_q);

void Y_p (int n, int p, int q, int mode, mpq_t coef_p);
void Y_v (int n, int p, int q, int mode, mpq_t coef_v);
void Y_q (int n, int p, int q, int mode, mpq_t coef_q);

void Z_p (int n, int p, int q, int mode, mpq_t coef_p);
void Z_v (int n, int p, int q, int mode, mpq_t coef_v);
void Z_q (int n, int p, int q, int mode, mpq_t coef_q);

void TQ_p (int n, int p, int q, mpq_t coef_p);
void TQ_v (int n, int p, int q, mpq_t coef_v);

/** for mobility functions **/
void x_u (int p, int q, int mode, mpq_t coef_u);
void y_u (int p, int q, int mode, mpq_t coef_u);
void y_q (int p, int q, int mode, mpq_t coef_q);
void x_q (int p, int q, mpq_t coef_q);

void x_E (int p, int q, int mode, mpq_t coef_p);
void y_E (int p, int q, int mode, mpq_t coef_p);
void z_E (int p, int q, int mode, mpq_t coef_p);


/** utility functions **/
int
comb (mpq_t comb, int n, int m);

/** table handling routines **/
/* table cache */
struct table_cache *
table_cache_init (void);
void
table_cache_free (struct table_cache *c);
void
table_cache_append (struct table_cache *c,
		    int n, int p, int q, mpq_t coef);
void
table_cache_set_i (struct table_cache *c,
		   int i,
		   int n, int p, int q, mpq_t coef);

int
search_results (char *file, struct table_cache *cache,
		int n, int p, int q, mpq_t coef);
void
append_result (char *file, struct table_cache *cache,
	       int n, int p, int q, mpq_t coef);



int
main (int argc, char** argv)
{
  extern struct table_cache *cache_p;
  extern struct table_cache *cache_v;
  extern struct table_cache *cache_q;
  extern int cache_max; // max size of cache
  extern int cache_i; // current index

  int i;
  int flag_out;
  char outfile [256];
  FILE *out = NULL;
  int nmax;
  int n0;
  int mode;
  int func;
  int flag_cache;


  /* default values */
  flag_out = 0;
  nmax = 11;
  n0   = 0;
  mode = 0;
  func = 0;
  flag_cache = 0;

  /* option analysis */
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-o") == 0 ||
	  strcmp (argv [i], "--out") == 0)
	{
	  strcpy (outfile, argv [++i]);
	  flag_out = 1;
	}
      else if (strcmp (argv [i], "-n") == 0 ||
	  strcmp (argv [i], "--nmax") == 0)
	{
	  nmax = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-n0") == 0 ||
	  strcmp (argv [i], "--n0") == 0)
	{
	  n0 = atoi (argv [++i]);
	}
      else if (strcmp (argv [i], "-c") == 0 ||
	  strcmp (argv [i], "--cache") == 0)
	{
	  flag_cache = 1;
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
	  else if (strcmp (argv [i], "xg") == 0)
	    {
	      func = 20; // xg
	    }
	  else if (strcmp (argv [i], "yg") == 0)
	    {
	      func = 21; // yg
	    }
	  else if (strcmp (argv [i], "yh") == 0)
	    {
	      func = 22; // yh
	    }
	  else if (strcmp (argv [i], "xm") == 0)
	    {
	      func = 23; // xm
	    }
	  else if (strcmp (argv [i], "ym") == 0)
	    {
	      func = 24; // ym
	    }
	  else if (strcmp (argv [i], "zm") == 0)
	    {
	      func = 25; // zm
	    }
	}
      else
	{
	  fprintf (stderr, "$Id: two-body.c,v 5.6 2007/04/12 20:19:59 ichiki Exp $\n");
	  fprintf (stderr, "USAGE\n");
	  fprintf (stderr, "%s [OPTIONS]\n", argv [0]);
	  fprintf (stderr, "\t-h or --help : show this message.\n");
          fprintf (stderr, "\t-o or --out  : output filename"
		   " (default: stdout)\n");
          fprintf (stderr, "\t-n or --nmax : max order (default: 11)\n");
          fprintf (stderr, "\t-n0 or --n0  : starting order (default: 0)\n");
          fprintf (stderr, "\t-c or --cache: set to use memory cache\n");
          fprintf (stderr, "\t-m or --mode : output format\n");
	  fprintf (stderr, "\t\t0 plain txt for f_k with lambda=1"
		   " (default)\n");
          fprintf (stderr, "\t\t1 plain txt each coef of l^q\n");
          fprintf (stderr, "\t\t2 C source\n");
          fprintf (stderr, "\t\t3 latex source for f_k with lambda=1\n");
          fprintf (stderr, "\t\t4 latex source for each coef of l^q\n");
          fprintf (stderr, "\t-f or --func : function name to calculate;\n");
	  fprintf (stderr, "\t\tXA, YA, YB, XC, YC,\n");
	  fprintf (stderr, "\t\tXG, YG, YH, XM, YM, ZM,\n");
	  fprintf (stderr, "\t\tXP, XQ, TQ,\n");
	  fprintf (stderr, "\t\txa, ya, yb, xc, yc.\n");
	  fprintf (stderr, "\t\txg, yg, yh, xm, ym, zm.\n");
	  exit (1);
	}
    }
  
  if (flag_out == 0)
    {
      out = stdout;
    }
  else
    {
      out = fopen (outfile, "w");
      if (out == NULL)
	{
	  fprintf (stderr, "cannot open output file %s.\n", outfile);
	  exit (1);
	}
    }

  if (flag_cache != 0)
    {
      cache_p = table_cache_init ();
      cache_v = table_cache_init ();
      cache_q = table_cache_init ();
      cache_max = 100000;
      cache_i = 0;
    }
  else
    {
      cache_p = NULL;
      cache_v = NULL;
      cache_q = NULL;
    }


  switch (func)
    {
    case 1:
      XA (out, n0, nmax, mode);
      break;
    case 2:
      YA (out, n0, nmax, mode);
      break;
    case 3:
      YB (out, n0, nmax, mode);
      break;
    case 4:
      XC (out, n0, nmax, mode);
      break;
    case 5:
      YC (out, n0, nmax, mode);
      break;
    case 6:
      XG (out, n0, nmax, mode);
      break;
    case 7:
      YG (out, n0, nmax, mode);
      break;
    case 8:
      YH (out, n0, nmax, mode);
      break;
    case 9:
      XM (out, n0, nmax, mode);
      break;
    case 10:
      YM (out, n0, nmax, mode);
      break;
    case 11:
      ZM (out, n0, nmax, mode);
      break;
    case 12:
      XP (out, n0, nmax, mode);
      break;
    case 13:
      XQ (out, n0, nmax, mode);
      break;
    case 14:
      TQ (out, n0, nmax, mode);
      break;
    case 15:
      xa (out, n0, nmax, mode);
      break;
    case 16:
      ya (out, n0, nmax, mode);
      break;
    case 17:
      yb (out, n0, nmax, mode);
      break;
    case 18:
      xc (out, n0, nmax, mode);
      break;
    case 19:
      yc (out, n0, nmax, mode);
      break;
    case 20:
      xg (out, n0, nmax, mode);
      break;
    case 21:
      yg (out, n0, nmax, mode);
      break;
    case 22:
      yh (out, n0, nmax, mode);
      break;
    case 23:
      xm (out, n0, nmax, mode);
      break;
    case 24:
      ym (out, n0, nmax, mode);
      break;
    case 25:
      zm (out, n0, nmax, mode);
      break;
    default:
      fprintf (stderr, "invalid function\n");
      exit (1);
      break;
    }


  if (flag_cache != 0)
    {
      table_cache_free (cache_p);
      table_cache_free (cache_v);
      table_cache_free (cache_q);
    }

  return 0;
}


/** I/O utility routines **/

/* for C source codes */
void
print_mpq_double (FILE *out, mpq_t x)
{
  double d;

  d = mpq_get_d (x);
  fprintf (out, "%.15e", d);
}

void
print_C_header (FILE *out, const char *label)
{
  fprintf (out, "void twobody_%s (int n, double l, double * f)\n", label);
  fprintf (out, "{\n");
}

/*
 * INPUT
 *  q :
 *  q0: order of q in the last term (0 for the first term)
 */
void
print_C_coef_q (FILE *out, mpq_t coef, int q, int q0)
{
  fprintf (out, "    + ");

  int i;
  for (i = q0; i < q; i ++)
    {
      fprintf (out, "l * ");
    }

  fprintf (out, "(");

  print_mpq_double (out, coef);

  fprintf (out, " // q = %d\n", q);
}

void
print_C_close_q (FILE *out, int nq, int k)
{
  fprintf (out, "    ");
  int i;
  for (i = 0; i < nq; i ++)
    {
      fprintf (out, ")");
    }
  if (nq == 0)
    {
      fprintf (out, "0.0;\n");
    }
  else
    {
      fprintf (out, ";\n");
    }
  fprintf (out, "  if (n == %d) return;\n", k);
}

void
print_C_footer (FILE *out)
{
  fprintf (out, "}\n\n");
}


/* for text */
void
print_mpq (FILE *out, mpq_t x)
{
  mpz_out_str (out, 10, mpq_numref (x));

  if (mpz_cmp_ui (mpq_denref (x), 1))
    {
      fprintf (out, "/");
      mpz_out_str (out, 10, mpq_denref (x));
    }
}


void
print_text_coef_q (FILE *out, mpq_t x, int q, int nq)
{
  if (nq > 0 && mpz_sgn (mpq_numref (x)) == 1)
    {
      fprintf (out, " +");
    }
  else
    {
      fprintf (out, " ");
    }
  print_mpq (out, x);

  if (q == 1)
    {
      fprintf (out, " l");
    }
  else
    {
      fprintf (out, " l^%d", q);
    }
}

void
print_text_fk (FILE *out, mpq_t f, int k, const char *label)
{
  fprintf (out, "%sf [%d] = ", label, k);
  print_mpq (out, f);
  fprintf (out, "\n");
}


/* for latex source */
void
print_mpz_formed (FILE *out, mpz_t x)
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

  fprintf (out, "%s", buf2);
}

/*
 * INPUT
 *  flag_plus : 0, don't print "+" for positive x.
 *              1, print "+" for positive x.
 */
void
print_mpq_formed (FILE *out, mpq_t x, int flag_plus)
{
  mpz_t tmp;
  mpz_init (tmp);

  if (flag_plus == 1
      && mpz_sgn (mpq_numref (x)) == 1)
    {
      fprintf (out, "+");
    }
  mpz_set (tmp, mpq_numref (x));
  if (mpz_sgn (mpq_numref (x)) == -1) // negative
    {
      mpz_neg (tmp, tmp);
      fprintf (out, "-");
    }
  if (mpz_cmp_ui (mpq_denref (x), 1))
    {
      /* fraction */
      fprintf (out, "\\frac{");
      //print_mpz_formed (mpq_numref (x));
      print_mpz_formed (out, tmp);
      fprintf (out, "}{");
      print_mpz_formed (out, mpq_denref (x));
      fprintf (out, "}");
    }
  else
    {
      //print_mpz_formed (mpq_numref (x));
      print_mpz_formed (out, tmp);
    }
  mpz_clear (tmp);
}

void
print_latex_coef_q (FILE *out, mpq_t x, int q, int nq)
{
  if (nq != 0 && nq % 4 == 0)
    {
      fprintf (out, "  \\nonumber\\\\\n  &&\n");
    }
  fprintf (out, "  ");
  if (nq == 0) print_mpq_formed (out, x, 0); // no "+"
  else         print_mpq_formed (out, x, 1); // with "+"

  if (q == 0)
    {
      fprintf (out, "\n");
    }
  else if (q == 1)
    {
      fprintf (out, " \\lambda\n");
    }
  else
    {
      fprintf (out, " \\lambda^{%d}\n", q);
    }
}

void
print_latex_fk (FILE *out, mpq_t f, int k, const char *label)
{
  fprintf (out, "\\begin{equation}\n  f^{\\rm %s}_{%d} = ", label, k);
  print_mpq_formed (out, f, 0);
  fprintf (out, "\n\\end{equation}\n");
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
print_fk_header (FILE *out, int flag, const char *label)
{
  if (flag == 2) // C source
    {
      print_C_header (out, label);
    }
}

/* 2 : some header statement for k */
void
print_fk_header_k (FILE *out, int flag, int k, const char *label)
{
  if (flag == 1) // text for each coef of l^q
    {
      fprintf (out, "%sf [%d] = ",
	       label, k);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      fprintf (out, "\\begin{eqnarray}\n  f^{\\rm %s}_{%d} &=&\n",
	       label, k);
    }
  else if (flag == 2) // C source
    {
      fprintf (out, "  f [%d] =\n", k);
    }
}

/* 3 : output the coeficieint for q */
void
print_fk_q (FILE *out, int flag, mpq_t coef, mpq_t f, int q, int q0, int nq)
{
  if (flag == 1) // text for each coef of l^q
    {
      print_text_coef_q (out, coef, q, nq);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      print_latex_coef_q (out, coef, q, nq);
    }
  else if (flag == 2) // C source
    {
      print_C_coef_q (out, coef, q, q0);
    }
}

/* 4 : some footer statement at the end of k */
void
print_fk_footer_k (FILE *out,
		   int flag, int k, mpq_t f, int nq, const char *label)
{
  if (flag == 0) // text for f_k with lambda=1
    {
      print_text_fk (out, f, k, label);
    }
  else if (flag == 1) // text for each coef of l^q
    {
      fprintf (out, "\n");
    }
  else if (flag == 3) // latex for f_k with lambda=1
    {
      print_latex_fk (out, f, k, label);
    }
  else if (flag == 4) // latex for each coef of l^q
    {
      if (nq == 0)
	{
	  if (nq == 0)
	    {
	      fprintf (out, "  0\n");
	    }
	}
      fprintf (out, "\\end{eqnarray}\n");
    }
  else if (flag == 2) // C source
    {
      print_C_close_q (out, nq, k);
    }
  fflush (out);
}

/* 5 : overall footer for the function */
void
print_fk_footer (FILE *out, int flag)
{
  if (flag == 2) // C source
    {
      print_C_footer (out);
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
XA (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "XA");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XA");
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
	      print_fk_q (out, flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XA");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YA (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "YA");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YA");
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
	      print_fk_q (out, flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YA");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YB (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "YB");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  // set two_k with n0
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k = 0, mpq_set_ui (two_k, 2, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YB");
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
	      print_fk_q (out, flag, coef_q, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YB");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XC (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "XC");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XC");
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
		  print_fk_q (out, flag, coef_q, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      else
		{
		  print_fk_q (out, flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XC");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YC (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "YC");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YC");
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
		  print_fk_q (out, flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_q, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YC");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XG (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "XG");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three4);
  mpq_set_ui (three4, 3, 4);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XG");
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
	      print_fk_q (out, flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XG");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

void
YG (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "YG");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three4);
  mpq_set_ui (three4, 3, 4);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YG");
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
	      print_fk_q (out, flag, coef_p, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YG");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three4);
}

void
YH (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "YH");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (mthree8);
  mpq_set_si (mthree8, -3, 8);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YH");
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
		  print_fk_q (out, flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YH");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (mthree8);
}

void
XM (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "XM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XM");
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
		  print_fk_q (out, flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
YM (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "YM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "YM");
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
		  print_fk_q (out, flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "YM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ZM (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_p;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "ZM");

  mpq_init (f);
  mpq_init (coef_p);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "ZM");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  Z_p (2, k - q, q, 10,/* ZM (0, 2, 2)*/ coef_p);

	  if (mpz_cmp_si (mpq_numref (coef_p), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_p, coef_p, two_k);
	      mpq_add (f, f, coef_p);

	      /* 3 : output the coeficieint for q */
	      if (k%2 == 0) // k == even
		{
		  print_fk_q (out, flag, coef_p, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_p, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "ZM");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_p);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
XP (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "XP");

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XP");
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
	      print_fk_q (out, flag, a, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XP");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

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
XQ (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "XQ");

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "XQ");
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
		  print_fk_q (out, flag, a, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, a, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "XQ");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
TQ (FILE *out, int n0, int nmax, int flag)
{
  int k, q, n;
  mpq_t f;
  mpq_t a, b;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "TQ");

  mpq_init (f);
  mpq_init (a);
  mpq_init (b);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "TQ");
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
		  print_fk_q (out, flag, a, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, a, f, q+1, q0, nq);
		  q0 = q + 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "TQ");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (two);
  mpq_clear (two_k);
}


/** mobility functions **/
void
xa (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_u;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "xa");

  mpq_init (f);
  mpq_init (coef_u);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "xa");
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
		  print_fk_q (out, flag, coef_u, f, q, q0, nq);
		  q0 = q;
		}
	      else // k == odd
		{
		  print_fk_q (out, flag, coef_u, f, q-1, q0, nq);
		  q0 = q - 1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "xa");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ya (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "ya");

  mpq_init (f);
  mpq_init (coef_u);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (minus1);

  mpq_set_ui (two, 2, 1);
  mpq_set_si (minus1, -1, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "ya");
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
		  print_fk_q (out, flag, coef_u, f, q, q0, nq);
		  q0 = q;
		}
	      else
		{
		  // is this right?
		  print_fk_q (out, flag, coef_u, f, q-1, q0, nq);
		  q0 = q-1;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "ya");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_u);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (minus1);
}

void
yb (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "yb");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (three);
  mpq_init (minus1);

  mpq_set_ui (two, 2, 1);
  mpq_set_ui (three, 3, 1); // for adjusting the result...
  mpq_set_si (minus1, -1, 1);

  // set two_k with n0
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k = 0, mpq_set_ui (two_k, 2, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "yb");
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
		  print_fk_q (out, flag, coef_q, f, q-1, q0, nq);
		  q0 = q-1;
		}
	      else
		{
		  print_fk_q (out, flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "yb");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (three);
  mpq_clear (minus1);
}

void
xc (FILE *out, int n0, int nmax, int flag)
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
  print_fk_header (out, flag, "xc");

  mpq_init (f);
  mpq_init (coef_q);
  //mpq_init (two);
  //mpq_init (two_k);

  mpq_init (minus1);
  mpq_set_si (minus1, -1, 1);

  /*
  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k = 0, mpq_set_ui (two_k, 2, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));
  */
  for (k = n0;
       k <= nmax;
       k++/*, mpq_mul (two_k, two_k, two)*/)
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "xc");
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
		  print_fk_q (out, flag, coef_q, f, q-2, q0, nq);
		  q0 = q-2;
		}
	      else
		{
		  print_fk_q (out, flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "xc");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  //mpq_clear (two);
  //mpq_clear (two_k);

  mpq_clear (minus1);
}

void
yc (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_q;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "yc");

  mpq_init (f);
  mpq_init (coef_q);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);
  //mpq_set_ui (three, 3, 1); // for adjusting the result...

  // set two_k with n0
  // note that 'two_k' is not 2^k but 2x2^k = 2^(k+1)
  for (k = 0, mpq_set_ui (two_k, 2, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "yc");
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
		  print_fk_q (out, flag, coef_q, f, q-2, q0, nq);
		  q0 = q-2;
		}
	      else
		{
		  print_fk_q (out, flag, coef_q, f, q, q0, nq);
		  q0 = q;
		}
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "yc");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_q);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
xg (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  mpq_t factor;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "xg");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (factor);
  mpq_set_si (factor, -3, 10);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "xg");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  x_E (k - q, q, 1,/* xg */ coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);

	      // for adjusting the factor
	      mpq_mul (coef_E, coef_E, factor);

	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "xg");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (factor);
}

void
yg (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  mpq_t factor;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "yg");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (factor);
  mpq_set_si (factor, -3, 10);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "yg");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  y_E (k - q, q, 2/* for yg */, coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);

	      // for adjusting the factor
	      mpq_mul (coef_E, coef_E, factor);

	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "yg");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (factor);
}

void
yh (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  mpq_t factor;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "yh");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_init (factor);
  mpq_set_si (factor, -9, 20);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "yh");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  y_E (k - q, q, 3/* for yh */, coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);

	      // for adjusting the factor
	      mpq_mul (coef_E, coef_E, factor);

	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "yh");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);

  mpq_clear (factor);
}

void
xm (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "xm");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k < n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "xm");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  x_E (k - q, q, 2,/* xm */ coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);
	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "xm");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
ym (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "ym");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k <= n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "ym");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  y_E (k - q, q, 4/* for ym */, coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);
	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "ym");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);
}

void
zm (FILE *out, int n0, int nmax, int flag)
{
  int k, q;
  mpq_t f;
  mpq_t coef_E;
  mpq_t two;
  mpq_t two_k;

  int nq = 0;
  int q0 = 0;

  /* 1 : output overall header for the function */
  print_fk_header (out, flag, "zm");

  mpq_init (f);
  mpq_init (coef_E);
  mpq_init (two);
  mpq_init (two_k);

  mpq_set_ui (two, 2, 1);

  // set two_k with n0
  for (k = 0, mpq_set_ui (two_k, 1, 1);
       k <= n0;
       k++, mpq_mul (two_k, two_k, two));

  for (k = n0;
       k <= nmax;
       k++, mpq_mul (two_k, two_k, two))
    {
      /* 2 : some header statement for k */
      print_fk_header_k (out, flag, k, "zm");
      mpq_set_ui (f, 0, 1);

      nq = 0;
      q0 = 0;
      for (q = 0; q <= k; q ++)
	{
	  z_E (k - q, q, 0/* for zm */, coef_E); // Epq

	  if (mpz_cmp_si (mpq_numref (coef_E), 0))
	    {
	      /* mul 2^k */
	      mpq_mul (coef_E, coef_E, two_k);
	      mpq_add (f, f, coef_E);

	      /* 3 : output the coeficieint for q */
	      print_fk_q (out, flag, coef_E, f, q, q0, nq);
	      q0 = q;
	      nq ++;
	    }
	}
      /* 4 : some footer statement at the end of k */
      print_fk_footer_k (out, flag, k, f, nq, "zm");
    }

  /* 5 : overall footer for the function */
  print_fk_footer (out, flag);

  mpq_clear (f);
  mpq_clear (coef_E);
  mpq_clear (two);
  mpq_clear (two_k);
}



/** recurrence relation solvers **/
/* calc coef Pnpq for X[AGM,PQ] and x[agm] by Eqs.(JO-3.9) == (JO-8.17)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility function xa
 *             where P1pq = delta_0p delta_0q
 *   mode == 1 for mobility function xg
 *             where P1pq = delta_0p delta_0q (F=1)
 *             where P2pq = 0                 (S=0)
 *   mode == 2 for mobility function xm
 *             where P1pq = 0                 (F=0)
 *             where P2pq = delta_0p delta_0q (S=1)
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
  extern struct table_cache *cache_p;
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
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 1: // mobility for xg

      /* initial condition */
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}
      if (n == 2) return;

      /* search buffer */
      sprintf (filename, "two-body.x_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 2: // mobility for xm

      /* initial condition */
      if (n == 1) return;
      if (n == 2)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.x_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
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
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
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

  append_result (filename, cache_p, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}
/* calc coef Vnpq for X[AGM, PQ] and x[agm] by Eqs.(JO-3.8) == (JO-8.15)
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility function xa
 *             where P1pq = delta_0p delta_0q
 *   mode == 1 for mobility function xg
 *             where P1pq = delta_0p delta_0q (F=1)
 *             where P2pq = 0                 (S=0)
 *   mode == 2 for mobility function xm
 *             where P1pq = 0                 (F=0)
 *             where P2pq = delta_0p delta_0q (S=1)
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
  extern struct table_cache *cache_v;
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
    case 1: // mobility for xg
    case 2: // mobility for xm

      /* Vnpq for mobility has no initial condition */

      /* search buffer */
      sprintf (filename, "two-body.x_v.%d", mode);
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
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
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
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

  append_result (filename, cache_v, n, p, q, coef_v);

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
  extern struct table_cache *cache_q;
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
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
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
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
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

  append_result (filename, cache_q, n, p, q, coef_q);

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
 *   mode == 2 for mobility functions yg
 *             where P1pq = delta_0p delta_0q  (F=1)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0                  (S=0)
 *   mode == 3 for mobility functions yh
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = delta_0p delta_0q  (T=1)
 *             and   P2pq = 0                  (S=0)
 *   mode == 4 for mobility functions ym
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
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
  extern struct table_cache *cache_p;
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
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 1: // yc : P1pq = 0 and Q1pq = delta_0p delta_0q

      /* initial condition */
      if (n == 1)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 2: // yg : P1pq = delta_0p delta_0q, Q1pq = P2pq = 0

      /* initial condition */
      if (n == 1)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}
      if (n == 2) return;

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 3: // yh : P1pq = P2pq = 0, Q1pq = delta_0p delta_0q

      /* initial condition */
      if (n == 1 || n == 2)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 4: // ym : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* initial condition */
      if (n == 1) return;
      if (n == 2)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
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
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
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

  append_result (filename, cache_p, n, p, q, coef_p);

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
 *   mode == 2 for mobility functions yg
 *             where P1pq = delta_0p delta_0q  (F=1)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0                  (S=0)
 *   mode == 3 for mobility functions yh
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = delta_0p delta_0q  (T=1)
 *             and   P2pq = 0                  (S=0)
 *   mode == 4 for mobility functions ym
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
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
  extern struct table_cache *cache_v;
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
    case 1: // yc : P1pq = 0 and Q1pq = delta_0p delta_0q
    case 2: // yg : P1pq = delta_0p delta_0q, Q1pq = P2pq = 0
    case 3: // yh : P1pq = P2pq = 0, Q1pq = delta_0p delta_0q
    case 4: // ym : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* Vnpq for mobility has no initial condition */

      /* search buffer */
      sprintf (filename, "two-body.y_v.%d", mode);
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
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
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
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

  append_result (filename, cache_v, n, p, q, coef_v);

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
 *   mode == 2 for mobility functions yg
 *             where P1pq = delta_0p delta_0q  (F=1)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0                  (S=0)
 *   mode == 3 for mobility functions yh
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = delta_0p delta_0q  (T=1)
 *             and   P2pq = 0                  (S=0)
 *   mode == 4 for mobility functions ym
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
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
  extern struct table_cache *cache_q;
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
    case 2: // yg : P1pq = delta_0p delta_0q, Q1pq = P2pq = 0
    case 4: // ym : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* initial condition */
      if (n == 1)
	{
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.y_q.%d", mode);
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    case 1: // yc : P1pq = 0 and Q1pq = delta_0p delta_0q
    case 3: // yh : P1pq = P2pq = 0, Q1pq = delta_0p delta_0q

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
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
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
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
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

  append_result (filename, cache_q, n, p, q, coef_q);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Pnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions zm
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
 *   mode == 10 for resistance function ZM
 *             where Qn00 = 0
 *                   Pn00 = delta_2n
 *                   Vn00 = delta_2n
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
Z_p (int n, int p, int q, int mode, mpq_t coef_p)
{
  extern struct table_cache *cache_p;
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0, tmp1;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (n < 0
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // zm : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* initial condition */
      if (n == 1) return;
      if (n == 2)
	{
	  if (p == 0 && q == 0)
	    {
	      mpq_set_ui (coef_p, 1, 1);
	    }
	  return;
	}

      /* search buffer */
      sprintf (filename, "two-body.z_p.%d", mode);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */

      break;

    case 10: // resistance function ZM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  // for all resistance functions
  if (mode == 10)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inp != 0 && n == inp)
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
      sprintf (filename, "two-body.Z_p.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      //Z_p (s, q - s, p - n + 1, inq, inp, inv, b);
      Z_p (s, q - s, p - n + 1, mode, b);
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
      //Z_p (s, q - s, p - n - 1, inq, inp, inv, b);
      Z_p (s, q - s, p - n - 1, mode, b);
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
      //Z_v (s, q - s - 2, p - n + 1, inq, inp, inv, b);
      Z_v (s, q - s - 2, p - n + 1, mode, b);
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
      //Z_q (s, q - s - 1, p - n + 1, inq, inp, inv, b);
      Z_q (s, q - s - 1, p - n + 1, mode, b);
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

  append_result (filename, cache_p, n, p, q, coef_p);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
}

/* calc coef Vnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions zm
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
 *   mode == 10 for resistance function ZM
 *             where Qn00 = 0
 *                   Pn00 = delta_2n
 *                   Vn00 = delta_2n
 * OUTPUT
 *   coef_v : (mpq_t)
 */
void
Z_v (int n, int p, int q, int mode, mpq_t coef_v)
{
  extern struct table_cache *cache_v;
  int s;

  mpq_t a;
  mpq_t tmp;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Vnpq */
  mpq_set_ui (coef_v, 0, 1);

  if (n < 0
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // zm : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* Vnpq for mobility has no initial condition */

      /* search buffer */
      sprintf (filename, "two-body.z_v.%d", mode);
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
	return; /* found */

      break;

    case 10: // resistance function ZM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  // for all resistance functions
  if (mode == 10)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inv != 0 && n == inv)
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
      sprintf (filename, "two-body.Z_v.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
	return; /* found */
    }

  mpq_init (a);
  mpq_init (tmp);

  /* Vnpq : 1st term Pnpq */
  //Z_p (n, p, q, inq, inp, inv, coef_v);
  Z_p (n, p, q, mode, coef_v);

  for (s = 1; s <= q; s ++)
    {
      /* Vnpq : 2nd term Ps(q-s)(p-n-1) */
      //Z_p (s, q - s, p - n - 1, inq, inp, inv, a);
      Z_p (s, q - s, p - n - 1, mode, a);
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

  append_result (filename, cache_v, n, p, q, coef_v);

  mpq_clear (a);
  mpq_clear (tmp);
}

/* calc coef Qnpq on ZM
 * INPUT
 *   n, p, q : (int)
 *   mode == 0 for mobility functions zm
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
 *   mode == 10 for resistance function ZM
 *             where Qn00 = 0
 *                   Pn00 = delta_2n
 *                   Vn00 = delta_2n
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
Z_q (int n, int p, int q, int mode, mpq_t coef_q)
{
  extern struct table_cache *cache_q;
  int s;

  mpq_t a, b;
  mpq_t tmp;

  char filename [40];
  int inq = 0;
  int inp = 0;
  int inv = 0;


  /* clear Vnpq */
  mpq_set_ui (coef_q, 0, 1);

  if (n < 0
      || p < 0
      || q < 0)
    {
      return;
    }

  switch (mode)
    {
    case 0: // zm : P1pq = 0, P2pq = delta_0p delta_0q, Q1pq = 0

      /* initial condition */
      if (n == 1) return;

      /* search buffer */
      sprintf (filename, "two-body.z_q.%d", mode);
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
	return; /* found */

      break;

    case 10: // resistance function ZM
      inq = 0;
      inp = 2;
      inv = 2;

      break;

    default:
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }

  // for all resistance functions
  if (mode == 10)
    {
      if (p == 0 && q == 0)
	{
	  /* initial condition */
	  if (inq != 0 && n == inq)
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
      sprintf (filename, "two-body.Z_q.%d.%d.%d", inq, inp, inv);
      if (search_results (filename, cache_q, n, p, q, coef_q) == 0)
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
      //Z_q (s, q - s - 1, p - n, inq, inp, inv, b);
      Z_q (s, q - s - 1, p - n, mode, b);

      if (mpz_cmp_si (mpq_numref (b), 0))
	{
	  /* mul s */
	  mpq_set_si (tmp, s, 1);
	  mpq_mul (b, b, tmp);

	  mpq_add (a, a, b);
	}
      /* Qnpq : 2nd term Ps(q-s)(p-n) */
      //Z_p (s, q - s, p - n, inq, inp, inv, b);
      Z_p (s, q - s, p - n, mode, b);

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

  append_result (filename, cache_q, n, p, q, coef_q);

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
  extern struct table_cache *cache_p;
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
  if (search_results (filename, cache_p, n, p, q, coef_p) == 0)
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

  append_result (filename, cache_p, n, p, q, coef_p);

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
  extern struct table_cache *cache_v;
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
  if (search_results (filename, cache_v, n, p, q, coef_v) == 0)
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

  append_result (filename, cache_v, n, p, q, coef_v);

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


/* calc coef Epq for xg (basically -Pnpq with n=2)
 * INPUT
 *   p, q : (int)
 *   mode : pass to X_[pv](). only the following values are valid;
 *   mode == 1 for mobility function xg
 *             where P1pq = delta_0p delta_0q (F=1)
 *             where P2pq = 0                 (S=0)
 *   mode == 2 for mobility function xm
 *             where P1pq = 0                 (F=0)
 *             where P2pq = delta_0p delta_0q (S=1)
 * OUTPUT
 *   coef_q : (mpq_t)
 */
void
x_E (int p, int q, int mode, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp;

  int n;

  if (mode != 1 && mode != 2)
    {
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }
  n = 2;

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  if (p == 0 && q == 0)
    {
      // from B.C. (P2pq = delta_0p delta_0q)
      if (mode == 2)
	{
	  mpq_set_ui (coef_p, 1, 1);
	}
      return;
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

  // coef_p = - coef_p
  mpq_set_si (tmp, -1, 1);
  mpq_mul (coef_p, coef_p, tmp);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
}

/* calc coef Epq for yg, yh (basically -Pnpq with n=2)
 * INPUT
 *   p, q : (int)
 *   mode : pass to Y_[pvq](). only the following values are valid;
 *   mode == 2 for mobility functions yg
 *             where P1pq = delta_0p delta_0q  (F=1)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0                  (S=0)
 *   mode == 3 for mobility functions yh
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = delta_0p delta_0q  (T=1)
 *             and   P2pq = 0                  (S=0)
 *   mode == 4 for mobility functions ym
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
y_E (int p, int q, int mode, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0;

  int n;

  if (mode != 2 &&
      mode != 3 &&
      mode != 4)
    {
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }
  n = 2;

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  if (p == 0 && q == 0)
    {
      // from B.C. (P2pq = delta_0p delta_0q)
      if (mode == 4)
	{
	  mpq_set_ui (coef_p, 1, 1);
	}
      return;
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

  // coef_p = - coef_p
  mpq_set_si (tmp, -1, 1);
  mpq_mul (coef_p, coef_p, tmp);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
}

/* calc coef Epq for zm (basically -Pnpq with n=2 by Z_p())
 * INPUT
 *   n, p, q : (int)
 *   mode : pass to Z_[pvq](). only the following values are valid;
 *   mode == 0 for mobility functions zm
 *             where P1pq = 0                  (F=0)
 *             and   Q1pq = 0                  (T=0)
 *             and   P2pq = 0delta_0p delta_0q (S=1)
 * OUTPUT
 *   coef_p : (mpq_t)
 */
void
z_E (int p, int q, int mode, mpq_t coef_p)
{
  int s;

  mpq_t a, b;
  mpq_t tmp, tmp0, tmp1;

  int n;

  if (mode != 0)
    {
      fprintf (stderr, "invalid mode\n");
      exit (1);
    }
  n = 2;

  /* clear Pnpq */
  mpq_set_ui (coef_p, 0, 1);

  if (p < 0
      || q < 0)
    {
      return;
    }

  if (p == 0 && q == 0)
    {
      // from B.C. (P2pq = delta_0p delta_0q)
      if (mode == 0)
	{
	  mpq_set_ui (coef_p, 1, 1);
	}
      return;
    }

  mpq_init (a);
  mpq_init (b);
  mpq_init (tmp);
  mpq_init (tmp0);
  mpq_init (tmp1);

  for (s = 1; s <= q; s ++)
    {
      /* Pnpq */
      mpq_set_ui (a, 0, 1);
      /* Pnpq : 1st term Ps(q-s)(p-n+1) */
      Z_p (s, q - s, p - n + 1, mode, b);
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
      Z_p (s, q - s, p - n - 1, mode, b);
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
      Z_v (s, q - s - 2, p - n + 1, mode, b);
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
      Z_q (s, q - s - 1, p - n + 1, mode, b);
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

  // coef_p = - coef_p
  mpq_set_si (tmp, -1, 1);
  mpq_mul (coef_p, coef_p, tmp);

  mpq_clear (a);
  mpq_clear (b);
  mpq_clear (tmp);
  mpq_clear (tmp0);
  mpq_clear (tmp1);
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


/** table handling routines **/
/* table cache */
struct table_cache *
table_cache_init (void)
{
  struct table_cache *c = NULL;

  c = (struct table_cache *) malloc (sizeof (struct table_cache));
  if (c == NULL)
    {
      fprintf (stderr, "allocation error in table_cache_init()\n");
      exit (1);
    }

  c->number = 0;
  c->n = NULL;
  c->p = NULL;
  c->q = NULL;
  c->coef = NULL;

  return (c);
}

void
table_cache_free (struct table_cache *c)
{
  if (c != NULL)
    {
      if (c->n != NULL) free (c->n);
      if (c->p != NULL) free (c->p);
      if (c->q != NULL) free (c->q);
      if (c->coef != NULL)
	{
	  int i;
	  for (i = 0; i < c->number; i ++)
	    {
	      mpq_clear (c->coef [i]);
	    }
	  free (c->coef);
	}
      free (c);
    }
}

void
table_cache_append (struct table_cache *c,
		    int n, int p, int q, mpq_t coef)
{
  c->number ++;
  c->n = (int *) realloc (c->n, sizeof (int) * (c->number));
  c->p = (int *) realloc (c->p, sizeof (int) * (c->number));
  c->q = (int *) realloc (c->q, sizeof (int) * (c->number));
  c->coef = (mpq_t *) realloc (c->coef, sizeof (mpq_t) * (c->number));

  int i = c->number - 1; // last element
  c->n [i] = n;
  c->p [i] = p;
  c->q [i] = q;
  mpq_init (c->coef [i]);
  mpq_set (c->coef [i], coef);
}

void
table_cache_set_i (struct table_cache *c,
		   int i,
		   int n, int p, int q, mpq_t coef)
{
  if (i >= c->number)
    {
      fprintf (stderr, "out of range\n");
      exit (1);
    }

  c->n [i] = n;
  c->p [i] = p;
  c->q [i] = q;
  mpq_set (c->coef [i], coef);
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
search_results (char *file, struct table_cache *cache,
		int n, int p, int q, mpq_t coef)
{
  extern int cache_max; // max size of cache
  extern int cache_i; // current index

  FILE *res;
  mpz_t num, den;
  int nn, pp, qq;
  char string [1024];
  int i;


  /* first search memory of table_cache */
  if (cache != NULL)
    {
      for (i = 0; i < cache->number; i ++)
	{
	  if (cache->n [i] == n &&
	      cache->p [i] == p &&
	      cache->q [i] == q)
	    {
	      /* found */
	      mpq_set (coef, cache->coef [i]);
	      return (0);
	    }
	}
    }

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

	      // asign table_cache
	      if (cache != NULL)
		{
		  if (cache->number < cache_max)
		    {
		      table_cache_append (cache, n, p, q, coef);
		    }
		  else
		    {
		      table_cache_set_i (cache, cache_i, n, p, q, coef);
		      cache_i ++;
		      if (cache_i >= cache_max)
			{
			  cache_i -= cache_max;
			}
		    }
		}

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
append_result (char *file, struct table_cache *cache,
	       int n, int p, int q, mpq_t coef)
{
  extern int cache_max; // max size of cache
  extern int cache_i; // current index

  FILE *res;

  // asign table_cache
  if (cache != NULL)
    {
      if (cache->number < cache_max)
	{
	  table_cache_append (cache, n, p, q, coef);
	}
      else
	{
	  table_cache_set_i (cache, cache_i, n, p, q, coef);
	  cache_i ++;
	  if (cache_i >= cache_max)
	    {
	      cache_i -= cache_max;
	    }
	}
    }

  res = fopen (file, "a");
  if (res != NULL)
    {
      flockfile (res); // lock

      fprintf (res, "%d %d %d ", n, p, q);
      print_mpq (res, coef);
      fprintf (res, "\n");

      funlockfile (res); // unlock

      fclose (res);
    }
  else
    {
      fprintf (stderr, "cannot open %s\n", file);
    }
}
