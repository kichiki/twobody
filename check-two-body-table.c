/* check the result table created by two-body.
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: $
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


/* check the table
 * INTPUT
 *   char *file : file name
 * OUTPUT
 */
void
check_table (char *file)
{
  int check = 0;

  FILE *table = NULL;

  mpq_t coef;
  mpz_t num, den;
  int nn, pp, qq;
  char string [1024];
  int i;


  /* search buffer */
  table = fopen (file, "r");
  if (table == NULL)
    {
      fprintf (stderr, "fail to open file %s\n", file);
      return;
    }

  int l;
  for (l = 1;; l++)
    {
      i = fscanf (table, "%d %d %d %s\n", &nn, &pp, &qq, string);
      if (i == EOF)
	{
	  break;
	}
      else if (i != 4)
	{
	  fprintf (stdout, "%s : (line %d) wrong format : %d %d %d %s\n",
		   file, l, nn, pp, qq, string);
	  check ++;
	  continue;
	}

      mpq_init (coef);
      mpz_init (num);
      mpz_init (den);

      for (i=0;
	   string [i] != '/'
	     && string [i] != '\0';
	   i ++);

      if (string [i] == '/')
	{
	  string [i] = '\0';
	  int status;
	  status = mpz_set_str (num, &string [0], 10);
	  if (status != 0)
	    {
	      fprintf (stdout, "%s : (line %d) fail to parse num :"
		       " %d %d %d %s/%s\n",
		       file, l, nn, pp, qq, string, string + i + 1);
	      check ++;
	      continue;
	    }
	  status = mpz_set_str (den, &string [i + 1], 10);
	  if (status != 0)
	    {
	      fprintf (stdout, "%s : (line %d) fail to parse den :"
		       " %d %d %d %s/%s\n",
		       file, l, nn, pp, qq, string, string + i + 1);
	      check ++;
	      continue;
	    }
	}
      else
	{
	  mpz_set_str (num, &string [0], 10);
	  mpz_set_ui (den, 1);
	}
      mpq_set_num (coef, num);
      mpq_set_den (coef, den);
      mpq_canonicalize (coef);

      if (mpz_cmp (num, mpq_numref (coef)) != 0
	  || mpz_cmp (den, mpq_denref (coef)) != 0)
	{
	  fprintf (stdout, "%s : (line %d) coef is not canonical :"
		   " %d %d %d %s\n",
		   file, l, nn, pp, qq, string);
	  check ++;
	  continue;
	}

      mpz_clear (num);
      mpz_clear (den);
      mpq_clear (coef);
    }
  fclose (table);

  if (check == 0)
    {
      fprintf (stdout, "%s: PASSED\n", file);
    }
}


int
check_uniqueness_sub (FILE *table, fpos_t *top, int l,
		      int n, int p, int q, const char *coef_str)
{
  int check = 0;

  int status;
  status = fsetpos (table, top);
  if (status != 0)
    {
      fprintf (stderr, "check_uniqueness_sub: fail on fsetpos()\n");
      exit (1);
    }

  int nn, pp, qq;
  char string [1024];
  l++; // l is the line for the reference, so top is (l+1)-th line.
  for (;; l++)
    {
      int i;
      i = fscanf (table, "%d %d %d %s\n", &nn, &pp, &qq, string);
      if (i == EOF)
	{
	  break;
	}
      else if (i != 4)
	{
	  fprintf (stdout, "(line %d) wrong format : %d %d %d %s\n",
		   l, nn, pp, qq, string);
	  continue;
	}

      if (n == nn &&
	  p == pp &&
	  q == qq)
	{
	  check ++;
	  if (strcmp (coef_str, string) == 0)
	    {
	      fprintf (stdout, "(line %d) duplicated entry :"
		       " %d %d %d %s\n",
		       l, nn, pp, qq, string);
	    }
	  else
	    {
	      fprintf (stdout, "(line %d) duplicated entry with different coef:"
		       " %d %d %d %s != %s\n",
		       l, nn, pp, qq, string, coef_str);
	    }
	}
    }

  return (check);
}

/* check the table
 * INTPUT
 *   char *file : file name
 * OUTPUT
 */
void
check_uniqueness (char *file)
{
  fprintf (stdout, "start the uniqueness check\n");

  FILE *table = NULL;

  int n, p, q;
  char string [1024];
  int i;


  /* search buffer */
  table = fopen (file, "r");
  if (table == NULL)
    {
      fprintf (stderr, "fail to open file %s\n", file);
      return;
    }

  fpos_t cur;
  int status;
  int l;
  for (l = 1;; l++)
    {
      i = fscanf (table, "%d %d %d %s\n", &n, &p, &q, string);
      if (i == EOF)
	{
	  break;
	}
      else if (i != 4)
	{
	  fprintf (stdout, "(line %d) wrong format : %d %d %d %s\n",
		   l, n, p, q, string);
	  continue;
	}

      fprintf (stdout, "[line %d : %d : %d %d %d]\n", l, (n+p+q), n, p, q);

      status = fgetpos (table, &cur);
      if (status != 0)
	{
	  fprintf (stderr, "check_uniqueness: fail on fgetpos()\n");
	  exit (1);
	}

      check_uniqueness_sub (table, &cur, l,
			    n, p, q, string);

      // reset the current position on FILE *table
      status = fsetpos (table, &cur);
      if (status != 0)
	{
	  fprintf (stderr, "check_uniqueness: fail on fsetpos()\n");
	  exit (1);
	}
    }
}


int
main (int argc, char** argv)
{
  char tablefile [256];

  /* option analysis */
  int i;
  int flag_file = 0;
  int flag_uniq = 0;
  for (i = 1; i < argc; i++)
    {
      if (strcmp (argv [i], "-f") == 0 ||
	  strcmp (argv [i], "--file") == 0)
	{
	  i++;
	  if (strlen (argv [i]) > 256)
	    {
	      fprintf (stderr, "too long file name\n");
	      exit (1);
	    }
	  strcpy (tablefile, argv [i]);
	  flag_file = 1;
	}
      else if (strcmp (argv [i], "-u") == 0 ||
	       strcmp (argv [i], "--uniq") == 0)
	{
	  flag_uniq = 1;
	}
      else
	{
	  fprintf (stderr, "$Id: $\n");
	  fprintf (stderr, "USAGE\n");
	  fprintf (stderr, "%s [OPTIONS]\n", argv [0]);
	  fprintf (stderr, "\t-f or --file : give the table file to check.\n");
	  fprintf (stderr, "\t-u or --uniq : check the uniqueness"
		   " (default: format check only)\n\n");
	  exit (1);
	}
    }

  if (flag_file == 1)
    {
      check_table (tablefile);

      if (flag_uniq == 1)
	{
	  check_uniqueness (tablefile);
	}
    }

  return 0;
}
