#define MAIN
#define _INCLUDE_XOPEN_SOURCE
/**************************************************************/
/* Jeffrey & Onish 1986                                       */
/**************************************************************/

#define N 100

#include <math.h>
#include <stdio.h>	/* printf() fprintf() */
#include <stdlib.h> /* exit() */

/* global variable */

/* function prototype */
int main( int argc, char** argv );
double combi(int n,int m);

/* main program */
int main( int argc, char** argv ){
  char
      flg;
  int
      i,
      n,p,q,s,
      nmax;

  double
      Q[N][N][N],
      t,f;

  /* option analysis */
  flg = 0;
  for( i=1;i<argc;i++ ){
    if( strcmp( argv[i],"-n" ) == 0){
      if ( i+1 < argc ) {
	nmax = atoi( argv[++i] );
	flg ++;
      }
    } else {
      fprintf(stderr,"USAGE\n");
      fprintf(stderr,"XC [options]\n");
      fprintf(stderr,"OPTIONS\n");
      fprintf(stderr,"  -n number of coefficients\n");
      exit( 1 );
    }
  }
  if(flg!=1){
      fprintf(stderr,"You Must Believe In Spring.\n");
      exit(1);
  }

  /* zero clear */
  for(n=0;n<=nmax;n++){
      for(p=0;p<=nmax;p++){
	  for(q=0;q<=nmax;q++){
	      Q[n][p][q] = 0.0;
	  }
      }
  }
  /* initial condition */
  Q[1][0][0] = 1.0;

  for(i=1;i<=nmax;i++){
      for(q=1;q<=nmax;q++){
	  p = i-q;
	  for(n=1;n<=nmax;n++){
	      if(p>=n){
		  for(s=1;s<=(q-1);s++){
		      Q[n][p][q] +=
			  combi(n+s,n)*(double)s/(double)(n+1)
			      *Q[s][q-1-s][p-n];
		  }
	      }
	  }
      }
  }

  t = 1.0;
  for(i=0;i<nmax;i++){
      f=0.0;
      for(q=0;q<=i;q++){
	  p = i-q;
	  f += Q[1][p][q];
      }
      fprintf(stdout,"%d %.15e x 2^%d\n",i,f,i);
      t*=2.0;
  }
}

double combi(int n,int m){
    int i,j;
    double x,y,z;

    if(n<0 || m<0 || n<m){
	return(0.0);
    }
    x=1.0;
    for(i=1;i<=n;i++){
	x*=(double)i;
    }
    y=1.0;
    for(i=1;i<=m;i++){
	y*=(double)i;
    }
    z=1.0;
    for(i=1;i<=(n-m);i++){
	z*=(double)i;
    }
    return(x/y/z);
}
