/*
 *     dynamic_arrays.c
 * 	Numerical Recipes utilities (nrutil) array types
 * (reference: http://www.nr.com/pubdom/nrutil.c.txt)
 *  Created on: Oct 9, 2013
 *      Author: Shangwu Xiong
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "user_defined.h"
/********************************************************************/
/*                                                                  */
/*   standard error                                                 */
/*                                                                  */
/********************************************************************/

void nrerror (char error_text[])
{       /* NR standard error handler */
        fprintf(stderr,"Numerical recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}
/*****************************************************************/
/*                                                               */
/*    free space for 1-D, 2-D, 3-D arrays                        */
/*                                                               */
/*****************************************************************/
void free_dvector(double *v, int nl, int nh)
{
        free((FREE_ARG) (v+nl-NR_END));
}
void free_ivector(int *v, int nl, int nh)
{
        free((FREE_ARG) (v+nl-NR_END));
}
void free_dmatrix(double **mm, int nrl, int nrh, int ncl, int nch)
{
        free((FREE_ARG) (mm[nrl]+ncl-NR_END));
        free((FREE_ARG) (mm+nrl-NR_END));
}
void free_d3tensor(double ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}
/********************************************************************/
/*                                                                  */
/*  dynamic allocation of 1-D array, double precision.              */
/*                                                                  */
/********************************************************************/
double *dvector(int nl, int nh)
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) nrerror("allocation failure in dvextor()");
        return v-nl+NR_END;
}

int *ivector(int nl, int nh)
{
        int *v;

        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) nrerror("allocation failure in dvextor()");
        return v-nl+NR_END;
}
/********************************************************************/
/*                                                                  */
/*  dynamic allocation of 2-D array, double precision.              */
/*                                                                  */
/********************************************************************/
double **dmatrix (int nrl, int nrh, int ncl, int nch)
{
        long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **mm;

        /* allocate pointers to rows */
        mm=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!mm) nrerror("allocation failure 1 in matrix()");
        mm += NR_END;
        mm -= nrl;

        /* allocate rows and set pointers to them */
        mm[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!mm[nrl]) nrerror("allocation failure 2 in matrix()");
        mm[nrl] += NR_END;
        mm[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) mm[i]=mm[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return mm;
}
/********************************************************************/
/*                                                                  */
/*  dynamic allocation of 3-D array, double precision.              */
/*                                                                  */
/********************************************************************/

double ***d3tensor (int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
       long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
       double ***t;

       t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
       if (!t) nrerror("allocation failure 1 in d3tensor()");
       t += NR_END;
       t-=nrl;

       t[nrl]=(double**) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
       if (!t[nrl]) nrerror("allocation failure 1 in d3tensor()");
       t[nrl] +=NR_END;
       t[nrl] -=ncl;

       t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
       if (!t[nrl][ncl]) nrerror("allocation failure 1 in d3tensor()");
       t[nrl][ncl] +=NR_END;
       t[nrl][ncl] -=ndl;

       for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
       for (i=nrl+1;i<=nrh;i++) {
            t[i]=t[i-1]+ncol;
            t[i][ncl]=t[i-1][ncl]+ncol*ndep;
            for (j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
       }

       return t;
}

