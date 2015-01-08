/*
 * Poisson_SPH.c
 *
 *This code was developed on the basis of theoretical works:
 *1)S.Kulasegaram, Development of particle based meshless method with applications in metal forming simulations.
 *               Doctor Thesis, University of Wales, 1999.
 *2)C.S.Li, S.W.Xiong, J.M.C.Rodrigues, P.A.F.Martins. ‘Simulation of Metal Forming Processes by Meshless Methods’.
 *   Northeastern University Press,Shenyang, China, 2004.8. (Chinese version, ISBN 7-81102-067-X).
 *
 *  Created on: Oct 9, 2013
 *  Author: Shangwu Xiong
 *
 */
#include <stdio.h>
#include <math.h>
#include "user_defined.h"
//#include "dynamic_arrays.c" //not needed because functions are declared as external//
	/*** memory allocation ***/

//extern void dynamic_allocation(void);
//extern void free_space(void);
extern double *dvector(int, int);
extern int *ivector(int, int);
extern double **dmatrix(int, int, int, int);
extern double ***d3tensor(int, int, int, int, int, int);
extern void free_dvector(double *, int, int);
extern void free_ivector(int *, int, int);
extern void free_dmatrix(double **, int, int, int, int);
extern void free_d3tensor(double ***, int, int, int, int, int, int);
extern void nrerror(char error_text[]);
/*******************************************************************/

//
 int *nkode,*nkbe,*ngama;
 double **coord,*fi,*vlm,*wxh,*ax, **dax,**beta,**dbetax,**dbetay;
 double *dist,**gama,*axs,**betas,*a2,*ffbe,**dff,*amtx,**gama0,*gamas;
//
//
//
 double *PosX, *PosY, *PosZ, **H,*Temp,*temp_old,**h_old,**pressure3;
 double **Phi_new,**Phi_old,**Hsurface,**F_new,**pressure1,**TempSurf;
 double **Deform_E1,**Deform_E2,**IC,*Posx_E,*Poy_E,*PosZ_E,**press2;
 int **lnods,**nodei,**nodlc,*noden,*nsvab,*nb_index;
 int **nlxi,**nlzi,**id;
 double  **deltd,**resid,*fthick1,*fthick0;
 double *ematrix1,**ematrix2;
//
 extern  void input1(int,int);
 //

 extern void Poisson_SPH_Equation(double **,double **,double **,double **, double**,double**,double **,\
		  double**, int *,	int *,int *,double *,double *,double *,double *,double *,double *,double *,\
				double *,double *,double **,double *);
 extern void Bearing_main1(double** ,double **,double **,double **, double**,double**,double **,\
 		int *,	int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,\
 		double *,double **,double *);

//*The parameters argc, argument count, and argv, argument vector
//int main(int argc, char *argv[])
 int main(void)
{   int i,m,n,ii9;
  setvbuf( stdout, NULL, _IONBF, 0 );
/*     printf("     Welcome to use Poisson_SPH \n\n\n");
     printf("       Start: for fun only (e.g. 1 2 3 4 5)\n");
	 printf( "argc = %d\n", argc );
	    for( i = 0; i < argc; ++i ) {
	        printf( "argv[ %d ] = %s\n", i, argv[ i ] );
	    };
     printf(" *****  End: fun ********\n\n");
*/
     printf( "\n\n\n\n '           WELCOME TO USE MESHLESS3	       '");
     printf(     "\n\n '--------------------------------------------'");
     printf(     "\n\n '               MESH FREE 2D PROGRAM         '");
     printf(     "\n\n '              SPH to solve Poisson equation '");
     printf(     "\n\n '                  2012/Oct/08               '");
     printf(     "\n\n '                                            '");
	 printf(     "\n\n '           Dr.Shangwu Xiong                 '");
     printf(     "\n\n '--------------------------------------------'");

     //c
             ii9=1;
     /*c
     c	  read input data
     c*/
             input1(5,ii9);
     //c
             npoin=npoin_x*npoin_y;
             nwork=(npoin+1)*npoin/2+1;
             npbe=(npoin_x+npoin_y)*2;
     //c
             ratio=(xmax-xmin)*(npoin_y-1.)/((ymax-ymin)*(npoin_x-1.));
             if(ratio<1.)ratio=1./ratio;
		   if(nmethod<=1)
		    {ratio=(sqrt(npoin*1.)-1)*4.*alfa*ratio;}
		   else
 			{ratio=(sqrt(npoin*1.)-1)*3.*alfa*ratio;}          		   
             nnu=npoin*ratio;
             if(nnu>nwork)nnu=nwork;


     		n=npoin_x;
     		m=npoin_y;
     //
         // 2-dimensional arrays ******
     		coord=dmatrix(1,2,1,npoin);
     		dax=dmatrix(1,2,1,npoin);
     		beta=dmatrix(1,2,1,npoin);
     		dbetax=dmatrix(1,2,1,npoin);
     		dbetay=dmatrix(1,2,1,npoin);
     		gama=dmatrix(1,2,1,npoin);
     		gama0=dmatrix(1,2,1,npoin);
     		betas=dmatrix(1,2,1,npoin);
     //************************************
     		nkode=ivector(1,npoin);
     		nkbe=ivector(1,npbe+1);
     		ngama=ivector(1,npoin);
     		fi=dvector(1,npoin);
     		vlm=dvector(1,npoin);
     		wxh=dvector(1,npoin);
     		ax=dvector(1,npoin);
     		axs=dvector(1,npoin);
     		gamas=dvector(1,npoin);
     		dist=dvector(1,npoin);
     		a2=dvector(1,npoin);
     		ffbe=dvector(1,npoin);
//     	    dff=dvector(1,npoin);
     	  if(nmethod>=2){
     		fthick0=dvector(1,npoin);
     		fthick1=dvector(1,npoin);};
     		dff=dmatrix(1,npoin,0,0);
     		amtx=dvector(1,nnu);
//
     	if(nmethod<=1)
     	{	Poisson_SPH_Equation(coord,dax,dbetax,dbetay, gama,gama0,beta,\
     		  betas, nkode, nkbe,ngama,fi,vlm,wxh,ax,gamas,axs,dist,\
     				a2,ffbe,dff,amtx);     		//c
     	}
     	else
     		{Bearing_main1(coord,dax,dbetax,dbetay, gama,gama0,beta, \
     				      nkode, nkbe,ngama,fi,vlm,wxh,ax,fthick1,fthick0,dist,a2,\
    				                      ffbe,dff,amtx);}
//
 		free_dmatrix(coord,1,2,1,npoin);
 		free_dmatrix(dax,1,2,1,npoin);
 		free_dmatrix(beta,1,2,1,npoin);
 		free_dmatrix(dbetax,1,2,1,npoin);
 		free_dmatrix(dbetay,1,2,1,npoin);
 		free_dmatrix(gama,1,2,1,npoin);
 		free_dmatrix(gama0,1,2,1,npoin);
 		free_dmatrix(betas,1,2,1,npoin);

 		free_ivector(nkode,1,npoin);
 		free_ivector(nkbe,1,npbe+1);
 		free_ivector(ngama,1,npoin);
 		free_dvector(fi,1,npoin);
 		free_dvector(vlm,1,npoin);
 		free_dvector(wxh,1,npoin);
 		free_dvector(ax,1,npoin);
 		free_dvector(axs,1,npoin);
 		free_dvector(gamas,1,npoin);
 		free_dvector(dist,1,npoin);
 		free_dvector(a2,1,npoin);
 		free_dvector(ffbe,1,npoin);
//     	free_dvector(dff,1,npoin);
 		if(nmethod>=2)
 		  {free_dvector(fthick0,1,npoin);
 		   free_dvector(fthick1,1,npoin);};
 		free_dmatrix(dff,1,npoin,0,0);
 		free_dvector(amtx,1,nnu);


     		     puts("\n");
     		     printf("\n '--------------------------------------------'");
     		     printf("\n '     MESHLESS3 IS ENDED SUCCESSFULLY        '");
     			 printf("\n '		 THANK YOU FOR USING MESHLESS3       '");
     		     printf("\n ''    ......... BYE,  BYE.........           '");
     		     printf("\n ''--------------------------------------------'");
     		     system("pause=nul");

	return 0;}

