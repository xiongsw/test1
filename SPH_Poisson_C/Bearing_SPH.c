/*
 * Bearing_SPH.c
 *(the main part is deleted and not open to the public)
 *  Created on: Oct 10, 2014
 *      Author: Dr. Shangwu Xiong
 */

#include "user_defined.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

extern double bandv0(int, int, int, double *, double **, int *,int);

void Bearing_main1(double **coord,double **dax,double **dbetax,double **dbetay, double**gama,double**gama0,double **beta,\
		int *nkode,	int *nkbe,int *ngama,double *fi,double *vlm,double *wxh,double *ax,double *fthick1,double *fthick0,double *dist,double *a2,\
		double *ffbe,double **dff,double *amtx)
{

//	using namespace std;
	FILE *fp_out;
//c********************************************************c
//c********************************************************c
//c
//c
//c
     int i,j,k,nelem_x,nelem_y,ix,iy,kk,m,ii,jj,nn,judge,nuvab,ntime,nttime,njudge,ik,ir=0,istep;
	 double ddx,ddy,xk,yk,c,xcm,xcn,ycm,ycn,xi,yi,xj,yj,xa,ya,xb,yb,xc,yc,xx,yy,hx,hy,h2,hxj,hyj,radiu1,radiu2,rmin,radius;
     double radiusx,radiusy,const1,yfl,djudge,volume,csxi,csxi2,csxi3,csyi,csyi2,csyi3,VxW,VyW,VxxW,VxyW,VyxW,VyyW,try0;
     double det,b11,b12,b21,b22,vdwxh1,vdwxh2,vwxhi,dwxh,vwdbe1,vwdbe2,csxia,csxia1,csxia2,csxia3,csyia,csyia1,csyia2,csyia3;
     double sign,dmwxh1,dmwxh2,wxhj,xdir,Vxddw1,Vxddw2,Vddw,Vrxddw1,Vrxddw2,error,ydir,VRDDW,dwxhr,ddwxh,ddwxhj,ddwxhs;
	 double ddwxhi,ddwxhSI,ddwxhSJ,dmwxh3,dmwxh4,csxib,csxib2,csxib3,csyib,csyib2,csyib3,rab2,ddw0,pK,pet,delta,csib,bta,wxhi;
     double qb,pi2,ftrue,wxhjx,wxhjy,dwxhx,dwxhy;
     double sita,sitaj,tmp,hj,Bearing_main1_return=0;
//c
//c
        omega=omega1;
	    printf(" eccen1[1]= %f eccen1[2]=%f omega=%f", eccen1[1], eccen1[2], omega);
//c

        ddx=(xmax-xmin)/(double)(npoin_x-1.);
        ddy=(ymax-ymin)/(double)(npoin_y-1.);
//c
        nelem_x=npoin_x-1;
        nelem_y=npoin_y-1;
//c
        for (i=1;i<=npoin_x;i++)     
        for (j=1;j<=npoin_y;j++)    
        {k=(i-1)*npoin_y+j;
        xk=(i-1.)*ddx+xmin;
        yk=(j-1.)*ddy+ymin;
        coord[1][k]=xk;
        coord[2][k]=yk;
		dff[k][0]=0.0;
		} 
//c
        vmean=((omega1+omega2)/2.)/60.*2.*pi*rshaft;
        for(i=1;i<=npoin_x;i++) 
        for(j=1;j<=npoin_y;j++) 
         {k=(i-1)*npoin_y+j;
          sita=coord[1][k]/rshaft;
//c film thickness at j-point
          hj=clrnce-eccen1[1]*sin(sita)-eccen1[2]*cos(sita);
          fthick1[k]=hj;
          if(i<=1){
            jj=(nelem_x-1)*npoin_y+j;
		    sitaj=coord[1][jj]/rshaft;}
          else
           {jj=k-npoin_y;
		    sitaj=coord[1][jj]/rshaft;}
//c
          tmp=clrnce-eccen1[1]*sin(sitaj)-eccen1[2]*cos(sitaj);
//c
          fi[k]=-12.*(hj-tmp)/ddx*vmean;
          gama0[1][k]=0.0;
          gama0[2][k]=0.0;
          nkode[k]=-1;
	     	}  
//c
        const1=0.0;
        for(i=1;i<=npoin_y;i++)
		{
          j=npoin-npoin_y+i;
          nkode[i]=4;
          nkode[j]=2;
          ii=npbe-i+1;
          jj=npoin_x+i;
          nkbe[ii]=i;
          nkbe[jj]=j;
//c
          ffbe[ii]=const1;
          ffbe[jj]=const1;
	 	}  
//c
        for (i=1;i<=npoin_x;i++)   
        { j=(i-1)*npoin_y+1;
          k=j+npoin_y-1;
          kk=npbe-npoin_y-i+1;
          nkode[j]=1;
          nkode[k]=3;
          nkbe[i]=j;
          nkbe[kk]=k;
          ffbe[i]=const1;
          ffbe[kk]=const1;
		}  
//c
        nkbe[npbe+1]=1;
		xj=coord[1][1]-coord[1][1+npoin_y];
        xk=coord[2][1]-coord[2][1+npoin_y];
        dist[1]=0.5*sqrt(xj*xj+xk*xk);
        for(i=2;i<=npbe;i++)  
         { j=nkbe[i];
           k=nkbe[i+1];
           xj=coord[1][j];
           yj=coord[2][j];
           xk=coord[1][k];
           yk=coord[2][k];
           radiu1=sqrt((xj-xk)*(xj-xk)+(yj-yk)*(yj-yk));
           k=nkbe[i-1];
           xk=coord[1][k];
           yk=coord[2][k];
           radiu2=sqrt((xj-xk)*(xj-xk)+(yj-yk)*(yj-yk));
           dist[i]=(radiu1+radiu2)*0.5;
		}  
//c
        hx=alfa*ddx;
        hy=alfa*ddy;
        h2=hx*hy;
//c
        eta=sqrt(10.*alfa/(7.*pi*h2));
//c
        hxj=2.*hx;
        hyj=2.*hy;
        for (i=1;i<=npoin; i++)  
        { judge=nkode[i];
          // 
		 if(judge<=0){
		   xa=coord[1][i];
		   ya=coord[2][i];

	        djudge=hxj-(xa-xmin);
	          if(djudge>0.) {nkode[i]=0;continue;}
	            djudge=hxj-(xmax-xa);
	          if(djudge>0.){nkode[i]=0;continue;}
	            djudge=hyj-(ya-ymin);
	          if(djudge>0.){nkode[i]=0;continue;}
	            djudge=hyj-(ymax-ya);
	          if(djudge>0.){nkode[i]=0;continue;}
		   }
//
          }
//c
//      *********here: the main part is deleted and not open to the public
//c
//c       output the final results
//c
 //      fp_out=fopen("try01.out","w");
//
	      if((fp_out=fopen("try01.out","w"))==NULL)
			   {
                             fprintf(stderr, "can not open the output file : try01.out");
                             exit(1);
                            };

       pi2=pi*pi;
//c
//c      ftrue: the analytic solution
      fprintf(fp_out," No.   Coord-x,  Coord-y   Pressure_SPh   Function_Poisson\n");
       for(j=1;j<=npoin_y;j++)    /// 
       for(i=1;i<=npoin_x;i++)   // 
        {
			k=(i-1)*npoin_y+j;
//c dff[k]: fluid pressure
//c fi[k]: function related with dh/dx
          fprintf(fp_out, " %d   %f   %f    %f    %f\n", k,coord[1][k]*180./(pi*rshaft),coord[2][k]/hlength,dff[k][0],fi[k]);
	   }

       fclose(fp_out);
//c
      fp_out=fopen("Xiong_SPH_P.dat","w");
//      write(10,*)'TITLE = "Fluid Pressure"'
//      write(10,*)'VARIABLES = "X", "Y", "P_bar"'
        fprintf(fp_out,"TITLE = 'Fluid Pressure'\n");
		fprintf(fp_out,"VARIABLES = X, Y, P_bar \n");
		fprintf(fp_out,"ZONE T=All, I=%d", npoin_x); fprintf(fp_out,", J=%d",npoin_y); fprintf(fp_out,", F=POINT\n");

      for(j=1;j<=npoin_y;j++)  
      for(i=1;i<=npoin_x;i++)  
       {k=(i-1)*npoin_y+j;
	    xa=coord[1][k]*180.0/pi;
	    ya=coord[2][k];
          Ab=dff[k][0];
          if(Ab<=0.0)Ab=0.0;
          fprintf(fp_out,"%f   %f   %12.5f\n", xa/rshaft,ya/hlength,Ab);
	  }

      fclose(fp_out);

//       return 0;
//    		   Bearing_main1_return;
//     	   }
}


