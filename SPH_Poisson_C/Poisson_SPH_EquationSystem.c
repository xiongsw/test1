/*
 * Poisson_SPH_EquationSystem.c
 *
 * shape function meets linear continuity condition
 * (tensor product for weight function)
 *  Created on: Oct 9, 2013
 *      Author: Dr. Shangwu Xiong
 */
#include "user_defined.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

 struct constant_name {
   double constant1;
   double constant2;
}  poisson_const;

extern double bandv0(int, int, int, double *, double **, int *,int);

void coordinates_poisson(double, double,double **);
void boundary(double, double,int *, int *, double *, double *, double **);
void function1_poisson (double *,double **);
void weight (double,double *wxhi);
void weight_deriv_i (double,double,double *wxhi,double *dwxhi);
void weight_deriv_r (double,double,double *wxhi,double *dwxhi,double*ddwxhi);
void simple_output (double *, double **, double **);

void Poisson_SPH_Equation(double **coord,double **dax,double **dbetax,double **dbetay, double**gama,double**gama0,double **beta,\
  double**betas, int *nkode,	int *nkbe,int *ngama,double *fi,double *vlm,double *wxh,double *ax,double *gamas,double *axs,double *dist,\
		double *a2,double *ffbe,double **dff,double *amtx)
{

//

//c********************************************************c
//c********************************************************c
//c
//c
//c
     int i,j,k,iy,kk,ii=0,jj=0,nn,judge=0,nuvab,ntime,nttime,njudge,ik,ir=0,istep;
	 double ddx,ddy,xa,ya,xb,yb,xc,yc,hx,hy,h2;
     double radiusx,radiusy,volume,csxi,csyi,VxW,VyW,VxxW,VxyW,VyxW,VyyW,try0;
     double det,b11,b12,b21,b22,vdwxh1,vdwxh2,vwxhi,vwdbe1,vwdbe2,csxia,csyia;
     double sign,dmwxh1,dmwxh2,wxhj,xdir,Vxddw1,Vxddw2,Vddw,Vrxddw1,Vrxddw2,ydir,dwxhr,ddwxh,ddwxhj,ddwxhs;
	 double ddwxhi,dmwxh3,dmwxh4,csxib,csyib,rab2,ddw0,pK,pet,delta,bta,wxhi;
     double qb,wxhjx,wxhjy,dwxhx,dwxhy,ddwxhsi,ddwxhsj,Vrddw;

//c
     poisson_const.constant1=0.0;
     poisson_const.constant2=0.0;
//c

        ddx=(xmax-xmin)/(double)(npoin_x-1.);
        ddy=(ymax-ymin)/(double)(npoin_y-1.);
//c
       coordinates_poisson(ddx, ddy,coord);
       //c
               hx=alfa*ddx;
               hy=alfa*ddy;
               h2=hx*hy;
       //c
               eta=sqrt(10.*alfa/(7.*pi*h2));
       //c

       for (k=1;k<=npoin;k++){ gama0[1][k]=0.0;
       gama0[2][k]=0.0;
       gama[1][k]=0.0;
       gama[2][k]=0.0;
       nkode[k]=-1;
       }
//c
       function1_poisson (fi,coord);
//
// define constant boundary conditions
        boundary(ddx, ddy, nkbe, nkode, ffbe, dist, coord);

//c
//compute 1/volume
       for(i=1;i<=npoin;i++)
         {xa=coord[1][i];
          ya=coord[2][i];
          volume=0.0;
          for(j=1;j<=npoin;j++)
             {xb=coord[1][j];
              yb=coord[2][j];
              radiusx=fabs(xa-xb);
                csxi=radiusx/hx;
              if(csxi>=2.)continue;
//
                radiusy=fabs(ya-yb);
                csyi=radiusy/hy;
              if(csyi>=2.)continue;
// weight function
            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);

            wxh[j]=wxhjx*wxhjy;
            volume=volume+wxh[j];
		  }
          vlm[i]=1./volume;
//
         }
//c
// start to obtain parameters for shape functions and its derivatives
	for(i=1;i<=npoin;i++) {  
       xa=coord[1][i];
       ya=coord[2][i];
       VxW=0.0;
       VyW=0.0;
       VxxW=0.0;
       VxyW=0.0;
       VyxW=0.0;
       VyyW=0.0;
       for (j=1;j<=npoin;j++){ 
         wxh[j]=0.0;
         xb=coord[1][j];
         yb=coord[2][j];
         radiusx=fabs(xa-xb);
         csxi=radiusx/hx;
//
         if(csxi>=2.)continue;
           radiusy=fabs(ya-yb);
           csyi=radiusy/hy;
//
         if(csyi>=2.)continue;
            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);

          wxh[j]=wxhjx*wxhjy;
//c
          VxW=VxW+vlm[j]*(xb-xa)*wxh[j];
          VyW=VyW+vlm[j]*(yb-ya)*wxh[j];
          VxxW=VxxW+vlm[j]*(xa-xb)*(xa-xb)*wxh[j];
          VxyW=VxyW+vlm[j]*(xa-xb)*(ya-yb)*wxh[j];
          VyxW=VyxW+vlm[j]*(ya-yb)*(xa-xb)*wxh[j];
          VyyW=VyyW+vlm[j]*(ya-yb)*(ya-yb)*wxh[j];
//
	     } // 
//c
          det=VxxW*VyyW-VxyW*VyxW;
          beta[1][i]=(VyyW*VxW-VyxW*VyW)/det;
          beta[2][i]=(-VxyW*VxW+VxxW*VyW)/det;
//c
//c
          b11=0.0;
          b12=0.0;
          b21=0.0;
          b22=0.0;
          vdwxh1=0.0;
          vdwxh2=0.0;
          vwxhi=0.0;
          for(j=1;j<=npoin;j++) {   
            wxh[j]=0.0;
            xb=coord[1][j];
            yb=coord[2][j];
            radiusx=fabs(xa-xb);
            csxi=radiusx/hx;
            if(csxi>=2.)continue;
			
              radiusy=fabs(ya-yb);
              csyi=radiusy/hy;
            if(csyi>=2.)continue;

            weight_deriv_i (csxi,hx,&wxhjx,&dwxhx);
            weight_deriv_i (csyi,hy,&wxhjy,&dwxhy);

            wxh[j]=wxhjx*wxhjy;
            dwxhx=dwxhx*wxhjy;
            dwxhy=dwxhy*wxhjx;
//c
            try0=1.+beta[1][i]*(xa-xb)+beta[2][i]*(ya-yb);
//c
            b11=b11-vlm[j]*wxh[j]*(1.+2.*beta[1][i]*(xa-xb)+beta[2][i]*(ya-yb))-vlm[j]*(xa-xb)*(dwxhx*(xa-xb))*try0;
	        b12=b12-vlm[j]*wxh[j]*beta[1][i]*(ya-yb)-vlm[j]*(ya-yb)*(dwxhx*(xa-xb))*try0;
            b21=b21-vlm[j]*wxh[j]*beta[2][i]*(xa-xb)-vlm[j]*(xa-xb)*(dwxhy*(ya-yb))*try0;
            b22=b22-vlm[j]*wxh[j]*(1.+2.*beta[2][i]*(ya-yb)+beta[1][i]*(xa-xb))-vlm[j]*(ya-yb)*(dwxhy*(ya-yb))*try0;
//c
            vdwxh1=vdwxh1+vlm[j]*dwxhx*(xa-xb)*try0;
            vdwxh2=vdwxh2+vlm[j]*dwxhy*(ya-yb)*try0;
            vwxhi=vwxhi+vlm[j]*wxh[j]*try0;
//c

          }  //
//c
           dbetax[1][i]=(VyyW*b11-VyxW*b12)/det;
           dbetax[2][i]=(-VxyW*b11+VxxW*b12)/det;
           dbetay[1][i]=(VyyW*b21-VyxW*b22)/det;
           dbetay[2][i]=(-VxyW*b21+VxxW*b22)/det;
//c
           ax[i]=1./vwxhi;
//c
//c
          vwdbe1=0.0;
          vwdbe2=0.0;
       for(j=1;j<=npoin;j++)    //
         { xb=coord[1][j];
           yb=coord[2][j];
           radiusx=fabs(xa-xb);
           csxi=radiusx/hx;
//
           if(csxi>=2.)continue;
            radiusy=fabs(ya-yb);
            csyi=radiusy/hy;
 //
           if(csyi>=2.)continue;
		   
            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);

        wxh[j]=wxhjx*wxhjy;
        vwdbe1=vwdbe1+vlm[j]*wxh[j]*(dbetax[1][i]*(xa-xb)+dbetax[2][i]*(ya-yb)+beta[1][i]);
        vwdbe2=vwdbe2+vlm[j]*wxh[j]*(dbetay[1][i]*(xa-xb)+dbetay[2][i]*(ya-yb)+beta[2][i]);

       }  ////
//c
       dax[1][i]=-ax[i]*(vdwxh1+vwdbe1)/vwxhi;
       dax[2][i]=-ax[i]*(vdwxh2+vwdbe2)/vwxhi;
//c

        } //
// done for shape functions and its derivatives
//start to obtain Integration Correction coefficients
        nuvab=0;
        for(i=1;i<=npoin;i++)   
        {                         
           if(nkode[i]<=0){
		 	nuvab=nuvab+1;
			ngama[nuvab]=i;}
        }

        nttime=6*alfa;
     for(ntime=1;ntime<=nttime;ntime++)
     {//printf("ntime=%d, nuvab=%d\n",ntime,nuvab);
       for(ii=1;ii<=nuvab;ii++)  //
        {i=ngama[ii];
        xa=coord[1][i];
        ya=coord[2][i];
        gama[1][i]=0.0;
        gama[2][i]=0.0;
        for (j=1;j<=npoin;j++)  //
       { xb=coord[1][j];
        yb=coord[2][j];
       radiusx=fabs(xa-xb);
       csxi=radiusx/hx;
  //     
        if(csxi>=2.)continue;
          radiusy=fabs(ya-yb);
          csyi=radiusy/hy;
//       
          if(csyi>=2.)continue;
            weight_deriv_i (csxi,hx,&wxhjx,&dwxhx);
            weight_deriv_i (csyi,hy,&wxhjy,&dwxhy);

         wxh[i]=wxhjx*wxhjy;
         dwxhx=dwxhx*wxhjy;
         dwxhy=dwxhy*wxhjx;
//c
       dmwxh1=(dwxhx*(xb-xa)*ax[j]+wxh[i]*dax[1][j])*(1.+\
    beta[1][j]*(xb-xa)+beta[2][j]*(yb-ya))+wxh[i]*ax[j]*\
	(dbetax[1][j]*(xb-xa)+dbetax[2][j]*(yb-ya)+beta[1][j]);
       dmwxh2=(dwxhy*(yb-ya)*ax[j]+wxh[i]*dax[2][j])*(1.+\
      beta[1][j]*(xb-xa)+beta[2][j]*(yb-ya))+wxh[i]*ax[j]*\
	(dbetay[1][j]*(xb-xa)+dbetay[2][j]*(yb-ya)+beta[2][j]);
//c
        wxhj=wxh[i]*ax[j]*(1.+beta[1][j]*(xb-xa)+beta[2][j]*(yb-ya));
        gama[1][i]=gama[1][i]+gama0[1][j]*vlm[j]*wxhj-vlm[j]*dmwxh1;
        gama[2][i]=gama[2][i]+gama0[2][j]*vlm[j]*wxhj-vlm[j]*dmwxh2;
//
           }  //
//     
        judge=nkode[i];
        if(judge<0)continue;

       for(j=1;j<=npbe;j++)   // 
       { jj=nkbe[j];
        xb=coord[1][jj];
        yb=coord[2][jj];
       radiusx=fabs(xa-xb);
       csxi=radiusx/hx; 
       if(csxi>=2.)continue;
	   
         radiusy=fabs(ya-yb);
         csyi=radiusy/hy;     
       if(csyi>=2.)continue;
	   
            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);

          wxh[jj]=wxhjx*wxhjy;
//c
        njudge=nkode[jj];
        sign=1.;
        Ab=dist[j];
        if(njudge<2||njudge>3)sign=-1.;
        nn=njudge-njudge/2*2;
	    ydir=nn*sign;
        xdir=(1.-nn)*sign;
        wxhj=wxh[jj]*ax[jj]*(1.+beta[1][jj]*(xb-xa)+beta[2][jj]*(yb-ya));
        gama[1][i]=gama[1][i]+Ab*xdir*wxhj;
        gama[2][i]=gama[2][i]+Ab*ydir*wxhj;

          }//
//c

        }//
//c
         for(i=1;i<=npoin;i++)   
           {gama0[1][i]=gama[1][i];
	 	    gama0[2][i]=gama[2][i];}
		    
		}
// done for Integration Correction coefficients
//
// start to compute stabilization coefficients
     for(i=1;i<=npoin;i++)
	 {xa=coord[1][i];
      ya=coord[2][i];
//c
      VxW=0.0;
      VyW=0.0;
      VxxW=0.0;
      VxyW=0.0;
      VyxW=0.0;
      VyyW=0.0;
//c
      Vxddw1=0.0;
      Vxddw2=0.0;
	 Vddw=0.0;
//c
      Vrddw=0.0;
      Vrxddw1=0.0;
      Vrxddw2=0.0;
//c
   for(j=1;j<=npoin;j++)
   {  xb=coord[1][j];
      yb=coord[2][j];
      radiusx=fabs(xa-xb);
      csxi=radiusx/hx;
      if(csxi>=2.)continue;
      radiusy=fabs(ya-yb);
      csyi=radiusy/hy;
      if(csyi>=2.)continue;

      weight_deriv_r (csxi,hx,&wxhjx,&dwxhx,&dwxhr);
      weight_deriv_r (csyi,hy,&wxhjy,&dwxhy,&ddwxh);
//c
      ddwxhj=dwxhr*wxhjy+ddwxh*wxhjx;
      VxW=VxW-vlm[j]*(xb-xa)*ddwxhj;
      VyW=VyW-vlm[j]*(yb-ya)*ddwxhj;
      VxxW=VxxW+vlm[j]*(xb-xa)*(xb-xa)*ddwxhj;
      VxyW=VxyW+vlm[j]*(xb-xa)*(yb-ya)*ddwxhj;
      VyxW=VyxW+vlm[j]*(yb-ya)*(xb-xa)*ddwxhj;
      VyyW=VyyW+vlm[j]*(yb-ya)*(yb-ya)*ddwxhj;
//c
      Vxddw1=Vxddw1+vlm[j]*(xb-xa)*ddwxhj;
      Vxddw2=Vxddw2+vlm[j]*(yb-ya)*ddwxhj;
      Vddw=Vddw+vlm[j]*ddwxhj;
//c
      rab2=radiusx*radiusx+radiusy*radiusy;
      Vrddw=Vrddw+vlm[j]*ddwxhj*rab2;
      Vrxddw1=Vrxddw1+vlm[j]*ddwxhj*rab2*(xb-xa);
      Vrxddw2=Vrxddw2+vlm[j]*ddwxhj*rab2*(yb-ya);
   }
//c
      det=VxxW*VyyW-VxyW*VyxW;
      betas[1][i]=(VyyW*VxW-VyxW*VyW)/det;
      betas[2][i]=(-VxyW*VxW+VxxW*VyW)/det;
//c
      ddw0=-eta*(3./(hx*hx)+3./(hy*hy))*eta;
      gamas[i]=(betas[1][i]*Vxddw1+betas[2][i]*Vxddw2+Vddw)/(-ddw0*vlm[i]);
/*c
c      Computations found that the values of axs at the corners
c          is too large compared with that at any other points
c          if alfa=1.0, which resulted in very bad results.
c      Therefore, we need restrict the minimium of axs.
c*/
      try0=Vrddw+betas[1][i]*Vrxddw1+betas[2][i]*Vrxddw2;
      if(fabs(try0)<0.01)
        {axs[i]=0.0;}
      else
       {axs[i]=4./try0;}
//      axs[i]=0.0;
//c

//      printf("i=%d, Vrddw=%f,  %f\n",i,Vrddw,axs[i]);
      }
//
//c
       pK=pi*2.0e9/3.*alfa;      //!pK=pi*200./3.*alfa
       pet=pi/3.*alfa*sqrt(hx*hy)/3.;
       if((imod<=1)&&(icase>1))pet=pet*3.;
       if(npoin<441)pet=pet*npoin/441.;
//c
       printf("nnu%d =ok\n",nnu);
//c
       istep=0;
       istep=istep+1;
//c
// start to build equation system 
//c
       ik=0;
       for(i=1;i<=npoin;i++)   
        {xa=coord[1][i];
        ya=coord[2][i];
        dff[i][0]=0.0;
        ngama[i]=0;
        ii=i-4*npoin_y;
        jj=i+4*npoin_y;
        if(ii<=1)ii=1;
        if(jj>=npoin)jj=npoin;
        for(k=ii;k<=jj;k++)a2[k]=0.0;
       for(j=1;j<=npoin;j++) 
       {xb=coord[1][j];
        yb=coord[2][j];
        radiusx=fabs(xa-xb);
        csxi=radiusx/hx;
//c
       if(csxi>=4.)continue;
         radiusy=fabs(ya-yb);
         csyi=radiusy/hy;
       if(csyi>=4.)continue;
//c
       if(csxi<2.&& csyi<2.)
	   {
            weight_deriv_r (csxi,hx,&wxhjx,&dwxhx,&dwxhr);
            weight_deriv_r (csyi,hy,&wxhjy,&dwxhy,&ddwxh);
//            weight (csxi,&wxhjx);
//            weight (csyi,&wxhjy);
           wxh[j]=wxhjx*wxhjy;
//c
           dff[i][0]=dff[i][0]+fi[j]*vlm[j]*vlm[i]*(wxh[j]*ax[j]*(1.+beta[1][j]*(xb-xa)+beta[2][j]*(yb-ya)));

           delta=1.0;
           if(i!=j)delta=0.0;
           ddwxhj=ddwxh*wxhjx+dwxhr*wxhjy;
           ddwxhs=axs[j]*(1.+gamas[j]*delta+betas[1][j]*(xb-xa)+betas[2][j]*(yb-ya))*ddwxhj;
           dff[i][0]=dff[i][0]-pet*h2*fi[j]*vlm[j]*vlm[i]*ddwxhs;

	      } // 
//c
//c
       try0=0.0;
       if(j<=i)
	   {
        for(k=1;k<=npoin;k++)  //
         { xc=coord[1][k];
           yc=coord[2][k];
//c
           radiusx=fabs(xc-xa);
           csxia=radiusx/hx;
         if(csxia>=2.)continue;
          radiusy=fabs(yc-ya);
          csyia=radiusy/hy;
         if(csyia>=2.)continue;

         radiusx=fabs(xc-xb);
         csxib=radiusx/hx;
        if(csxib>=2.)continue;
         radiusy=fabs(yc-yb);
         csyib=radiusy/hy;
        if(csyib>=2.)continue;
//c
        weight_deriv_r (csxia,hx,&wxhjx,&dwxhx,&dwxhr);
        weight_deriv_r (csyia,hy,&wxhjy,&dwxhy,&ddwxh);

        wxh[i]=wxhjx*wxhjy;
        dwxhx=dwxhx*wxhjy;
        dwxhy=dwxhy*wxhjx;
//c
       delta=1.0;
       if(i!=k)delta=0.;
       bta=1.+beta[1][k]*(xc-xa)+beta[2][k]*(yc-ya);
       dmwxh1=(dwxhx*(xc-xa)*ax[k]+wxh[i]*dax[1][k])*bta+wxh[i]*ax[k]*\
	        (dbetax[1][k]*(xc-xa)+dbetax[2][k]*(yc-ya)+beta[1][k])+\
	        gama[1][i]/vlm[i]*delta-gama[1][k]*wxh[i]*ax[k]*bta;
       dmwxh2=(dwxhy*(yc-ya)*ax[k]+wxh[i]*dax[2][k])*bta+wxh[i]*ax[k]*\
	        (dbetay[1][k]*(xc-xa)+dbetay[2][k]*(yc-ya)+beta[2][k])+\
	        gama[2][i]/vlm[i]*delta-gama[2][k]*wxh[i]*ax[k]*bta;
//c
       ddwxhi=ddwxh*wxhjx+dwxhr*wxhjy;
       ddwxhsi=axs[k]*(1.+gamas[k]*delta+betas[1][k]*(xc-xa)+betas[2][k]*(yc-ya))*ddwxhi;
//c

            weight_deriv_r (csxib,hx,&wxhjx,&dwxhx,&dwxhr);
            weight_deriv_r (csyib,hy,&wxhjy,&dwxhy,&ddwxh);
		
        wxh[j]=wxhjx*wxhjy;
        dwxhx=dwxhx*wxhjy;
        dwxhy=dwxhy*wxhjx;
//c
       delta=1.0;
       if(j!=k)delta=0.;
       bta=1.+beta[1][k]*(xc-xb)+beta[2][k]*(yc-yb);
       dmwxh3=(dwxhx*(xc-xb)*ax[k]+wxh[j]*dax[1][k])*bta+wxh[j]*ax[k]*\
	        (dbetax[1][k]*(xc-xb)+dbetax[2][k]*(yc-yb)+beta[1][k])+\
	        gama[1][j]/vlm[j]*delta-gama[1][k]*wxh[j]*ax[k]*bta;
       dmwxh4=(dwxhy*(yc-yb)*ax[k]+wxh[j]*dax[2][k])*bta+wxh[j]*ax[k]*\
	        (dbetay[1][k]*(xc-xb)+dbetay[2][k]*(yc-yb)+beta[2][k])+\
	        gama[2][j]/vlm[j]*delta-gama[2][k]*wxh[j]*ax[k]*bta;
//c
       ddwxhj=ddwxh*wxhjx+dwxhr*wxhjy;
       ddwxhsj=axs[k]*(1.+gamas[k]*delta+betas[1][k]*(xc-xb)+betas[2][k]*(yc-yb))*ddwxhj;
//
       try0=try0+(dmwxh1*dmwxh3+dmwxh2*dmwxh4)*vlm[k]+ddwxhsi*ddwxhsj*pet*h2*vlm[k];
       }//
//c
          a2[j]=a2[j]+try0*vlm[i]*vlm[j];
	  } // condition: if(j<=i)ends  
//c
         try0=0.0;
//c
        judge=nkode[i];
        if(judge<0)continue;  //  
          njudge=nkode[j];
        if(njudge<0)continue; //  

        for(k=1;k<=npbe;k++) // 
        {kk=nkbe[k];
         xc=coord[1][kk];
         yc=coord[2][kk];
//c
        radiusx=fabs(xc-xa);
        csxia=radiusx/hx;
        if(csxia>=2.)continue;
        radiusy=fabs(yc-ya);
        csyia=radiusy/hy;
        if(csyia>=2.)continue;

        radiusx=fabs(xc-xb);
        csxib=radiusx/hx;
        if(csxib>=2.)continue;
        radiusy=fabs(yc-yb);
        csyib=radiusy/hy;
        if(csyib>=2.)continue;
//c
            weight (csxia,&wxhjx);
            weight (csyia,&wxhjy);
         wxh[i]=wxhjx*wxhjy;
//c
            weight (csxib,&wxhjx);
            weight (csyib,&wxhjy);
         wxh[j]=wxhjx*wxhjy;
//c
        Ab=dist[k];
//c
        wxhi=wxh[i]*ax[kk]*(1.+beta[1][kk]*(xc-xa)+beta[2][kk]*(yc-ya));
        wxhj=wxh[j]*ax[kk]*(1.+beta[1][kk]*(xc-xb)+beta[2][kk]*(yc-yb));
        try0=try0+pK*Ab*wxhi*wxhj;
       } //
        a2[j]=a2[j]+try0*vlm[i]*vlm[j];
//c
      } //
//
      if(judge>=0)
       { for (j=1;j<=npbe; j++)  //
        {jj=nkbe[j];
        xb=coord[1][jj];
        yb=coord[2][jj];
        radiusx=fabs(xa-xb);
        csxi=radiusx/hx;
        if(csxi>=2.)continue;
        radiusy=fabs(ya-yb);
        csyi=radiusy/hy;
        if(csyi>=2.)continue;

            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);

          wxh[jj]=wxhjx*wxhjy;
//c
        njudge=nkode[jj];
        sign=1.;
        Ab=dist[j];
        if(njudge<2||njudge>3)sign=-1.;
        nn=njudge-njudge/2*2;
	    ydir=nn*sign;
        xdir=(1.-nn)*sign;
        qb=.5*(cos(pi*xb)*sin(pi*yb)*xdir+sin(pi*xb)*cos(pi*yb)*ydir)/pi;
        dff[i][0]=dff[i][0]+qb*Ab*vlm[i]*wxh[jj]*ax[jj]*(1.+beta[1][jj]*(xb-xa)+beta[2][jj]*(yb-ya));
//c
       }//
//c
        for(j=1;j<=npbe;j++)   //
        { jj=nkbe[j];
          xb=coord[1][jj];
          yb=coord[2][jj];
        radiusx=fabs(xa-xb);
        csxi=radiusx/hx;
        if(csxi>=2.)continue;
          radiusy=fabs(ya-yb);
          csyi=radiusy/hy;
        if(csyi>=2.)continue;
		
            weight (csxi,&wxhjx);
            weight (csyi,&wxhjy);
          wxh[jj]=wxhjx*wxhjy;
		  
        Ab=dist[j];
        dff[i][0]=dff[i][0]+ffbe[j]*Ab*vlm[i]*wxh[jj]*ax[jj]*(1.+beta[1][jj]*(xb-xa)+beta[2][jj]*(yb-ya))*pK;
       } //
	  } // condition "if(judge>=0)" ends
//	  
         iy=0;
         ii=i-4*npoin_y;
         if(ii<=1)ii=1;
         for(j=ii;j<=i;j++)  // 
         { if(a2[j]!=0.)iy=1;
           if(iy>0)
             {ik=ik+1;
              amtx[ik]=a2[j];
              if(j==i)ngama[j]=ik;};
          }//
//c
   }//
//c
       printf("system equation is ready and start to solve  it\n\n");

	     try0=bandv0(npoin,0,ik, amtx,dff,ngama,ir);
//c
	     printf(" system equation is solved\n\n");
//c       output the final results
	     simple_output (fi, coord, dff);
}

void coordinates_poisson(double deltxi, double deltyi,double **coord)
{int i,j,k; double xk,yk;
    for (i=1;i<=npoin_x;i++)     
     for (j=1;j<=npoin_y;j++)     
      { k=(i-1)*npoin_y+j;
        xk=(i-1.)*deltxi+xmin;
        yk=(j-1.)*deltyi+ymin;
        coord[1][k]=xk;
        coord[2][k]=yk;};

};
void function1_poisson (double *fi,double **coord){
	int i,j,k;
	for (i=1;i<=npoin_x;i++)
    for (j=1;j<=npoin_y;j++)
     { k=(i-1)*npoin_y+j;
// a known function fi
         fi[k]=sin(pi*coord[1][k])*sin(pi*coord[2][k])+poisson_const.constant1;
    	 };};

void boundary(double deltxi, double deltyi,int *nkbe, int *nkode, double *ffbe, double *dist, double **coord)
{
	// define constant boundary conditions
	int i,j,k,ii,jj,kk,judge;
	double xa,ya,xj,xk,yj,yk,radiu1,radiu2,djudge,hx,hy,hxj,hyj;

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
	          ffbe[ii]=poisson_const.constant1;
	          ffbe[jj]=poisson_const.constant1;
		 	}
	//c
	        for (i=1;i<=npoin_x;i++)   //
	        { j=(i-1)*npoin_y+1;
	          k=j+npoin_y-1;
	          kk=npbe-npoin_y-i+1;
	          nkode[j]=1;
	          nkode[k]=3;
	          nkbe[i]=j;
	          nkbe[kk]=k;
	          ffbe[i]= poisson_const.constant1;
	          ffbe[kk]= poisson_const.constant1;
			}
	//c
	        nkbe[npbe+1]=1;
			xj=coord[1][1]-coord[1][1+npoin_y];
	        xk=coord[2][1]-coord[2][1+npoin_y];
	        dist[1]=0.5*sqrt(xj*xj+xk*xk);
	        for(i=2;i<=npbe;i++)  // distance of 2 neighboring boundary nodes
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
            hx=alfa*deltxi;
            hy=alfa*deltyi;
            hxj=2.*hx;
            hyj=2.*hy;

	        for (i=1;i<=npoin; i++)
	        { judge=nkode[i];

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

	          }
};


   void weight(double csi, double *wxhi)
 {double csi2, csi3;
   csi2=csi*csi;
   csi3=csi2*csi;
       if(csi<=1.){*wxhi=eta*(1.-1.5*csi2+0.75*csi3);
		    		   }
        else
          {*wxhi=eta*(2.-3.*csi+1.5*csi2-0.25*csi3);
           		  }
   }

   
   void weight_deriv_i(double csi, double hi, double *wxhi,double *dwxhi)
 {double csi2, csi3;
   csi2=csi*csi;
   csi3=csi2*csi;
       if(csi<=1.)
	       {*wxhi=eta*(1.-1.5*csi2+0.75*csi3);
		    *dwxhi=eta*(-3.+2.25*csi)/(hi*hi);
		   }
        else
          {*wxhi=eta*(2.-3.*csi+1.5*csi2-0.25*csi3);
           *dwxhi=eta*(-3./csi+3.-0.75*csi)/(hi*hi);
		  }
   }
	 
   void weight_deriv_r(double csi, double hi, double *wxhi,double *dwxhi,double*ddwxhi)
 {double csi2, csi3;
   csi2=csi*csi;
   csi3=csi2*csi;
       if(csi<=1.)
	       {*wxhi=eta*(1.-1.5*csi2+0.75*csi3);
		    *dwxhi=eta*(-3.+2.25*csi)/(hi*hi);
		    *ddwxhi=eta*(-3.+4.5*csi)/(hi*hi);
		   }
        else if(csi<=2.)
          {*wxhi=eta*(2.-3.*csi+1.5*csi2-0.25*csi3);
           *dwxhi=eta*(-3./csi+3.-0.75*csi)/(hi*hi);
		   *ddwxhi=eta*(3.-1.5*csi)/(hi*hi);
          }
        else
        {*wxhi=0.0;*dwxhi=0.;*ddwxhi=0.0;};
   }

   void simple_output (double *fi, double **coord, double **dff){

   	FILE *fp_out;
   	double xa,ya,pi2,ftrue,error1;
   	int i,j,k;
	      if((fp_out=fopen("try01.out","w"))==NULL)
			   {
                   fprintf(stderr, "can not open the output file : try01.out");
                   exit(1);
                  };

       pi2=pi*pi;
//c
//c      ftrue: the analytic solution
       fprintf(fp_out,"TITLE = 'solution of Poisson Equation'\n");
       fprintf(fp_out," No.   Coord-x,  Coord-y   Result_SPh   Exact_Solution\n");
       for(j=1;j<=npoin_y;j++)    ///
       for(i=1;i<=npoin_x;i++)   //
        {
			k=(i-1)*npoin_y+j;
		       ftrue=0.5*fi[k]/pi2+poisson_const.constant1;
		       error1=ftrue-dff[k][0];

//c dff[k]: solution
//c fi[k]: prescribed function in Poisson equation
          fprintf(fp_out, " %d   %f   %f    %12.4f    %f\n", k,coord[1][k],coord[2][k],dff[k][0],ftrue);
	   }
//c
       fclose(fp_out);
//c
      if((fp_out=fopen("Xiong_SPH_P.dat","w"))==NULL)
                   {fprintf(stderr, "can not open the output file : Xiong_SPH_P");
                   exit(1);
                  };

        fprintf(fp_out,"TITLE = 'solution of Poisson Equation'\n");
		fprintf(fp_out,"VARIABLES = X, Y, Result \n");
		fprintf(fp_out,"ZONE T=All, I=%d", npoin_x); fprintf(fp_out,", J=%d",npoin_y); fprintf(fp_out,", F=POINT\n");

      for(j=1;j<=npoin_y;j++)
      for(i=1;i<=npoin_x;i++)
       {k=(i-1)*npoin_y+j;
	    xa=coord[1][k];
	    ya=coord[2][k];
          Ab=dff[k][0];
          if(Ab<=0.0)Ab=0.0;
          fprintf(fp_out,"%f   %f   %12.5f\n", xa,ya,Ab);
	  }

      fclose(fp_out);


   }
