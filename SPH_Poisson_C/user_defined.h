/*
 * user_defined.h
 *
 *  Created on: Oct 9, 2014
 *      Author: Shangwu Xiong
 */

#ifndef USER_DEFINED_H_
#define USER_DEFINED_H_

#define pi 3.14159265

#define NR_END          1
#define FREE_ARG        char*
#define MAX_ARRAY       2000


#define mnode 4
#define mextf 480
#define mcont 101
#define mhist 10
#define mcavck 10



//{
     double xmin,xmax,ymin,ymax,alfa,ratio;
     int npoin_x,npoin_y,npoin,nwork,npbe,icase,imod,nnu,ii9,npoinx,npoiny;
//
      char filename1[40];

      double rshaft1,rshaft,clrnce,hlength,hmin,anglemis;

      int nsmth,nwavx,nwavy,n12,nfine,kincf,nok;
      double arough,aa,Ab,ala,alb,al,gamma9,dwave,sigma,mrough;
      double Pc_c,Pc_r,Pc_c9,Pc_r9,Pc_A1,Pc_A2,Pc_af1,Pc_af2,Pc_af3;


      double emod1,emod2,pois1,pois2,fcoef,emod,pois,eta0,eta;

      int maxnp,maxel;

      int numel,numnp,nnode,ndime,ndofn,ntlelm,nzlelm,ntnode,nznode,nint,oilsupply,idimen;


      int nrigid,infopw,irepeat,ncontact,nstep,ndynamic;


     double press0,ex0,ey0,pressinlet,rdinlet;


      double omega1,omega2,omega,vmean;


      int NP_Vcurve,NP_Fcurve,ntype_v;

      int nextf;

      double pomeg,extfx[mextf],extfy[mextf],dtime0,timef[mextf],timeomega[mextf],omegashift,cfx,cfy;


      int ntinp,ntout,nthis,ntmsg,nthouse,ntres,ntool,iboundary;


      double ccontact,ccc,pa,pc;


      int maxeq,nefn,neq,ntoteq;
      double  eccen1[3],eccen2[3];


      int kcontact,ncont;
      double frictx,fricty,frict,wa,wamax,watime;

      int nincmax,kinc,nattmax,nattmp,nitemax,iter,nmethod,ncavbc,icavbc,iterold,nsecond;
      double dtime,time,time0,stime,qcontr,dtmin,dtmax,olderr,newerr;

     int nchange,ncavck,ncavcki,nchg[mcavck];

      char fileres[20];
      int nrestart,kinc0,nresinc;

      int nlsmax,ilfail;
      double permls,permet,ampmx,etmxa,etmna,slp;

      double tolres,tolerr,tolpre,fmaxp,PressureMax;

      int njacob;

      int nfreq,notime;
      double otime[10];

      int nhist,indexh[mhist],nfreqh;
      double htime,plimit,hist[mhist];
      char hname[mhist];

#endif /* USER_DEFINED_H_ */
