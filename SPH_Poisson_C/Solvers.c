//#include <stdlib.h>
//#include "nrutil.c"

/* An efficient procedure for solving B = A.X is the LU-decomposition.
The LU-decomposition method first "decomposes" matrix A into A = L.U, where L and U are lower triangular and upper triangular matrices, respectively. 
 More precisely, if A is a n�n matrix, L and U are also n�n matrices with forms 
The lower triangular matrix L has zeros in all entries above its diagonal and the upper triangular matrix U has zeros in all entries below its diagonal.
 If the LU-decomposition of A = L.U is found, the original equation becomes B = (L.U).X. This equation can be rewritten as B = L.(U.X). 
Since L and B are known, solving for B = L.Y gives Y = U.X. Then, since U and Y are known, solving for X from Y = U.X yields the desired result.
 In this way, the original problem of solving for X from B = A.X is decomposed into two steps: 
1. Solving for Y from B = L.Y 
2. Solving for X from Y = U.X 

*/
//c********************************************************c
  double bandv0(int n, int m, int np, double *amtx, double **dff, int *jd,int ir)
   {
//c********************************************************c
//c
        int i,j,k,l,i0,j0,mi,mj,im1,jm1,ij,ik,jk,ii,jj,kk,mij;
        double bandv_return=1;
//c
//c
     for(i=1;i<=n;i++){
        i0=jd[i]-i;
        if(i!=1)
        {mi=jd[i-1]-i0+1;
        for(j=mi;j<=i;j++){
         j0=jd[j]-j;
         mj=1;
         if(j>1) mj=jd[j-1]-j0+1;
         mij=mi;
         if(mj>mi) mij=mj;
         ij=i0+j;
         jm1=j-1;
         for (k=mij;k<=jm1;k++)
          {
//
            if(mij>jm1) continue;
            ik=i0+k;
            kk=jd[k];
            jk=j0+k;
            amtx[ij]=amtx[ij]-amtx[ik]*amtx[kk]*amtx[jk];
//
          };
        if(j==i) continue;
          jj=jd[j];
          amtx[ij]=amtx[ij]/amtx[jj];
//
        for (k=0;k<=m;k++){
            dff[i][k]=dff[i][k]-amtx[ij]*amtx[jj]*dff[j][k];}         
//
        }
       }
//	   
        ii=i0+i;
        if(amtx[ii]==0.0) { ir=i; return ir*1.0;};
//
        for(k=0;k<=m;k++){
          dff[i][k]=dff[i][k]/amtx[ii];
        }
//
     }
//                
       for(l=2;l<=n;l++) {
        i=n-l+2;
        i0=jd[i]-i;
        mi=jd[i-1]-i0+1;
        im1=i-1;
        for (j=mi;j<=im1;j++){
          if(mi>im1) continue;
          ij=i0+j;
//
          for(k=0;k<=m;k++){
            dff[j][k]=dff[j][k]-amtx[ij]*dff[i][k];

           }

         }
       }
        
        ir=0;
//
        return bandv_return;                

//c
   }
//
   
