c******************************************************* c
c      ----    unt. No.1  of FAR2D ---- 1998.03.17  ---- c
c                                                        c
c------------------------------------------------------- c
c program for solving steady-state Plane Strain Rolling  c
c                                                        c
c                  by Rigid-Plastic FEM                  c
c            (slightly compressible material models)     c
c                                                        c
c  (2x2 Gauss(interior region), 5 Gauss(boundary line)   c
c  (1x1 Gauss(Reduced integration for imcompressiblity)  c
c                                                        c
c   edited by Dr.Xiong Shangwu(xiongsw@hotmail.com)      c
c********************************************************c
c
c			 
        program FEMroll2
c
	  character title*70
c
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
	  common /compr/ c_factor,beta,wet,force,pc,ps
        common /totalcpu/cputime,cpubt
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
c	  dinamic dimension of vectors
c
	  integer, allocatable:: ia(:)
	  real, allocatable:: a(:)
c
        write(*,*)
        write(*,*)
        write(*,*)
        write (*,998)
998     format(10(/),
     1          15X,'           WELCOME TO USE ROLL2	           ',//,
     2          15X,'--------------------------------------------',/,
     3          15X,'          FINITE ELEMENT PROGRAM            ',//,
     4          15X,'           OF 2D ROLLING PROCESS            ',//,
     6          15X,'                                            ',/,
     7          15X,'             1998.5.05 AT IST               ',/,
     8          15X,'--------------------------------------------')
        write(*,*)
        write(*,*)
        write(*,*)'       Please Press  <Enter>   to  Continue.......  '
112     read(*,111)n1
111     format(i4)
c
	  write(*,*)
c
c	  read input data
c
	  call input
c
        cpubt=secnds(0.0)
c
        open(7,file='FEMforce1.dat',form='formatted',status='unknown')
        write(7,*)'Npoin_Y, Npoin_X, Force, Torque,   Total-CPU-time(s)'
c
777     npoin_x=nelem_x+1
        npoin_y=nelem_y+1
        nelem=nelem_x*nelem_y
        npoin=npoin_x*npoin_y
c
c	  pointers for the vectors 
c 
	  m1=1
	  m2=m1+nelem*4
	  m3=m2+npoin*4
	  m4=m3+npoin*4
	  m5=m4+npoin
	  m6=m5+npoin
	  m7=m6+npoin*2
	  m8=m7+npoin*2
	  m9=m8+nelem
	  m10=m9+npoin*2
	  m11=m10+4
c
	  n1=1
	  n2=n1+npoin*2
	  n3=n2+npoin
	  n4=n3+npoin_x
	  n5=n4+nelem
	  n6=n5+nelem*ngaus_f
	  n7=n6+npoin*2
	  n8=n7+npoin*2
	  n9=n8+npoin*2
	  n10=n9+(npoin*2)*(npoin_y+2)*2
	  n11=n10+nelem*3*4
	  n12=n11+nelem*4
	  n13=n12+nelem
	  n14=n13+nelem
	  n15=n14+ngaus_f
	  n16=n15+ngaus_f
	  n17=n16+nelem
	  n18=n17+npoin*2
	  n19=n18+nelem*4
	  n20=n19+nelem*4
	  n21=n20+npoin*4
	  n22=n21+npoin
cc
c	  information for the user
c
	  write (*,1000)
	  memor=4*m11+4*n22
	  write (*,1001) m11,n22,memor
1000	  format (//,
     1   5x,'INPUT TABLE No. 0 *** MEMORY REQUIREMENTS ***')
1001    format (/,
     1   5x,'INTEGERS ................... =',i8,/,
	2   5x,'REALS ...................... =',i8,//,
	3   5x,'TOTAL MEMORY ............... =',i8,' bytes')
c
c	  dinamic allocation of ia and a
c
	  allocate (ia(m11),a(n22))
c
c	  initialize all the vectors 
c	
	  do 1 i=1,m11
	  ia(i)=0
1	  continue
c
        do 2 i=1,n22
	  a(i)=0.0
2       continue
c
	  call main (ia(m1),
     1             ia(m2),ia(m3),ia(m4),
	2			 ia(m5),ia(m6),ia(m7),
	3	         ia(m8),ia(m9),ia(m10),
	4             a(n1),a(n2),
     5             a(n3),
     6             a(n4),a(n5),
     7	         a(n6),a(n7),
     8	         a(n8),a(n9),a(n10),
     9             a(n11),a(n12),a(n13),a(n14),
     1             a(n15),a(n16),a(n17),a(n18),a(n19),
     2             a(n20),a(n21))
c
c	  deallocate ia and a
c
	  deallocate (ia,a)
c
        nelem_x=nelem_x+1
c
c if one would like to see the influence of number of nodes at x direction
c (please recover the following 2 commands)
c        nelem_xmax=34
c        if(nelem_x.le.nelem_xmax)goto 777
        nelem_x=3
        nelem_y=nelem_y+1
c
c if one would like to see the influence of number of nodes at y direction
c  (please recorver the following 2 commands)
c        nelem_ymax=4
c        if(nelem_y.le.nelem_ymax)goto 777
        close(7)
	  end

c*************************************************** c
c                                                    c
c     subroutine No.1  of FAR2D -- 1998.03 ---       c
c     function: main control program                 c
c	called by roll2, call other subroutines        c
c                                                    c
c****************************************************c
c			 
	  subroutine main (lnods,
     1                   nodei,nodlc,noden,
	2			       ndofn,nfixn,nsuvab,
     3				   nefric,jd,nofric,
     4                   coord,tangent,
     5                   dtl,
     6                   stress,shear_f,
	7			       a1,a2,vnode,a,
     8                   strnrte,strt_egas,strt_qv,strt_eqm,
     9                   vegaus_x,vefric,strn_eq,
     1                   vnode0,s2,s3,
     2                   s2_nodes,area)
c
c	  this subroutine controls roll2
c
	  character title*70
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
	  common /compr/ c_factor,beta,wet,force,pc,ps
        common /totalcpu/cputime,cpubt
        common /bem1/  nbem,nodlast(5),nsurf,inpt,imp,nbem1_s
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
	  dimension lnods(4,*),
     1            nodei(4,*),nodlc(4,*),noden(*),
	2			ndofn(*),nfixn(2,*),nsuvab(*),
	3			nefric(*),jd(*),nofric(4),
     4            coord(2,*),tangent(*),
     5            dtl(*),
     6            stress(*),shear_f(5,*),
	7			a1(*),a2(*),vnode(2,*),a(*),
	8            strnrte(3,4,*),strt_egas(4,*),
     9            strt_qv(*),strt_eqm(*),
     1            vegaus_x(*),vefric(*),strn_eq(*),
	2	        vnode0(2,*),s2(4,*),s3(4,*),
	3            s2_nodes(4,*),area(*)
	  dimension err_w(3),err_u(3)
c
        beta=1.0
	  m6=0
	  up=1.0
c
c	  build fem model
c
c	  topology
c
	  call topology (lnods)
c
c	  find neighbours
c
	  call neighbour (lnods,nodei,nodlc,noden)
c
c	  assign co-ordinates to nodal points
c
	  call coordinates (coord,tangent,dtl)
c
c	  set up the unknown variables
c	      
        call unknowns (lnods,ndofn,nfixn,nsuvab,
     1                 tangent)

c        
c	  set up the contact elements	
c
        call contact_elements (nefric,nofric)
c
c	  establish initial flow stress and interfacial friction stress
c
	  call flwsts (stress,shear_f)
c
c	  set up initial velocity field by G-function
c
        call initial_velocity (lnods,nodei,nodlc,noden,
	1                         nfixn,nsuvab,nefric,nofric,jd,
     2                         coord,tangent,stress,shear_f,
     3                         a1,a2,vnode,a)
c
c	  calculate deformation energy for initial velocity field
c
        call energy(lnods,nefric,coord,stress,shear_f,
     1              vnode,strnrte,strt_egas,strt_qv,
     2              strt_eqm,vegaus_x,vefric)
c
c	  modify flow stress
c
16	  call stress1(lnods,nefric,coord,stress,shear_f,vnode,
     1                     strt_eqm,strn_eq,vegaus_x)
c
c       calculate deformation energy
c
99      call  energy(lnods,nefric,coord,stress,shear_f,
     1               vnode,strnrte,strt_egas,strt_qv,
     2               strt_eqm,vegaus_x,vefric)

c        optimize velocity field
c
77      call minfi (lnods,nodei,nodlc,noden,
     1              nsuvab,nefric,nofric,jd,
     2              coord,tangent,stress,shear_f,
     3              a1,a2,vnode,a,
     4              strnrte,strt_egas,strt_qv,strt_eqm,strn_eq,
     5              vegaus_x,vefric,vnode0)
c
c	  write(*,*)'up=',up,beta
	  write(*,*)'up=',up,m6
	  m6=m6+1
	  if(m6.gt.3)then
	  err_w(1)=err_w(2)
	  err_w(2)=err_w(3)
	  err_w(3)=wet
	  err_u(1)=err_u(2)
	  err_u(2)=err_u(3)
	  err_u(3)=up
	  call bet(err_w,err_u)
	  else
	  err_u(m6)=up
	  err_w(m6)=wet
	  endif
c
	  if (up.ge.0.0001)goto 16
c
c       calculate stress and deviation of equilibrium of forces
c 
        call strs(lnods,coord,nodei,nodlc,noden,stress,
     1             strnrte,strt_egas,strt_qv,strt_eqm,
     2             s2,s3)
c
c       export the calculated results by this program
c
199     if(ps.gt.0.0001)goto 16
	  call prt(lnods,coord,stress,vnode,strt_egas,
     1                 strt_qv,strt_eqm,
     2	       	     strn_eq,s2,s3)
c
c       calculate the seperating force at nodal points
c
        call nodarea (lnods,coord,s2,s2_nodes,area)
c
c       Post-processing for animation.exe 
c       (which was developed by STM/IST)
c
	  call neutral2 (lnods,coord,tangent,stress,vnode,strnrte,
	1                strt_eqm,strt_qv,strn_eq,s2,s2_nodes)
c
992     write(*,*)
        write(*,*)   '****** far-2d ended ******'
        write(*,*)
	  write(*,*)' Thank you greatly for using Roll2'
        write(*,*)
	  write(*,*)'......      Bye     ............. '
        write(*,*)
        write(*,*)
        write(*,*)'        Please Press <Enter> to End .......  '
112     read(*,111)n1
111     format(i4)
c        if(n1.lt.0)goto 112
	    write(*,*)
        end

c************************************************c
c                                                c
c  ---- sub. no.2  of far2d ---- 1998.02 ----    c
c   function: read the input data                c
c   called by: main , read from: fi.dat          c
c                                                c
c************************************************c
        subroutine input
c
c	  this subroutine reads the input data
c	
c	  title - name of the project
c
c	  nelem     - total number of elements
c	  nelem_x   - number of elements in x-axis
c	  nelem_y   - number of elements in y-axis
c	  nelem_xi	- number of elements at the inlet of the plastic zone 
c       nelem_xf	- number of elements at the exit of the plastic zone
c
c	  npoin     - total number of nodal points
c       npoin_x   - number of nodal points in x-axis
c	  npoin_y   - number of nodal points in y-axis
c
c	  nnode     - number of nodes per element (4)
c
c       radius        - radius of the roll	
c	  vel_roll      - linear velocity of the roll
c       thick_initial - initial half thickness of the workpiece
c	  thick_final   - final half thickness of the workpiece
c       width         - width of the workpiece (plane strain width)
c	  f_tension     - front tension
c	  b_tension     - back tension
c
c	  stress = ga*((gb+gc*ef_strain)**gn)*(ef_strain_rate**gm)+g0	
c
c       fric_m - Prandtl friction factor (m)
c       fric_f - Coulomb friction coeficient (f)
c       
c       c_factor - compressibility factor
c
	  character title*70
c
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
	  common /compr/ c_factor,beta,wet,force,pc,ps   
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
        data sgaus/-0.57735026918963E0,0.57735026918963E0/,
     1       sgaus_f/-0.9061798459E0,-0.5384693101E0,0.0,
     2               0.5384693101E0,0.9061798459E0/,
	3       wgaus/2*1.0E0/,wgaus_f/0.2369268851E0,0.4786286705E0,
     4              0.5688888889E0,0.4786286705E0,0.2369268851E0/
c
c	  define constants
c
	  ngaus=4
	  ngaus_f=5       
	  nnode=4
c
c	  open the input data file
c
	  idat=10	
	  open(idat,file='FEMroll2d.dat')
        read (idat,'(a70)') title
        read (idat,*) nelem_x,nelem_y,nelem_xi,nelem_xf
        read (idat,*) radius,vel_roll,thick_initial,thick_final,width,
 	1                f_tension,b_tension
        read (idat,*) ga,gb,gc,gn,gm,g0
        read (idat,*) fric_m
	  read (idat,*) c_factor
c
c	  close the input data file
c
        close (unit=10)
c        
c
	  e11=0.0025
        if(nelem_x.le.15)then
	  e2=0.00005*(nelem_x-nelem_xi-nelem_xf)
        else
        e2=0.00005*15.
        endif
c
c In fact one can give a positive integer to nelem_xi and nelem_xf,
c  this program works well. But in order to compare the results
c   with that calculated by RKPMroll_2D.for, we have to set them zero.
c
        nelem_xi=0
        nelem_xf=0
	  return
	  end	       

c************************************************c
c                                                c
c  ---- sub. no.3  of far2d ---- 1998.02 ----    c
c   function: define no. of node for each ele.   c
c   called by: main                              c
c                                                c
c************************************************c
        subroutine topology (lnods)
c
c	  this subroutine builds the topology of the fem model
c	  
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
        dimension lnods(4,*)
c
        do 100 ielem=1,nelem_x
        do 100 jelem=1,nelem_y
        knode=(ielem-1)*npoin_y+jelem
        kelem=(ielem-1)*nelem_y+jelem
        lnods (1,kelem)=knode
        lnods (2,kelem)=knode+npoin_y
        lnods (3,kelem)=knode+npoin_y+1
        lnods (4,kelem)=knode+1
100	  continue
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. no.4  of far2d ---- 1998.02 ----    c
c   function: inquire the node and ele.          c
c   called by: main , call iclear                c
c                                                c
c************************************************c
        subroutine neighbour (lnods,
     1                        nodei,nodlc,noden)
c
c	  this subroutine detects and stores the neighbours
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
c
	  dimension lnods(4,*),
     1            nodei(4,*),nodlc(4,*),noden(*)
c
c       knelem,nnelem:       Index and Maxium number of elements containing
c                             ipoin-th node(here nnelem=4, no used).
c       noden(ipoin):        Number of Elements Containing ipoin-th Node
c       nodei(knelem,ipoin): Index of jnelem-th Element containing ipoin-th Node
c       nodlc(knelem,ipoin): Local Index of ipoin-th node in knelem-th Element
c
        do 200 ipoin=1,npoin
        knelem=0
        do 190 jelem=1,nelem
        do 190 inode=1,nnode
        jpoin=lnods(inode,jelem)
	  if(jpoin.ne.ipoin) go to 190
        knelem=knelem+1
        nodei(knelem,ipoin)=jelem
        nodlc(knelem,ipoin)=inode
190     continue
        noden(ipoin)=knelem
200     continue
c
        return
        end

c************************************************c
c                                                c
c  ---- sub. No.5  of FAR2D ---- 1998.02 ----    c
c   function: set the unknowns sequence          c
c   called by: main                              c
c                                                c
c************************************************c
        subroutine unknowns (lnods,ndofn,nfixn,nsuvab,tangent)
c						   
c	  this subroutine sets the unknowns
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
c
	  dimension lnods(4,*),ndofn(*),nfixn(2,*),nsuvab(*),tangent(*)
c
c      ndofn(ipoin):       Number of Degrees of Freedom per Node
c      nfixn(idofn,ipoin): Fixity code value for the Idofn-th degree of
c                          freedom of ipoin-th node,'1' means free,
c	                     '0' restrained
c				
	  do 1 ipoin=1,npoin
	  ndofn(ipoin)=1
	  nfixn(1,ipoin)=1
	  nfixn(2,ipoin)=0
c
	  if (tangent(ipoin).le.0.0) goto 1
c
c       free in x,y direcitons for nodal points
c          outside the arc of contact
c
        ndofn(ipoin)=2
	  nfixn(2,ipoin)=1
1	  continue
c
c       restrained in y direction for nodal points at the x-axis
c
	  do 10 ielem=2,nelem_x
	  kpoin=ielem*npoin_y-nelem_y
	  ndofn(kpoin)=1
	  nfixn(2,kpoin)=0
10	  continue
c
c       restrained in x,y directions for nodes at the front
c                  and back of the workpiece
c
	  lh=npoin-npoin_y
	  do 15 ipoin=1,npoin_y
	  ndofn(ipoin)=0
	  nfixn(1,ipoin)=0
	  nfixn(2,ipoin)=0
	  ndofn(ipoin+lh)=0
	  nfixn(1,ipoin+lh)=0
	  nfixn(2,ipoin+lh)=0
15	  continue
c
c       only one unknown variable for nodes at the front
c            and back of the workpiece
c
	  ndofn(1)=1
	  nfixn(1,1)=1
	  ndofn(npoin-npoin_y+1)=1
	  nfixn(1,npoin-npoin_y+1)=1
c
c       set the tangential angle for all nodal points
c           outside of the contact arc
c
	  do 150 ipoin=1,npoin
	  if(tangent(ipoin).gt.0.0) tangent(ipoin)=0.0
150     continue
c
c       nuvab:  Number of Unknown VAriaBles
c
	  nuvab=0
	  do 300 ipoin=1,npoin
	  nuvab=nuvab+ndofn(ipoin)
300	  continue
c
	  write (*,1000)
	  write (*,1001) npoin,nuvab
1000	  format (//,
     1   5x,'OUTPUT TABLE No. 0 *** NUMBER OF UNKNOWN VARIABLES ***')
1001    format (/,
     1   5x,'NUMBER OF NODAL POINTS...... =',i8,/,
	2   5x,'NUMBER OF UNKNOWN VARIABLES. =',i8)
c
c       set up sequence of the unknown variables over all nodes
c
 	  do 50 ipoin=1,npoin_x
	  do 49 kpoin=1,npoin_y
	  node=(ipoin-1)*npoin_y+kpoin
c
	  if (ndofn(node).eq.0) goto 49
	  n1=0
	  do 40 ipoin1=1,node-1
	  n1=n1+ndofn(ipoin1)
40	  continue
c
c       nsuvab takes odd number for unknown variables in x direction
c              even for those in y direction
c
	  do 45 lpoin=1,ndofn(node)
	  nsuvab(n1+lpoin)=2*(node-1)+lpoin
45	  continue
c
49	  continue
50	  continue
c
c	  calculate the semi-width of the Hessian Matrix
c
	  ik=(lnods(3,1)-lnods(1,1)+1)*2*nuvab
c
	  write (*,1002) ik
1002    format (5x,'SEMI-WIDTH HESSIAN MATRIX... =',i8)
c
	  return
	  end
c************************************************c
c                                                c
c  ---- sub. No.6  of FAR2D ---- 1998.02 ----    c
c   function: set shape function,its derivative  c
c             for friction region                c
c   called by: initial_velocity and minfi et al  c
c                                                c
c************************************************c
        subroutine shader_f(shap_f,deriv_f,exitp)
c
c	  this subroutine sets the shape functions, Gauss points, etc.
c
	  dimension shap_f(2),deriv_f(2)
c
c	  2 node elements
c
	  t=exitp
	  shap_f(1)=(1+t)*0.5
	  shap_f(2)=(1-t)*0.5
c
	  deriv_f(1)=0.5
	  deriv_f(2)=-0.5
	  return
        end

c************************************************c
c                                                c
c  ---- sub. no.7  of far2d ---- 1998.02 ----    c
c   function: record node and ele. on surface    c
c   called by: main                              c
c                                                c
c************************************************c
       subroutine contact_elements (nefric,nofric)
c
c	 this subroutine sets the elements in contact with the roll
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
	2                 nnode,m6
c
	  dimension nefric(*),nofric(4)
c
c	  nofric: Code value of local index for nodes with friction surface
c       nefric: Code value of elements in contact with the roll
c	 	
c	  initialize nofric
c
	  do 2 inode=1,nnode
	  nofric(inode)=0
2	  continue
c
        do 1 inode=3,4
        nofric(inode)=1
1	  continue
c
c       set up contact elements  
c
        do 3 ielem=nelem_xi+1,nelem_x-nelem_xf
        melem=ielem*nelem_y
        nefric(melem)=1
3	  continue
c
	  return
        end

c************************************************c
c                                                c
c  ---- sub. no.8  of far2d ---- 1998.02 ----    c
c   function: define the nodes co-ordinates,     c
c   called by: main                              c
c                                                c
c************************************************c
        subroutine coordinates (coord,tangent,
	1                          dtl)
c
c	  this subroutine defines the co-ordinates of the nodal points
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
c
	  dimension coord(2,*),tangent(*),
     1            dtl(*)
c
c	  length of the arc of contact
c
        al=sqrt(radius**2-(radius-(thick_initial-thick_final))**2)
c
        dl=al/(nelem_x-nelem_xi-nelem_xf)
        dh0=thick_initial/nelem_y
        dh1=thick_final/nelem_y
c        
	  do 5 ipoin=2,npoin_x
	  dtl(ipoin)=dl
5	  continue
c
c	  insert one thin element to deal with 
c	    singular point at the entry of roll
c
	  if(nelem_xi.gt.2) dtl(nelem_xi+1)=0.01*dl
c
	  all=0.0
	  do 22 ipoin=2,npoin_x
        all=all+dtl(ipoin)
22	  continue
c
	  ale=0.0
	  do 25 ipoin=2,npoin_x-nelem_xf
	  ale=ale+dtl(ipoin)
25	  continue
c
c	  define the co-ordinates
c
	  do 300 ipoin=2,npoin_x
        do 300 jpoin=1,npoin_y
        kpoin=(ipoin-1)*npoin_y+jpoin
        coord(1,kpoin)=coord(1,kpoin-npoin_y)+dtl(ipoin)
300     continue
c
	  npoin_xi=nelem_xi+1
	  npoin_xf=nelem_xf+1
c
        do 301 ipoin=1,npoin_xi
	  jpoin=ipoin*npoin_y
        coord(2,jpoin)=thick_initial
301     continue
c
        do 302 ipoin=npoin_xi+1,npoin_x-npoin_xf
        jpoin=ipoin*npoin_y
        thick_y=radius+thick_final-sqrt(radius**2-
     1          (ale-coord(1,jpoin))**2)
        coord(2,jpoin)=thick_y
302     continue
c
        do 303 ipoin=npoin_x-npoin_xf+1,npoin_x
        jpoin=ipoin*npoin_y
        coord(2,jpoin)=thick_final
303     continue
c
c
	  do 320 ipoin=1,npoin_x
	  lpoin=ipoin*npoin_y
	  dhy=coord(2,lpoin)/nelem_y
	  do 320 jpoin=1,npoin_y
	  kpoin=(ipoin-1)*npoin_y+jpoin
	  coord(2,kpoin)=dhy*(jpoin-1)
320	  continue
c
c	  define the tangents
c
	  do 555 ipoin=1,npoin
	  tangent(ipoin)=1.0
555	  continue
c
        do 304 ielem=nelem_xi+1,nelem_x-nelem_xf+1
        i1=ielem*npoin_y
        tangent(i1)=(coord(1,i1)-ale)/(radius+thick_final-coord(2,i1))
	  if(abs(tangent(i1)).lt.0.0001) tangent(i1)=0.0
304     continue
c
	  return
	  end
	  
c************************************************c
c                                                c
c  ---- sub. no.9  of far2d ---- 1998.02 ----    c
c   function: set up initial flow stress         c
c   called by: main                              c
c                                                c
c************************************************c
	  subroutine flwsts (stress,shear_f)
c
c	  this subroutine establishes the initial flow stress
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
c
	  dimension stress(*),shear_f(5,*)
c	  	
        al=sqrt(radius**2-(radius-(thick_initial-thick_final))**2)
	  alfa=al/radius
c
c	  excitation values for strain and strain_rate
c
	  ef_strain=(thick_initial-thick_final)/thick_initial
	  ef_strain_rate=vel_roll*alfa/(thick_initial+thick_final)

c	  compute the stress and the interfacial shear stress

	  do 500 ielem=1,nelem
	  stress(ielem)=ga*(gb+gc*ef_strain)**gn*ef_strain_rate**gm+g0
	  do 500 igaus_f=1,ngaus_f
	  shear_f(igaus_f,ielem)=stress(ielem)*0.57735027*fric_m
500	  continue
c
	  return
        end

c************************************************c
c                                                c
c  ---- sub. no.10  of far2d ---- 1998.02 ----   c
c   function: set initial velosity field         c
c   called by: main,                             c
c                                                c
c************************************************c
        subroutine initial_velocity (lnods,nodei,nodlc,noden,
	1                               nfixn,nsuvab,nefric,nofric,jd,
     2                               coord,tangent,stress,shear_f,
     3                               a1,a2,vnode,a)
c
c     this subroutine sets the initial velocity field
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /compr/ c_factor,beta,wet,force,pc,ps   
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
	  dimension lnods(4,*),
     1            nodei(4,*),nodlc(4,*),noden(*),
	2			nfixn(2,*),nefric(*),nsuvab(*),jd(*),nofric(4),
	3	        coord(2,*),tangent(*),stress(*),shear_f(5,*),
	4			a1(*),a2(*),vnode(2,*),a(*)
	  dimension bm(3,3),bxl(3,8),ri(6),bvolum(2,4),coorde_f(2),
     1            shap_f(2),deriv_f(2),secant(5),djcb_f(5)
	  dimension bmatx(3,8),bmean(2,4),shap(4),djcb(4),
     1            deriv(2,4),cartd(2,4),coorde(2,4),gpcod(2,4)
c
c       set up constant Matrix BM
	  ik=0
	  do 5 i=1,2
	  do 6 j=1,2
	  bm(i,j)=-2.0/9.0
6	  continue
c
	  bm(i,i)=4.0/9.0
5	  continue
	  bm(3,3)=1.0/3.0
c
c       nuvab: number of unknown variables
c
	  ss=0.0
	  do 10 iuvab=1,nuvab
c
	  a1(iuvab)=0.0
c
	  do 12 jj=1,2*npoin
	  a2(jj)=0.0
12	  continue
c
c       ipoin: index of nodal points
c       ik3:   '-1 ' for unknown variables in x direction
c              '0'   for those in y direction
c
        ipoin=(nsuvab(iuvab)+1)/2
	  ik3=2*ipoin-nsuvab(iuvab)-1
c
c       common nodal points
c
	  me=1
c
c       nodal points at the front and back of the workpiece
c
	  if (iuvab.eq.1.or.iuvab.eq.nuvab) me=npoin_y
c
	  do 20 mj=1,me
c
c       nodal points at the front and back of the workpiece
c
	  if (iuvab.eq.1)  ipoin=mj
	  if (iuvab.eq.nuvab) ipoin=npoin-npoin_y+mj
c
c       knelem: number of elements containing ipoin-th node
c
	  knelem=noden(ipoin)
c
	  do 146 inelem=1,knelem
c
c       jelem: Index of inelem-th elements containing ipoin-th node
c	  inode: Local Index of ipoin-th node in inelem-th element
c
	  jelem=nodei(inelem,ipoin)
	  inode=nodlc(inelem,ipoin)
c
c       take the derivative of volume strain rate of ielem-th element 
c       as 1x1 Gauss point during calculating sub-matrix [B]
c
c       initialize bmean
	  do 37 jnode=1,nnode
        bmean(1,jnode)=0.0
        bmean(2,jnode)=0.0
37      continue
c
c       1x1 Gauss point
	  kgasp=1
        exigp=0.0
        etagp=0.0
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 77 inode1=1,nnode
	  lnode=lnods(inode1,jelem)
	  do 77 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
77	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c
c       calculate [B] Matrix at Gauss points
c
	  call bmatps(cartd,bmatx)
c
c	  calculate bmean
	  do 39 jnode=1,nnode
	  kb1=(jnode-1)*2
        bmean(1,jnode)=bmatx(1,kb1+1)+bmatx(2,kb1+1)
	  bmean(2,jnode)=bmatx(2,kb1+2)+bmatx(1,kb1+2)
39	  continue
c
c
	  kgaus=0
        do 24 i=1,2
        exigp=sgaus(i)
        do 42 j=1,2
        etagp=sgaus(j)
        kgaus=kgaus+1
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 66 inode1=1,nnode
	  lnode=lnods(inode1,jelem)
	  do 66 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
66	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgaus,shap)
c
c       calculate Determinant of Jacobian matrix sampled at gauss
c       points with the element
	  djcb(kgaus)=djaco
c
c       calculate [B] Matrix at Gauss points
c
	  call bmatps(cartd,bmatx)
c
c
c       bvolum: takes mean value of derivatives for sub-matrix [B]
c               for volume strain rate of jelem-th element
c       bxl:    special [B] matrix considering velocity boundary condition
c
	  do 81 idim=1,2
	  do 81 knode=1,nnode
	  bvolum(idim,knode)=0.0
81      continue
c
	  do 80 jl=1,3
	  do 80 il1=1,8
	  bxl(jl,il1)=0.0
80      continue
c
	  do 26 jnode=1,nnode
	  jpoin=lnods(jnode,jelem)
c
c       nodal points at the front and back of the workpiece
c
        if (jpoin.le.npoin_y) jpoin=1
        if (jpoin.gt.npoin-npoin_y) jpoin=npoin-npoin_y+1
c
	  if(2*jpoin-1.gt.nsuvab(iuvab))goto 26
c
c       calculate special [B] Matrix
c
	  do 28 ike=-1,0
c
c       only for node with one or more unknown variables
c
	  if (nfixn(1-ike,jpoin).lt.1) goto 28
c
	  kb1=(jnode-1)*2
c
c       when ike=0, sub-matrix [B] for x direciton
c       considering velocity boundary conditions on surface
	  if(ike.eq.0) then
        bx1=bmatx(1,kb1+1)
        by1=bmatx(2,kb1+2)*tangent(jpoin)
        bxy1=bmatx(3,kb1+1)+bmatx(3,kb1+2)*tangent(jpoin)
        bvolx1=bmean(1,jnode)+bmean(2,jnode)*tangent(jpoin)
c
        bvol1=bvolx1
	  else
c
c       when ik3=-1, sub-matrix [B] for y direciton
        bx1=bmatx(1,kb1+2)
        by1=bmatx(2,kb1+2)
        bxy1=bmatx(3,kb1+2)
	  bvoly1=bmean(2,jnode)
c
	  bvol1=bvoly1
        endif
c
c
   	  ilo=(jnode-1)*2+1-ike
c
c       calculate Special [B] Matrix
	  bxl(1,ilo)=bx1
	  bxl(2,ilo)=by1
	  bxl(3,ilo)=bxy1
c
c       for volume strain ratio
	  bvolum(1-ike,jnode)=bvol1
c
28	  continue
26	  continue
c
c       ri(): the matrix of [B] times [BM] for the inode-th node
c
	  kb1=2*(inode-1)+1-ik3
	  do 68 iq=1,3
	  g1=0.0		 
	  ri(iq)=0.0
	  do 70 ip=1,3
	  g1=bxl(ip,kb1)*bm(ip,iq)+g1
70	  continue
c
	  ri(iq)=g1
68	  continue
c
	  wss=2.0*(djcb(kgaus)*stress(jelem))**2
c
	  do 72 jnode=1,nnode
	  jpoin=lnods(jnode,jelem)
c
c       nodal points at the front and back of the workpiece
c
        if (jpoin.le.npoin_y) jpoin=1
        if (jpoin.gt.npoin-npoin_y) jpoin=npoin-npoin_y+1
c
	  if(2*jpoin-1.gt.nsuvab(iuvab))goto 72
c
	  do 59 ike=-1,0
	  if (nfixn(1-ike,jpoin).lt.1) goto 59
	  kb2=(jnode-1)*2+1-ike
c
c       calculate the sum of {ri} times [bxl] for jnode-th node
	  g1=0.0
	  do 73 ip=1,3
	  g1=g1+ri(ip)*bxl(ip,kb2)
73	  continue
c
c       add the item of volume strain ratio
	  g1=g1+bvolum(1-ik3,inode)*bvolum(1-ike,jnode)*c_factor
c
	  iu=2*(jpoin-1)+1-ike
c
c       calculate the left coeffience for large matrix
	  a2(iu)=a2(iu)+g1*wss/(1.0+gm)
59	  continue
72	  continue
42      continue
24	  continue
c
	  if (ik3.ne.0)goto 146
	  if (nofric(inode).ne.1.or.nefric(jelem).ne.1)goto 146
c
c		for contacted nodal points
c
	  inode_f=inode/2
	  nnode_f=2
	  lnode3=lnods(3,jelem)
	  lnode4=lnods(4,jelem)
	  coorde_f(1)=coord(1,lnode3)
	  coorde_f(2)=coord(1,lnode4)
c
	  do 85 igaus_f=1,ngaus_f
	  exitp=sgaus_f(igaus_f)
	  call shader_f(shap_f,deriv_f,exitp)
        call jacob_f(coorde_f,shap_f,deriv_f,igaus_f,
     1               djaco_f,cosi)
 	  secant(igaus_f)=1./cosi
	  djcb_f(igaus_f)=djaco_f
c
	  sigi=secant(igaus_f)
	  fss=2.0*(shear_f(igaus_f,jelem)*djcb_f(igaus_f)*
     1      wgaus_f(igaus_f))**2
c
	  do 88 jnode_f=1,nnode_f
	  node=2+jnode_f
	  jpoin=lnods(node,jelem)
        if (jpoin.le.npoin_y) jpoin=1
        if (jpoin.gt.npoin-npoin_y)jpoin=npoin-npoin_y+1
	  iu=2*jpoin-1
	  if(iu.gt.nsuvab(iuvab))goto 88
c
c       add friction item	to the left coeffiecence for large matrix
	  a2(iu)=a2(iu)+fss*shap_f(inode_f)*sigi*
     1         shap_f(jnode_f)*sigi
88	  continue
c
c	  calulate the right item of large matrix(constant items)
c
	  a1(iuvab)=a1(iuvab)+fss*vel_roll*shap_f(inode_f)*sigi
85	  continue
c
146	  continue
20	  continue
c
c	  to deal with front and back tension
c
c       back tension
	  if (iuvab.eq.1) a1(iuvab)=a1(iuvab)-b_tension/npoin_y
c
c	  front tension
	  if (iuvab.eq.nuvab)a1(iuvab)=a1(iuvab)+f_tension/npoin_y
c
c	  to store Hessian Matrix by the 1D array(a())
c
        iy=0
	  do 170 ji=1,iuvab
	  it=nsuvab(ji)
	  if (a2(it).ne.0.0) iy=1
	  if (iy.eq.0) go to 170
	  ik=ik+1
	  a(ik)=a2(it)
	  if(it.eq.nsuvab(iuvab))jd(iuvab)=ik
170     continue
10      continue
c
	  write(*,*)'ik=', ik
c
c	  to solve the large equations
c
	  call bandv(nuvab,1,ik,a,a1,jd,ir)
c
	  do 155 iuvab=1,nuvab
	  ipoin=(nsuvab(iuvab)+1)/2
	  ik3=2*ipoin-nsuvab(iuvab)-1
	  if(ik3.eq.0)  vnode(1,ipoin)=a1(iuvab)
	  if(ik3.eq.-1) vnode(2,ipoin)=a1(iuvab)
155	  continue
c
c	  call the velocity condition
c
	  call v_boundary(tangent,vnode)
c
998	  write(*,*)'primary sets is ok?   yes!' 
c	  
	  return
	  end
c************************************************c
c                                                c
c  ---- sub. no.11  of far2d ---- 1998.02 ----   c
c   function: area of contact,cosine of angle    c
c   called by: minif,initial_velocity,et al      c
c                                                c
c************************************************c
        subroutine jacob_f(coorde_f,shap_f,deriv_f,kgasp,
     1               djaco_f,cosi)

c
c	  this subroutine sets the area of the elements in contact 
c       with the roll
c
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
c
	  dimension coorde_f(2),shap_f(2),deriv_f(2),gpcod_f(5)
c
c	  coordinate of the gauss points              
c
	  gpcod_f(kgasp)=0.0
	  do 10 inode1=1,2
        gpcod_f(kgasp)=gpcod_f(kgasp)+shap_f(inode1)*coorde_f(inode1)
10	  continue
c
c	  Jacobian matrix for friction region at gauss points
c
	  xjacm_f=0.0
	  do 20 inode1=1,2
	  xjacm_f=xjacm_f+deriv_f(inode1)*coorde_f(inode1)
20	  continue
c
c       determinant of Jacobian matrix for friction region at gauss points
c
	  djaco_f=xjacm_f
c
c       the cosine value of contact angle at guass points
c
        sini=(ale-gpcod_f(kgasp))/radius
        cosi=sqrt(1.0-sini**2)
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. No.12 of FAR2D ---- 1998.02 ----    c
c   function: calculate matrix [B]               c
c   called by: initial_velocity and minfi et al  c
c                                                c
c************************************************c

	  subroutine bmatps(cartd,bmatx)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  dimension cartd(2,4),bmatx(3,8)
c
c       this subroutine calculates the strain ratio matrix B
c       for plane strain problems using the Cartesian shape
c       function derivatives
c
        kgash=0
	  do 10 inode=1,nnode
	  lgash=kgash+1
	  kgash=lgash+1
	  bmatx(1,lgash)=cartd(1,inode)
	  bmatx(1,kgash)=0.0
	  bmatx(2,lgash)=0.0
	  bmatx(2,kgash)=cartd(2,inode)
	  bmatx(3,lgash)=cartd(2,inode)
	  bmatx(3,kgash)=cartd(1,inode)
10      continue
c        
	  return
	  end
c************************************************c
c                                                c
c  ---- sub. no.13 of far2d ---- 1998.02 ----    c
c   function: compute the energy rate functional c
c   called by: main and ff                       c
c                                                c
c************************************************c
        subroutine energy(lnods,nefric,coord,stress,shear_f,
     1                    vnode,
     2                    strnrte,strt_egas,strt_qv,
     3                    strt_eqm,vegaus_x,vefric)
c
c	  this subroutine is to calculate total energy of workpiece
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /compr/ c_factor,beta,wet,force,pc,ps   
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
	  dimension lnods(4,*),
     1            nefric(*),
     2            coord(2,*),stress(*),shear_f(5,*),
	3			vnode(2,*),	
     4            strnrte(3,4,*),strt_egas(4,*),
     5            strt_qv(*),strt_eqm(*),
     6            vegaus_x(*),vefric(*)
	  dimension shap_f(2),deriv_f(2),secant(5),djcb_f(5),coorde_f(2),
     1            bmatx(3,8),cartd(2,4),deriv(2,4),coorde(2,4),
     2            djcb(4),gpcod(2,4),shap(4),bmean1(2)
c
c       strnrte(3,ngaus,ielem): Strain ratio components of ielem-th element
c                               at gauss points
c       strt_egas(ngaus,ielem): Equivalent strain ratio of ielem-th element
c
c       strt_qv(ielem):         Mean equivalent volume strain ratio of
c                               ielem-th element
c       strt_eqm(ielem):	      Mean equivalent strain ratio of 
c                               ielem-th element
c       vefric(ielem,mgaus):    relative velocity at gauss points in the
c                               surface element
c	  vegas_x(ielem,mgaus):   velocity component of element at gauss
c                               points in x direction
c
        wp=0.0
        wt=0.0
        wf=0.0
	  xs11=0.0
	  kthin=nelem_xi*nelem_y
        neg=neg+1
c
c		calculate strain rate at Gauss points
c
        do 706 ielem=1,nelem
c
c       initialize strain ratio components of element at gauss points
        do 700 igaus=1,ngaus
        do 700 n=1,3
        strnrte(n,igaus,ielem)=0.0
700	  continue
c
c       initialize mean volume strain ratio of element
        strt_qv(ielem)=0.
c
c       calculate mean volume strain ratio of element
c       taken as 1x1 gauss
        kgasp=1
        exigp=0.0
        etagp=0.0
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 77 inode1=1,nnode
	  lnode=lnods(inode1,ielem)
	  do 77 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
77	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c
c       calculate [B] Matrix at Gauss points
	  call bmatps(cartd,bmatx)
c
c       bmean1: strain ratio components of element at 1x1 gauss point
c
c       initialize bmean1
	  bmean1(1)=0.0
	  bmean1(2)=0.0
c
c       calculate bmean1
	  kb1=0
	  do 32 inode=1,nnode
        kb1=(inode-1)*2
	  ipoin=lnods(inode,ielem)
c
	  do 34 idim=1,2
	  kb1=kb1+1
	  do 34 kk=1,2
	  bmean1(kk)=bmean1(kk)+vnode(idim,ipoin)*bmatx(kk,kb1)
34      continue
c
32      continue
c
c       calculate mean volume strain ratio of element
        strt_qv(ielem)=bmean1(1)+bmean1(2)
c
c       calculate strain ratio components of element at gauss points
        kgasp=0
        do 40 i=1,2
        exigp=sgaus(i)
        do 40 j=1,2
        etagp=sgaus(j)
        kgasp=kgasp+1
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 88 inode1=1,nnode
	  lnode=lnods(inode1,ielem)
	  do 88 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
88	  continue
c
c       compute Cartesian shape function derivatives and DJACO
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c
c       calculate Determinant of Jacobian matrix sampled at gauss
c       points with the element
        djcb(kgasp)=djaco
c
c       calculate [B] Matrix at Gauss points
	  call bmatps(cartd,bmatx)
c
c       calculate strain ratio components of element at gauss points
	  kb1=0
	  do 42 inode=1,nnode
        kb1=(inode-1)*2
	  ipoin=lnods(inode,ielem)
	  do 44 idim=1,2
	  kb1=kb1+1
	  do 44 kk=1,3
	  strnrte(kk,kgasp,ielem)=strnrte(kk,kgasp,ielem)+
     1                          vnode(idim,ipoin)*bmatx(kk,kb1)
44      continue
42      continue
40      continue
c
c       initialize mean equivalent strain ratio of the elememt
        strt_eqm(ielem)=0.
        do 508 mgaus=1,ngaus
	  g1=strnrte(1,mgaus,ielem)**2
	  g2=strnrte(2,mgaus,ielem)**2
	  g3=(strnrte(1,mgaus,ielem)-strnrte(2,mgaus,ielem))**2+g1+g2
c
c       calculate equivalent strain ratio of element at gauss points
        strt_egas(mgaus,ielem)=sqrt((g3+1.5*strnrte(3,mgaus,ielem)**2)
	1                         *2./9.+strt_qv(ielem)**2*c_factor)
c
c       calculate mean equivalent strain ratio of element
	  strt_eqm(ielem)=strt_eqm(ielem)+strt_egas(mgaus,ielem)*0.25
c
c	  calculate plastic deformation energy
c
        wp=wp+stress(ielem)*strt_egas(mgaus,ielem)*
     1     djcb(mgaus)/(gm+1.0)
508     continue
c	  
c       calculate mean equivalent strain ratio of element,where
c       e11 for small deformation region
c
	  strt_eqm(ielem)=sqrt(strt_eqm(ielem)**2+e11*e11)
c
	  if(ielem.le.kthin)goto 706
	  if (strt_eqm(ielem).gt.xs11)xs11=strt_eqm(ielem)
c
706     continue
c
	  e11=xs11*0.001
c
c	  calculate friction energy 
c
	  do 709 jelem=1,nelem
c
        if(nefric(jelem).ne.1) goto 709
c
	  lnode3=lnods(3,jelem)
	  lnode4=lnods(4,jelem)
	  coorde_f(1)=coord(1,lnode3)
	  coorde_f(2)=coord(1,lnode4)
c
        do 719 mgaus_f=1,ngaus_f
	  exitp=sgaus_f(mgaus_f)
	  call shader_f(shap_f,deriv_f,exitp)
	  call jacob_f(coorde_f,shap_f,deriv_f,mgaus_f,
     1               djaco_f,cosi)
	  secant(mgaus_f)=1.0/cosi
	  djcb_f(mgaus_f)=djaco_f
c
c       calculate velocity component of element on surface
c       at gauss points in x direciton
        vegaus_x(mgaus_f)=shap_f(1)*vnode(1,lnode3)+
     1                    shap_f(2)*vnode(1,lnode4)
c
c       calculate relative velocity of elements on surface
c       at gauss points,where e2 for improving converg at
c       neutral point
        vefric(mgaus_f)=sqrt((vegaus_x(mgaus_f)*secant(mgaus_f)-
     1                  vel_roll)**2+e2)
        wf=shear_f(mgaus_f,jelem)*djcb_f(mgaus_f)*vefric(mgaus_f)*
     1     wgaus_f(mgaus_f)+wf
719     continue
c
709     continue
c
c	  calculate tension energy
c
        wt=wt+b_tension*vnode(1,1)-f_tension*vnode(1,npoin)
c
c	  calculate total energy
c        
	  w=wf+wp+wt
c
	  return
        end

c************************************************c
c                                                c
c  ---- sub. no.14 of far2d ---- 1998.02 ----    c
c   function: get  minimum value of functional   c
c   called by: main , call: a12, No.11,12 et al  c
c                                                c
c************************************************c
        subroutine minfi(lnods,nodei,nodlc,noden,
     1              nsuvab,nefric,nofric,jd,
     2              coord,tangent,stress,shear_f,
     3              a1,a2,vnode,a,
     4              strnrte,strt_egas,strt_qv,strt_eqm,strn_eq,
     5              vegaus_x,vefric,vnode0)
c
c     this subroutine is to optimize the velocity field
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /compr/ c_factor,beta,wet,force,pc,ps   
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
	  dimension lnods(4,*),
     1            nodei(4,*),nodlc(4,*),noden(*),
	2			nsuvab(*),nefric(*),jd(*),nofric(4),
     3            coord(2,*),tangent(*),
     4            stress(*),shear_f(5,*),
	5			a1(*),a2(*),vnode(2,*),a(*),
	6            strnrte(3,4,*),strt_egas(4,*),
     7            strt_qv(*),strt_eqm(*),strn_eq(*),
     8            vegaus_x(*),vefric(*),vnode0(2,*)
	  dimension bmatx(3,8),bmean(2,4),cartd(2,4),
     1            deriv(2,4),shap(4),coorde(2,4),gpcod(2,4),djcb(4)
	  dimension shap_f(2),deriv_f(2),secant(5),djcb_f(5),coorde_f(2)
c
        w1=w
        ik=0
c        m6=m6+1
c
        do 140 iuvab=1,nuvab
        a1(iuvab)=0.0
c
	  do 555 ii=1,2*npoin
	  a2(ii)=0.0
555     continue
c
c       ipoin: index of nodal points
c       ik3:   '-1 ' for unknown variables in x direction
c              '0'   for those in y direction
c
        ipoin=(nsuvab(iuvab)+1)/2
        ik3=2*ipoin-nsuvab(iuvab)-1
c
c       common nodal points
c
        me=1
c
c       nodal points at the front and back of the workpiece
c	  because there are only one unknown variable
c
        if(iuvab.eq.1.or.iuvab.eq.nuvab)  me=npoin_y
c
        do 1460 mj=1,me
c
c       nodal points at the front and back of the workpiece
c
        if(iuvab.eq.1) ipoin=mj
        if(iuvab.eq.nuvab) ipoin=npoin-npoin_y+mj
c
c       knelem: number of elements containing i1-th node
c
        knelem=noden(ipoin)
c
        do 146 inelem=1,knelem
c
c       ielem: Index of inelem-th elements containing ipoin-th node
c	  inode: Local Index of ipoin-th node in inelem-th element
c
        ielem=nodei(inelem,ipoin)
        inode=nodlc(inelem,ipoin)
c
c       take the derivative of volume strain rate of ielem-th element 
c       as 1x1 Gauss point during calculating sub-matrix [B]
c
c       initialize bmean
	  do 37 jnode=1,nnode
        bmean(1,jnode)=0.0
        bmean(2,jnode)=0.0
37      continue
c
c       1x1 Gauss point
	  kgasp=1
        exigp=0.0
        etagp=0.0
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 77 inode1=1,nnode
	  lnode=lnods(inode1,ielem)
	  do 77 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
77	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c
c       calculate [B] Matrix at Gauss points
	  call bmatps(cartd,bmatx)
c
c       calculate bmean
	  do 39 jnode=1,nnode
	  kb1=(jnode-1)*2
        bmean(1,jnode)=bmatx(1,kb1+1)+bmatx(2,kb1+1)
	  bmean(2,jnode)=bmatx(2,kb1+2)+bmatx(1,kb1+2)
39	  continue
c
c ***********************************************************
c *	  Calculate Derivative of Plastic Deformation Energy  *
c ***********************************************************
c        
	  kgasp=0
        do 148 i=1,2
        exigp=sgaus(i)
        do 42 j=1,2
        etagp=sgaus(j)
        kgasp=kgasp+1
c
c------------------------------------------------------
c       Cal. 1st derivative of plastic deformation energy
c------------------------------------------------------
c
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 88 inode1=1,nnode
	  lnode=lnods(inode1,ielem)
	  do 88 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
88	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c	  calculate Determinant of Jacobian Matrix
	  djcb(kgasp)=djaco
c
c       calculate [B] Matrix at Gauss points
	  call bmatps(cartd,bmatx)
c
c       bci: the devirative of pp with respect to velocity components in x or
c            y direction if pp=the second power of equivalent strain ratio
c
	  kb1=(inode-1)*2
c
c       when ik3=0, bci for x direciton
c       add velocity boundary condition on contacted surface
	  if(ik3.eq.0) then
        bx1=bmatx(1,kb1+1)
        by1=bmatx(2,kb1+2)*tangent(ipoin)
        bxy1=bmatx(3,kb1+1)+bmatx(3,kb1+2)*tangent(ipoin)
        bvolx1=bmean(1,inode)+bmean(2,inode)*tangent(ipoin)
	  bci1=2./9.*((2.0*bx1-by1)*strnrte(1,kgasp,ielem)+(2.0*by1-bx1)*
     1       strnrte(2,kgasp,ielem))+1./3.*bxy1*strnrte(3,kgasp,ielem)+
     2       c_factor*bvolx1*strt_qv(ielem)
        bvol1=bvolx1
	  else
c
c       when ik3=-1, bci for y direciton
        bx1=bmatx(1,kb1+2)
        by1=bmatx(2,kb1+2)
        bxy1=bmatx(3,kb1+2)
	  bvoly1=bmean(2,inode)
        bci1=2./9.*((2.0*bx1-by1)*strnrte(1,kgasp,ielem)+(2.0*by1-bx1)*
     1       strnrte(2,kgasp,ielem))+1./3.*bxy1*strnrte(3,kgasp,ielem)+
     2       c_factor*bvoly1*strt_qv(ielem)
	  bvol1=bvoly1
        endif
c
c	  cal. 1st derivative of plastic deformation energy
c
        sey=stress(ielem)*djcb(kgasp)/strt_egas(kgasp,ielem)
        a1(iuvab)=a1(iuvab)-sey*bci1
c
c-------------------------------------------------------
c      Cal. 2rd derivative of plastic deformation energy
c-------------------------------------------------
c
        do 856 jnode=1,nnode
        jpoin=lnods(jnode,ielem)
c
c       consider if jpoin belongs to the front or back of workpiece
c
        if(jpoin.le.npoin_y) jpoin=1
        if(jpoin.gt.npoin-npoin_y) jpoin=npoin-npoin_y+1
c
c       only for the unknown variables with index not large than i-th
c	  according to symmetry of Hessian Matrix
c 
        if(2*jpoin-1.gt.nsuvab(iuvab)) go to 856
c
c       set up bcd to compute sub-matrix of [B]
c          for partial derivatives of strain ratio with respect to velocity
c              components in x direction
c
        kb2=(jnode-1)*2
        bx2=bmatx(1,kb2+1)
        by2=bmatx(2,kb2+2)*tangent(jpoin)
        bxy2=bmatx(3,kb2+1)+bmatx(3,kb2+2)*tangent(jpoin)
        bvolx2=bmean(1,jnode)+bmean(2,jnode)*tangent(jpoin)
	  bci2=2./9.*((2.0*bx2-by2)*strnrte(1,kgasp,ielem)+(2.0*by2-bx2)*
     1       strnrte(2,kgasp,ielem))+1./3.*bxy2*strnrte(3,kgasp,ielem)+
     2       c_factor*bvolx2*strt_qv(ielem)
	  bvol2=bvolx2
c
        bcp=bx1*bx2+by1*by2+c_factor*bvol1*bvol2+1./3.*bxy1*bxy2
c
        z=strt_egas(kgasp,ielem)**2
c
c	Cal. one of 2rd derivative of plastic deformation energy
        m1=2*jpoin-1
        a2(m1)=a2(m1)+sey*(bcp+(gm-1.)*bci1*bci2/z)
c
c      the following is for another 2rd derivative
c
        bx2=bmatx(1,kb2+2)
        by2=bmatx(2,kb2+2)
        bxy2=bmatx(3,kb2+2)
	  bvoly2=bmean(2,jnode)
        bci2=2./9.*((2.0*bx2-by2)*strnrte(1,kgasp,ielem)+(2.0*by2-bx2)*
     1       strnrte(2,kgasp,ielem))+1./3.*bxy2*strnrte(3,kgasp,ielem)+
     2       c_factor*bvoly2*strt_qv(ielem)
	  bvol2=bvoly2
c
        ccp=bx1*bx2+by1*by2+c_factor*bvol1*bvol2+1./3.*bxy1*bxy2
c
c	  Cal. another of 2rd derivative of plastic deformation energy
        m2=2*jpoin
        a2(m2)=a2(m2)+sey*(ccp+(gm-1.)*bci1*bci2/z)
856     continue
42	  continue
148     continue
c
c ****************************************************
c *      Calculate the Derivative of Friction Energy *
c ****************************************************
c
c       only for nodes contacted with the roll and 
c          unknown variales in x direction
c
        if(ik3.ne.0) goto 146
        if(nofric(inode).ne.1.or.nefric(ielem).ne.1) goto 146
c
        inode_f=inode/2
	  nnode_f=2
	  lnode3=lnods(3,ielem)
	  lnode4=lnods(4,ielem)
	  coorde_f(1)=coord(1,lnode3)
	  coorde_f(2)=coord(1,lnode4)
	  do 149 igaus_f=1,ngaus_f
	  exitp=sgaus_f(igaus_f)
	  call shader_f(shap_f,deriv_f,exitp)
        call jacob_f(coorde_f,shap_f,deriv_f,igaus_f,
     1               djaco_f,cosi)
 	  secant(igaus_f)=1./cosi
	  djcb_f(igaus_f)=djaco_f
c
c       calculate velocity component of element on surface
c       at gauss points in x direciton
        vegaus_x(igaus_f)=shap_f(1)*vnode(1,lnode3)+
     1                    shap_f(2)*vnode(1,lnode4)
c
c       calculate relative velocity of elements on surface
c       at gauss points,where e2 for improving converg at
c       neutral point
        vefric(igaus_f)=sqrt((vegaus_x(igaus_f)*
	1                        secant(igaus_f)-vel_roll)**2+e2)
c
        sigi=secant(igaus_f)
        fss=shear_f(igaus_f,ielem)*djcb_f(igaus_f)*wgaus_f(igaus_f)*
     1      shap_f(inode_f)/vefric(igaus_f)
        pu=vegaus_x(igaus_f)*sigi-vel_roll
c--------------------------------------------------
c       calculate 1st derivative of friction energy
c--------------------------------------------------
        a1(iuvab)=a1(iuvab)-pu*fss*sigi
c
        do 147 jnode_f=1,nnode_f
        lnode=2+jnode_f
        jpoin=lnods(lnode,ielem)
c
c       nodal points at the front and back of the workpiece
c
        if(jpoin.le.npoin_y) jpoin=1
        if(jpoin.gt.npoin-npoin_y) jpoin=npoin-npoin_y+1
c
c       only for the unknown variables with index not large than i-th
c	  according to symmetry of Hessian Matrix
c
        if(2*jpoin-1.gt.nsuvab(iuvab)) go to 147
c
        sigj=secant(igaus_f)
        pdu=sigi*sigj*e2/vefric(igaus_f)**2*shap_f(jnode_f)
c------------------------------------------
c     calculate  2nd deviratives of friction energy
c-----------------------------------------
        a2(2*jpoin-1)=a2(2*jpoin-1)+pdu*fss
c
147     continue
c
149     continue
c
146     continue
c
1460    continue
c
c ***********************************************
c *     Calculate Derivatives of Tension Energy *
c ***********************************************
c
c       back tension
c
        if(ipoin.le.npoin_y) a1(iuvab)=a1(iuvab)-b_tension/npoin_y
c
c       front tension
c
        if(ipoin.gt.npoin-npoin_y)a1(iuvab)=a1(iuvab)+f_tension/npoin_y
c
c   *******************************************************
c       Storage as one-dimesion matrix for Hessian Matrix	*
c   *******************************************************
456     iy=0
        do 170 juvab=1,iuvab
        it=nsuvab(juvab)
        if(a2(it).ne.0.0) iy=1
        if(iy.eq.0) go to 170
        ik=ik+1
        a(ik)=a2(it)
c
c       set up index of element at main angal line during storage
c
        if(it.eq.nsuvab(iuvab)) jd(iuvab)=ik
170     continue
c
140     continue
c
c    *************************************** 
c    *     call resolve subroutine BANDV	 *
c    **************************************
        call bandv(nuvab,1,ik,a,a1,jd,ir)
c
c      search for the maxium value of modified velocity
        a3=0.
        do 129 iuvab=1,nuvab
        a4=abs(a1(iuvab))
129     if(a4.gt.a3) a3=a4
c
	  if(a3.lt.0.5)a3=0.5
c
	  if(up.gt.0.0003)then
c
c       choose intitial coefficience
        ek=beta/a3
        if(ek.gt.0.8) ek=0.8
	  if(ek.le.1.0e-4)ek=1.0e-4
c
c       Store old values of velocity fields
        do 128 ipoin=1,npoin
        do 128 j=1,2
128     vnode0(j,ipoin)=vnode(j,ipoin)
c
c *****************************************
c *       call Golden Section Method	    *
c *****************************************
        call a12(0.0,ek,lnods,nsuvab,nefric,coord,tangent,
     1              stress,shear_f,
     2              a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3              vegaus_x,vefric,vnode0)
c
	  else
c ***********************************
c *    call another search method   *
c ***********************************
        a3=beta/a3*1.025
	  if(a3.gt.1.0)a3=1.0
c
175     do 180 iuvab=1,nuvab
        ipoin=(nsuvab(iuvab)+1)/2
        ik3=2*ipoin-nsuvab(iuvab)-1
        if(ik3.eq.0) vnode(1,ipoin)=vnode(1,ipoin)+a3*a1(iuvab)
        if(ik3.eq.-1) vnode(2,ipoin)=vnode(2,ipoin)+a3*a1(iuvab)
180     continue
c
c      use boundary conditions of velocity field
c
        call v_boundary(tangent,vnode)
c
c      calculate total energy based on new velocity field
c
         call energy(lnods,nefric,coord,stress,shear_f,
     1              vnode,strnrte,strt_egas,strt_qv,
     2              strt_eqm,vegaus_x,vefric)
c
	  if(abs(1.0-w1/w).le.1.0e-6) goto 199
c
	  if(w.ge.w1) goto 324
	  call stress1(lnods,nefric,coord,stress,shear_f,vnode,
     1                     strt_eqm,strn_eq,vegaus_x)
	  w1=w
	  a3=a3*1.08
	  goto 175
c
324     if(abs(a3).le.1.0e-5) goto 199
	  a3=-0.618*a3
	  goto 175
	  endif
c
c       calculated relative modified value of velocity
199     uv=0.
        duv=0.
        do 258 i=1,npoin
258     uv=vnode(1,i)**2+vnode(2,i)**2+uv
        do 259 i=1,nuvab
259     duv=a1(i)**2+duv
        up=sqrt(duv/uv)
c
        wet=abs((w-w1)/w)
	  do 998 i=1,2*npoin
	  jd(i)=0
	  a1(i)=0.0
	  a2(i)=0.0
998     continue
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. no.15 of far2d ---- 1998.02 ----    c
c   function: calculate flow & friction stress   c
c   called by: main                              c
c                                                c
c************************************************c
	  subroutine stress1(lnods,nefric,coord,stress,shear_f,vnode,
     1                     strt_eqm,strn_eq,vegaus_x)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
	  dimension lnods(4,*),nefric(*),
     4            coord(2,*),
     6            stress(*),shear_f(5,*),
	9			vnode(2,*),
	1            strt_eqm(*),strn_eq(*),
     2            vegaus_x(*)
	  dimension shap_f(2),deriv_f(2),secant(5),coorde_f(2)
c
c      this subroutine calculates flow stress of element and
c      friction stress of element on surface
c
c      strn_eq:  equivalent strain of element
c
	  do 261 ielem=1,nelem_x
	  do 261 jelem=1,nelem_y
	  kelem=(ielem-1)*nelem_y+jelem
	  strn_eq(kelem)=0.0
	  stress(kelem)=0.0
	  shear_f(1,kelem)=0.0
	  shear_f(2,kelem)=0.0
c
c	  calculate mean velocity of element in x direction
	  uv=0.0
	  do 262 il=1,nnode
	  nj=lnods(il,kelem)
	  uv=uv+vnode(1,nj)*0.25
262	  continue
c
c       calculate mean length of element in x direction
	  g1=0.0
	  do 50 m=1,4,3
	  n1=lnods(m,kelem)
	  n2=n1+npoin_y
	  g1=g1+(coord(1,n2)-coord(1,n1))*0.5
50	  continue
c
c       calculate equivanent strain of element
c
	  if(ielem.gt.1)then
	  strn_eq(kelem)=strn_eq(kelem-nelem_y)+strt_eqm(kelem)*g1/uv
	  else
	  strn_eq(kelem)=strt_eqm(kelem)*g1/uv
	  endif
c
c       calculate flow stress of element
c
	  stress(kelem)=ga*(gb+gc*strn_eq(kelem))**gn*strt_eqm(kelem)**gm
     1                +g0
c
	  if(nefric(kelem).lt.1)goto 261
c
c       calculate relative velocity on surface
	  g3=stress(kelem)/sqrt(3.0)*fric_m
	  lnode3=lnods(3,kelem)
	  lnode4=lnods(4,kelem)
	  coorde_f(1)=coord(1,lnode3)
	  coorde_f(2)=coord(1,lnode4)
	  do 84 igaus_f=1,ngaus_f
	  exitp=sgaus_f(igaus_f)
	  call shader_f(shap_f,deriv_f,exitp)
	  call jacob_f(coorde_f,shap_f,deriv_f,igaus_f,
     1               djaco_f,cosi)
	  secant(igaus_f)=1./cosi
c
c       calculate velocity component of element on surface
c       at gauss points in x direciton
        vegaus_x(igaus_f)=shap_f(1)*vnode(1,lnode3)+
     1                    shap_f(2)*vnode(1,lnode4)
c
c       calculate absolute relative velocity of elements on surface
c       at gauss points
c
	  g1=abs(vegaus_x(igaus_f)*secant(igaus_f)-vel_roll)
c	  g1=abs(g1)+abs(g2)
c
c
c       calculate friction stress accoring to Kobayashi model
c
	  shear_f(igaus_f,kelem)=g3*(2./3.1415926*atan(g1*5000.))
84      continue
261	  continue
c
	  return
	  end	

c************************************************c
c                                                c
c  ---- sub. No.16 of FAR2D ---- 1990.10 ----    c
c   function: solve large symmmetrial linear eq. c
c   called by: minfi                             c
c                                                c
c************************************************c
        subroutine bandv(n,m,np,a,a1,jd,ir)
        dimension a(np),a1(n,m),jd(n)
c
c       this subroutine solve the velocity increacements
c
        do 1 i=1,n
        i0=jd(i)-i
        if(i.eq.1) go to 4
        mi=jd(i-1)-i0+1
        do 2 j=mi,i
        j0=jd(j)-j
        mj=1
        if(j.gt.1) mj=jd(j-1)-j0+1
        mij=mi
        if(mj.gt.mi) mij=mj
        ij=i0+j
        jm1=j-1
        do 3 k=mij,jm1
        if(mij.gt.jm1) go to 3
        ik=i0+k
        kk=jd(k)
        jk=j0+k
        a(ij)=a(ij)-a(ik)*a(kk)*a(jk)
3       continue
        if(j.eq.i) go to 4
        jj=jd(j)
        a(ij)=a(ij)/a(jj)
        do 2 k=1,m
2       a1(i,k)=a1(i,k)-a(ij)*a(jj)*a1(j,k)
4       ii=i0+i
        if(a(ii).eq.0.0) go to 6
        do 1 k=1,m
1       a1(i,k)=a1(i,k)/a(ii)
        do 5 l=2,n
        i=n-l+2
        i0=jd(i)-i
        mi=jd(i-1)-i0+1
        im1=i-1
        do 5 j=mi,im1
        if(mi.gt.im1) go to 5
        ij=i0+j
        do 7 k=1,m
7       a1(j,k)=a1(j,k)-a(ij)*a1(i,k)
5       continue
        ir=0
        go to 44
6       ir=i
44      continue
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. No.17 of FAR2D ---- 1990.10 ----    c
c   function: one search with golden section md. c
c   called by: minfi, call ff                    c
c                                                c
c************************************************c
        subroutine a12(x11,dx,lnods,nsuvab,nefric,coord,tangent,
     1              stress,shear_f,
     2              a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3              vegaus_x,vefric,vnode0)
	  dimension lnods(4,*),nsuvab(*),
     1            nefric(*),
     2            coord(2,*),tangent(*),stress(*),shear_f(5,*),
	3			vnode(2,*),a1(*),	
     4            strnrte(3,4,*),strt_egas(4,*),
     5            strt_qv(*),strt_eqm(*),
     6            vegaus_x(*),vefric(*),vnode0(2,*)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
c
c      this subroutine uses Golden Section Method to opitimize velocity
c
c      the first step determins search region[a,b]
        x1=x11
        x0=x1
        f0=ff(x0,lnods,nsuvab,nefric,coord,tangent,
     1          stress,shear_f,
     2          a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3          vegaus_x,vefric,vnode0)
	  f1=ff(x0+dx,lnods,nsuvab,nefric,coord,tangent,
     1             stress,shear_f,
     2             a1,vnode,strnrte,strt_egas,strt_qv,
     3             strt_eqm,vegaus_x,vefric,vnode0)
        if(f1-f0) 1,1,2
2       x2=x0+dx
        dx=-dx
1       x1=x0+dx
        f1=ff(x1,lnods,nsuvab,nefric,coord,tangent,
     1           stress,shear_f,
     2           a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3           vegaus_x,vefric,vnode0)
        if(f1-f0) 3,3,4
3       dx=2.0*dx
        x2=x0
        x0=x1
        f0=f1
        go to 1
4       b=x2
        a=x2
c
c       gets search region[a,b]
        if(x2.gt.x1) a=x1
        if(x2.lt.x1) b=x1
c
c       the second step shorts search region[a,b] further
c
c       t is constant factor
        t=(sqrt(5.0)-1.0)*0.5
c
        x1=b+t*(a-b)
        f1=ff(x1,lnods,nsuvab,nefric,coord,tangent,
     1           stress,shear_f,
     2           a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3           vegaus_x,vefric,vnode0)
        x2=a+t*(b-a)
        f2=ff(x2,lnods,nsuvab,nefric,coord,tangent,
     1           stress,shear_f,
     2           a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3           vegaus_x,vefric,vnode0)
5       if(f1-f2) 6,7,7
6       b=x2
        x2=x1
        f2=f1
        x1=b+t*(a-b)
        f1=ff(x1,lnods,nsuvab,nefric,coord,tangent,
     1           stress,shear_f,
     2           a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3           vegaus_x,vefric,vnode0)
        go to 9
7       a=x1
        x1=x2
        f1=f2
        x2=a+t*(b-a)
        f2=ff(x2,lnods,nsuvab,nefric,coord,tangent,
     1           stress,shear_f,
     2           a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3           vegaus_x,vefric,vnode0)
c
c       judgement of precision of opitimization
9       if((b-a).gt.1.025e-6) go to 5
        if((f2-f1).gt.1.025e-6) go to 5
        x0=0.5*(a+b)
        dx=ff(x0,lnods,nsuvab,nefric,coord,tangent,
     1          stress,shear_f,
     2          a1,vnode,strnrte,strt_egas,strt_qv,strt_eqm,
     3          vegaus_x,vefric,vnode0)
c
        return
        end
c************************************************c
c                                                c
c  ---- sub. no.18 of far2d ---- 1998.02 ----    c
c   function: compute functional based on new v  c
c   called by: a12,   call: ub, enegry           c
c                                                c
c************************************************c
        function ff(a3,lnods,nsuvab,nefric,coord,tangent,
     1              stress,shear_f,
     2              a1,vnode,strnrte,strt_egas,strt_qv,
     3              strt_eqm,vegaus_x,vefric,vnode0)
	  dimension lnods(4,*),nsuvab(*),
     1            nefric(*),
     2            coord(2,*),tangent(*),stress(*),shear_f(5,*),
	3			vnode(2,*),a1(*),	
     4            strnrte(3,4,*),strt_egas(4,*),
     5            strt_qv(*),strt_eqm(*),
     6            vegaus_x(*),vefric(*),vnode0(2,*)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
c
c      this fuction calculates energy based on new velocity field
c
c
c       add small increment to old velocity field
        do 180 iuvab=1,nuvab
        ipoin=(nsuvab(iuvab)+1)/2
        ik3=2*ipoin-nsuvab(iuvab)-1
        if(ik3.eq.0) vnode(1,ipoin)=vnode0(1,ipoin)+a3*a1(iuvab)
        if(ik3.eq.-1) vnode(2,ipoin)=vnode0(2,ipoin)+a3*a1(iuvab)
180     continue
c
c      use boundary conditions of velocity field
c
        call v_boundary(tangent,vnode)
c
c      calculate total energy based on new velocity field
c
         call energy(lnods,nefric,coord,stress,shear_f,
     1              vnode,strnrte,strt_egas,strt_qv,
     2              strt_eqm,vegaus_x,vefric)
        ff=w
c
        return
        end
c************************************************c
c                                                c
c  ---- sub. No.19 of FAR2D ---- 1998.02 ----    c
c   function: compute boundary condition of v.   c
c   called by: ff                                c
c                                                c
c************************************************c
        subroutine v_boundary(tangent,vnode)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  dimension tangent(*),vnode(2,*)
c
c       this subroutine computes velocity boundary conditions
c
c
c       the nodal points at the back of workpiece
        do 190 ipoin=2,npoin_y
        vnode(1,ipoin)=vnode(1,1)
	  vnode(2,ipoin)=0.0
190     continue
c
c       contacted nodal points
        do 1 ipoin=nelem_xi+1,npoin_x-nelem_xf
	  jpoin=ipoin*npoin_y
        vnode(2,jpoin)=vnode(1,jpoin)*tangent(jpoin)
1       continue
c
c       nodal points at the front of the workpiece
	  do 10 ipoin=npoin-npoin_y+1,npoin
	  vnode(1,ipoin)=vnode(1,npoin-npoin_y+1)
	  vnode(2,ipoin)=0.0
10	  continue
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. no.20 of far2d ---- 1998.02 ----    c
c   function: compute the stress field           c
c   called by: main                              c
c                                                c
c************************************************c
        subroutine strs(lnods,coord,nodei,nodlc,noden,stress,
     1                  strnrte,strt_egas,strt_qv,strt_eqm,
     2                  s2,s3)
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /compr/ c_factor,beta,wet,force,pc,ps
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
	  dimension nodei(4,*),nodlc(4,*),noden(*)
    	  dimension lnods(4,*),coord(2,*),stress(*),strnrte(3,4,*),
     1            strt_egas(4,*),strt_qv(*),strt_eqm(*),
     2            s2(4,*),s3(4,*)
        dimension c(3),cartd(2,4),deriv(2,4),shap(4),coorde(2,4),
     1            gpcod(2,4),djcb(4)
c
c       this subroutine calculates stress field and
c          non-equilibrium force of nodal points
c
	  write(*,*) '------in strs-------'
c
        z3=2./3.
c
        do 10 jelem=1,nelem
c
        do 11 n=1,4
11      s2(n,jelem)=0.0
        fep=(c_factor-2./9.)*strt_qv(jelem)
        z4=stress(jelem)/strt_eqm(jelem)
c
	  do 27 ii=1,3
	  c(ii)=0.0
27	  continue
c
        do 70 l=1,3
        do 70 kgaus=1,ngaus
        c(l)=strnrte(l,kgaus,jelem)*0.25+c(l)
70      continue
c
c       s2(n,jelem): stress components of jelem-th element
c                    n=1, normal stress in x direction
c	               n=2, normal stress in y direction
c                    n=3, shear stress
c                    n=4, volume stress
c
        s2(1,jelem)=z4*(c(1)*z3+fep)
        s2(2,jelem)=z4*(c(2)*z3+fep)
        s2(3,jelem)=z4*c(3)*z3*0.5
        s2(4,jelem)=s2(1,jelem)+s2(2,jelem)
c
c       s3(igaus,jelem): normal stress component of jelem-th element
c                        in y direction at Gauss points
c
        do 10 igaus=1,ngaus
        z5=stress(jelem)/strt_egas(igaus,jelem)
        s3(igaus,jelem)=z5*(strnrte(2,igaus,jelem)*z3+fep)
10      continue
c
c       calculate non-equilibrium force of nodal points: ps
c        
c       tx:non-equilibrium force of nodal points in x direction
c       ty:non-equilibrium force of nodal points in y direction
c       ps:total non-equilibrium force of nodal points
c
	  tx=0.0
        ty=0.0
c
        do 30 ipoin=1,npoin
        knelem=noden(ipoin)
c
        do 30 inelem=1,knelem
        jelem=nodei(inelem,ipoin)
        inode=nodlc(inelem,ipoin)
c
        kgasp=0
        do 40 i=1,2
        exigp=sgaus(i)
        do 40 j=1,2
        etagp=sgaus(j)
        kgasp=kgasp+1
c
c       calculate shape function and its derivatives
	  call shader (deriv,etagp,exigp,shap)
c
	  do 77 inode1=1,nnode
	  lnode=lnods(inode1,jelem)
	  do 77 idime=1,2
	  coorde(idime,inode1)=coord(idime,lnode)
77	  continue
c
c       compute Cartesian shape function derivatives 
        call jacob (cartd,coorde,deriv,djaco,
     1              gpcod,kgasp,shap)
c
c       calculate Determinant of Jacobian Matrix
	  djcb(kgasp)=djaco
c
c       calculate tx
	  tx=tx+(cartd(1,inode)*s2(1,jelem)+
     1     cartd(2,inode)*s2(3,jelem))*djcb(kgasp)
c
c	  calculate ty
        ty=ty+(cartd(1,inode)*s2(3,jelem)+
     1     cartd(2,inode)*s2(2,jelem))*djcb(kgasp)
c
40      continue
30      continue
c
c       calculate non-equilibrium force of nodal points: ps
        ps=sqrt(tx*tx+ty*ty)
c
	  write(*,*) 'ps=',ps
	  write(*,*) '------out strs-----'
c
	  return
        end
c************************************************c
c                                                c
c  ---- sub. no.21 of far2d ---- 1998.02 ----    c
c   function: output the computation results     c
c   called by: main , write on: fo.dat           c
c                                                c
c************************************************c
        subroutine prt(lnods,coord,stress,vnode,strt_egas,
     1                 strt_qv,strt_eqm,
     2	       	     strn_eq,s2,s3)
	  character title*70
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
	  common /compr/ c_factor,beta,wet,force,pc,ps
        common /totalcpu/cputime,cpubt
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
	  dimension lnods(4,*),coord(2,*),
     1            stress(*),vnode(2,*),
	2            strt_egas(4,*),strt_qv(*),strt_eqm(*),strn_eq(*),
	3            s2(4,*),s3(4,*)
	  dimension coorde_f(2),shap_f(2),deriv_f(2)
	  open(11,file='FEMroll2d.out')
	  n100=0
	  neg=0.0
        force=0.0
        f2=0.0
        wwf=wf/w
        wwp=wp/w
        al=sqrt(radius**2-(radius-(thick_initial-thick_final))**2)
	  do 32 ielem=nelem_xi+1,nelem_x-nelem_xf
        jelem=ielem*nelem_y
	  fric_area=0.0
	  lnode3=lnods(3,jelem)
	  lnode4=lnods(4,jelem)
	  coorde_f(1)=coord(1,lnode3)
	  coorde_f(2)=coord(1,lnode4)
	  do 34 igaus_f=1,ngaus_f
	  exitp=sgaus_f(igaus_f)
	  call shader_f(shap_f,deriv_f,exitp)
	  call jacob_f(coorde_f,shap_f,deriv_f,igaus_f,
     1               djaco_f,cosi)
	  fric_area=fric_area+djaco_f*wgaus_f(igaus_f)
34      continue
        force=force+s2(2,jelem)*fric_area
	  f2=f2+fric_area
32      continue
        pc=-force/f2
        force=-force
        s11=0.
        do 321 ielem=1,nelem
321     s11=s11+stress(ielem)
        s11=s11/nelem
        pcs=pc/s11
        hl=al/(thick_initial+thick_final)
        s12=0.
        do 302 ipoin=npoin-npoin_y+1,npoin
        s12=s12+vnode(1,ipoin)
302	  continue
        vx=s12/npoin_y/vel_roll
        sv=vx-1.
        dh=thick_initial-thick_final
        hr=dh/thick_initial*100.
	  torque=w/vel_roll*radius
        write(11,11) title
        write(11,12)
        write(11,13)thick_initial*2.,thick_final*2.,dh*2.,hr,
     1  radius,vel_roll,fric_m,e2,al,hl,radius/(thick_initial*2.),
     2  torque,ga,gc,gm,gn,g0,b_tension,f_tension,c_factor
        write(11,43)nelem_x,nelem_y,nelem,npoin
        write(11,12)
        write(11,14) vx,sv,wwf,wwp,pc
     1  ,force,f2,pcs,npoin,n100,m6,neg,up,ps
        write(11,15)w,wp,wf,wt
c
        cputime=secnds(cpubt)
        write(11,101)cputime
        write(*,101)cputime
101     format(/'Total CPU time: ',f20.2,' sec.'/)
c
        write(11,12)
        write(11,16)(jpoin,(vnode(i,jpoin),i=1,2),jpoin=1,npoin)
        write(11,12)
        write(11,17)(jelem,(s2(i,jelem),i=1,4),jelem=1,nelem)
        write(11,12)
        write(11,18)(jelem,(s3(igaus,jelem),igaus=1,ngaus),
     1               jelem=1,nelem)
        write(11,12)
        write(11,19)(ipoin,(coord(j,ipoin),j=1,2),ipoin=1,npoin)
        write(11,12)
        write(11,20)(ielem,(strt_egas(jgaus,ielem),jgaus=1,ngaus),
     1               ielem=1,nelem)
        write(11,12)
        write(11,21)(jelem,strt_qv(jelem),jelem=1,nelem)
        write(11,12)
        write(11,22)(jelem,strt_eqm(jelem),jelem=1,nelem)
        write(11,12)
        write(11,23)(jelem,strn_eq(jelem),jelem=1,nelem)
        write(11,12)
        write(11,24)(jelem,stress(jelem),jelem=1,nelem)
11      format(1x,a)
12      format(1x,3('.................'))
13      format(2x,'thick_in=',f7.2,2x,'h1=',f7.2,
     1  2x,'dh=',f7.2,2x,'dh%=',f7.3/
     2  2x,'rr=',f7.3,2x,'vel_roll=',f7.2,
     3  2x,'tm=',f7.3,2x,'e2=',f7.3/
     4  2x,'al=',f8.3,3x,'l/h=',f8.3,
     5  2x,'r/h=',f8.2,2x,'Torque=',e9.4/
     6  2x,'ga= ',f8.3,2x,'gc= ',f8.3,
     7  2x,'gm= ',f8.3,2x,'gn= ',f8.3/
     8  2x,'g0= ',f8.3,2x,'t0= ',f8.3,
     9  2x,'t1= ',f8.3,2x,'f = ',f8.0)
43      format(2x,'nelem_x= ',i3,2x,'nelem_y= ',i4,
     1  2x,'nelem= ',i4  ,2x,'npoin= ',i4)
14      format(2x,'vx= ',f8.3,2x,'sv= ',f8.3,
     1  2x,'f/w=',f8.3,2x,'p/w=',f8.3/
     2  2x,'pc= ',f8.3,2x,'force=',f8.3,
     3  2x,'f2= ',f8.1,2x,'qs= ',f8.3/
     4  2x,'npoin= ',i8,  2x,'n10=',i8,
     5  2x,'m6= ',i8,  2x,'neg=',i8//
     6  8x,'dup=',e11.4,5x,'ps=',e11.4)
15      format(/1x,'    w           wp',
     1  5x,'wf          wt',/4e12.4)
16      format(1x,2('i      u         v      '),
     1  /2(i4,2f9.2,1x))
17      format(1x,'i     cgma x  cgma y',
     1  2x,'tao xy  cgm p'/1(i4,4f8.2,1x))
18      format(/2x,'s3',4x,'1       2       3',
     1  5x,'4    '/1(i4,4f8.1,1x))
19      format(1x,2(' i     x        y      '),
     1  /2(i4,2f9.4,1x))
20      format(/2x,'strt_egas',5x,'1       2       3',
     1  5x,'4    '/1(i4,4f8.1,1x))
21      format(1x,'strt_qv'/4(i4,2x,f8.5,3x))
22      format(1x,'strt_eqm'/4(i4,2x,f8.3,3x))
23      format(1x,'strn_eq'/4(i4,2x,f8.3,3x))
24      format(1x,'s1 '/4(i4,2x,f8.3,3x))
        close(11)
c
        write(7,777)npoin_y,npoin_x,force,torque,cputime
777     format(2(1x,i5),2(1x,f10.3),1x,f18.2)
	  return
        end

c************************************************c
c                                                c
c  ---- sub. no.23 of far2d ---- 1998.02 ----    c
c   function: initial factor for Golden Method   c
c   called by: main                              c
c                                                c
c************************************************c
	  subroutine bet(err_w,err_u)
	  common /compr/ c_factor,beta,wet,force,pc,ps   
	  dimension err_w(3),err_u(3)
	  if(err_w(3).lt.err_w(1).and.err_u(3).lt.err_u(1))goto 100
	    if(err_w(3).gt.err_w(1).and.err_u(3).gt.err_u(1))then
	      if(err_w(1).gt.err_w(2).and.err_u(1).gt.err_w(2))then
	        beta=beta*0.25
	      else
	        beta=beta*0.35
	      endif
	    else
	      beta=beta*0.5
	    endif
	  goto 130
100	  continue
c
	  if(err_w(2).gt.err_w(1).and.err_u(2).gt.err_u(1))then
	     beta=beta*1.3
	   else
	   statu=+1.
	   beta=beta*0.9
	   endif
	   if(beta.gt.1.0)beta=1.0
130	   continue
         if(beta.lt.1.0e-2)beta=0.125
	   return
	   end
c************************************************c
c                                                c
c  ---- sub. no.24 of far2d ---- 1998.02 ----    c
c   function: calculate stress of nodal points   c
c   called by: main                              c
c                                                c
c************************************************c
        subroutine nodarea (lnods,coord,s2,s2_nodes,area)
c
c       this subroutine transfer the stresses from the center of
c       the element to the nodal points
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
        common /bem1/  nbem,nodlast(5),nsurf,inpt,imp,nbem1_s
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
c
	  dimension lnods(4,*),
     1            coord(2,*),s2(4,*),
	2            s2_nodes(4,*),area(*)
	  dimension cartd(2,4),coorde(2,4),deriv(2,4),gpcod(2,4),shap(4)
c
c       initialize
c
        do 10 ipoin=1,npoin
        area(ipoin)=0.0
c
        s2_nodes(1,ipoin)=0.0
        s2_nodes(2,ipoin)=0.0
        s2_nodes(3,ipoin)=0.0
        s2_nodes(4,ipoin)=0.0
10      continue
c
c       compute nodal value from adjacent elements
c
        do 20 ielem=1,nelem
        area_gaus=0.0
c
        lnode1=lnods(1,ielem)
        lnode2=lnods(2,ielem)
        lnode3=lnods(3,ielem)
        lnode4=lnods(4,ielem)
        coorde(1,1)=coord(1,lnode1)
        coorde(2,1)=coord(2,lnode1)
        coorde(1,2)=coord(1,lnode2)
        coorde(2,2)=coord(2,lnode2)
        coorde(1,3)=coord(1,lnode3)
        coorde(2,3)=coord(2,lnode3)
        coorde(1,4)=coord(1,lnode4)
        coorde(2,4)=coord(2,lnode4)
c
c       compute area using 2 x 2 gauss point rule
c
        kgasp=0
        do 40 i=1,2
        exigp=sgaus(i)
        do 40 j=1,2
        etagp=sgaus(j)
        kgasp=kgasp+1
c
        call shader (deriv,etagp,exigp,shap)
        call jacob  (cartd,coorde,deriv,djaco,
     1               gpcod,kgasp,shap)
c
        are=djaco*wgaus(i)*wgaus(j)
        area_gaus=area_gaus+are
40      continue
c
        do 20 inode=1,nnode
        lnode=lnods(inode,ielem)
        area(lnode)=area(lnode)+area_gaus
	  s2_nodes(1,lnode)=s2_nodes(1,lnode)+s2(1,ielem)*area_gaus
	  s2_nodes(2,lnode)=s2_nodes(2,lnode)+s2(2,ielem)*area_gaus
	  s2_nodes(3,lnode)=s2_nodes(3,lnode)+s2(3,ielem)*area_gaus
	  s2_nodes(4,lnode)=s2_nodes(4,lnode)+s2(4,ielem)*area_gaus
20      continue
c
c       set-up nodal value
c
        do 50 ipoin=1,npoin
	  s2_nodes(1,ipoin)=s2_nodes(1,ipoin)/area(ipoin)
	  s2_nodes(2,ipoin)=s2_nodes(2,ipoin)/area(ipoin)
	  s2_nodes(3,ipoin)=s2_nodes(3,ipoin)/area(ipoin)
	  s2_nodes(4,ipoin)=s2_nodes(4,ipoin)/area(ipoin)
50      continue
c
c       output the calculated nodal stresses along the contact arc
c
	  open(15,file='FEMpressure.dat')
	  do 66 ipoin=1,npoin_x
	  jpoin=ipoin*npoin_y
	  g1=-s2_nodes(2,jpoin)
	  g2=(ale-coord(1,jpoin))/(radius-coord(2,jpoin)+thick_final)
        g2=atan(g2)*180./3.1415926
        write(15,68)jpoin,g1,g2
66	  continue
68      format(2x,i5,2(2x,f8.3))
        close(15)
c
	  nsurf=2
c
c       external surface of the roll
c
c	  nbem1_s=10
	  nbem1_s=40
	  nbem1=nelem_x-nelem_xi-nelem_xf+4*2+nbem1_s*2
c
c       internal small circle (fixed)
c
	  nbem2=36
	  nodlast(1)=nbem1
	  nodlast(2)=nbem1+nbem2
c
c       total number of boundary elements
c
	  nbem=nodlast(2)
c
        return
        end

c*************************************************c
c                                                 c
c  ---- sub. no.25 of far2d ---- 1998.02 ----     c
c  function: set up shape function,its derivative c
c  called by: main                                c
c                                                 c
c*************************************************c
        subroutine shader (deriv,etagp,exigp,shap)
c
c       this subroutine set-up the shape functions 
c       and its derivatives  
c
c       shader supports nodarea
c
        dimension deriv(2,4),shap(4)

        s=exigp
        t=etagp
        st=s*t
c
c       shape functions (4 node element)
c 
        shap(1)=(1.-t-s+st)*0.25
        shap(2)=(1.-t+s-st)*0.25
        shap(3)=(1.+t+s+st)*0.25
        shap(4)=(1.+t-s-st)*0.25
c
c       derivatives
c
        deriv(1,1)=(-1.+t)*0.25
        deriv(1,2)=(+1.-t)*0.25
        deriv(1,3)=(+1.+t)*0.25
        deriv(1,4)=(-1.-t)*0.25
        deriv(2,1)=(-1.+s)*0.25
        deriv(2,2)=(-1.-s)*0.25
        deriv(2,3)=(+1.+s)*0.25
        deriv(2,4)=(+1.-s)*0.25
c
        return
        end
c*************************************************c
c                                                 c
c  ---- sub. no.26 of far2d ---- 1998.02 ----     c
c   function: builds jacobian matrix,cartesian    c
c             derivatives of the shape functions  c
c   called by: minfi,energy,intial_velocity et al c
c                                                 c
c*************************************************c
        subroutine jacob  (cartd,coorde,deriv,djaco,
     1                     gpcod,kgasp,shap)
c
c       this subroutine builds the jacobian matrix and the cartesian
c       derivatives of the shape functions
c
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  dimension cartd(2,4),coorde(2,4),deriv(2,4),gpcod(2,4),shap(4),
     1            xjaci(2,2),xjacm(2,2)
c
c       coordinate of the gauss points
c
        do 10 idime=1,2
        gpcod(idime,kgasp)=0.
        do 10 inode=1,nnode
        gpcod(idime,kgasp)=gpcod(idime,kgasp)+coorde(idime,inode)
     1                     *shap(inode)
10      continue
c
c       Jacobian matrix - XJACM
c       
        do 20 idime=1,2
        do 20 jdime=1,2
        xjacm(idime,jdime)=0.
        do 20 inode=1,nnode
        xjacm(idime,jdime)=xjacm(idime,jdime)+
     1                     deriv(idime,inode)*coorde(jdime,inode)
20      continue
c
c       determinant and inverse of the Jacobian matrix
c
        djaco=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
c
        xjaci(1,1)=xjacm(2,2)/djaco
        xjaci(2,2)=xjacm(1,1)/djaco
        xjaci(1,2)=-xjacm(1,2)/djaco
        xjaci(2,1)=-xjacm(2,1)/djaco
c
c       Cartesian derivatives	of the shape functions
c
        do 50 idime=1,2
        do 50 inode=1,nnode
        cartd(idime,inode)=0.0
        do 50 jdime=1,2
        cartd(idime,inode)=cartd(idime,inode)+
     1                     xjaci(idime,jdime)*deriv(jdime,inode)
50      continue
c
        return
	  end
c
c*************************************************c
c                                                 c
c  ---- sub. no.27 of far2d ---- 2000.02 ----     c
c   function: set the boundary elements and their c
c             coordinates at the surface of roll  c
c   called by: neutral2                           c
c                                                 c
c*************************************************c
c
	  subroutine prep (coord,xdie,ydie)
c
c       ---------------------------------------------------------------
c       this subroutine set the boundary elements and their
c       coordinates at the surface of the roll
c       ---------------------------------------------------------------
c
c       BEM
c       nsurf:    Number of different surface(in this case=2)
c       nodlast(i):The last node of surface i-th
c       nbem1_s:   Number of Boundary Elements for the half remainder
c                 of the external surface arc

	  character title*70
c
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /compr/ c_factor,beta,wet,force,pc,ps
        common /bem1/  nbem,nodlast(5),nsurf,inpt,imp,nbem1_s
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
        dimension coord(2,*),xdie(nbem),ydie(nbem)
c
        ale1=ale
c
	  nbem1=nodlast(1)
	  nbem2=nodlast(2)-nbem1
	  nii=nelem_x-nelem_xi-nelem_xf
        nbemi=nii+1
c
c       start from the entry of the roll in a counter clockwise direction
c
        al=sqrt(radius**2-(radius-(thick_initial-thick_final))**2)
	  alfa=asin(al/radius)
c
	  dangli=alfa/nii
	  ipoin=(nelem_xi+1)*npoin_y
	  xdie(1)=coord(1,ipoin)+radius*sin(dangli*0.05)
	  ydie(1)=radius+thick_final-sqrt(radius**2-(ale1-xdie(1))**2)
	  jpoin=ipoin+nii*npoin_y
	  xdie(nbemi)=coord(1,jpoin)-radius*sin(dangli*0.05)
	  ydie(nbemi)=radius+thick_final-sqrt(radius**2-(ale1-
     1              xdie(nbemi))**2)
c
	  do 10 ibem=2,nbemi-1
	  jpoin=(nelem_xi+ibem)*npoin_y
c
c       calculate the coordinates of boundary elements
c
	  xdie(ibem)=coord(1,jpoin)
	  ydie(ibem)=coord(2,jpoin)
10      continue
c
c       set 4 small boundary elements adjacent to the exit & entry of
c       the roll respectively
c
	  angli=0.0
	  anglj=0.0
	  do 15 i=1,4
c
c       adjacent to the exit of the roll in a counter clockwise direction
c
        angli=dangli*(i-1)**3**0.05+angli
	  ibem=nbemi+i
	  xdie(ibem)=ale1+radius*sin(angli)
	  ydie(ibem)=radius+thick_final-radius*cos(angli)
c        
c       adjacent to the entry of the roll in a clockwise direction
c
        jbem=nbem1-i+1
	  anglj=dangli*(i-1)**2**0.1+anglj
        xdie(jbem)=ale1-radius*sin(anglj+alfa)
	  ydie(jbem)=radius+thick_final-radius*cos(anglj+alfa)
15      continue
c
	  num=nbemi+4
c
c       calculate the half remainder of angle
c
        yg=thick_final+radius
        alfric1=atan((ale1-xdie(nbem1-3))/(yg-ydie(nbem1-3)))
	  alfric2=atan((xdie(nbemi+4)-ale1)/(yg-ydie(nbemi+4)))
	  alfric=alfric1+alfric2
	  angl=(2.0*3.1415926-alfric)*0.5
c
c       dii is the difference of the equal differential series,
c       supposed the angles of the boundary points takes a 
c       distribution of an equal differential series
c
	  dii=2.5
	  angl1=1.0
	  angls=angl1+(nbem1_s-1.)*dii
	  sum=nbem1_s*(angl1+angls)*0.5
	  angli=0.0
	  anglj=0.0
c
	  do 18 i=1,nbem1_s
c
	  ibem=num+i
	  angli=angli+angl/sum*(angl1+(i-1)*dii)
c
	  xdie(ibem)=ale1+radius*sin(angli+alfric2)
	  ydie(ibem)=radius+thick_final-radius*cos(angli+alfric2)
c
	  jbem=nbem1-4-(i-1)
	  anglj=angl/sum*(angl1+(i-1)*dii)+anglj
c
        xdie(jbem)=ale1-radius*sin(anglj+alfric1)
	  ydie(jbem)=radius+thick_final-radius*cos(anglj+alfric1)
c
18      continue
c
        return
	  end
c
c*************************************************c
c                                                 c
c  ---- sub. no.28 of far2d ---- 2000.02 ----     c
c   function: writes the neutral file roll.neu    c
c   called by: main                               c
c                                                 c
c*************************************************c
	  subroutine neutral2 (lnods,coord,tangent,stress,vnode,strnrte,
	1                      strt_eqm,strt_qv,strn_eq,s2,s2_nodes)
c
c       ---------------------------------------------------------------				   
c       this subroutine writes the neutral file roll.neu
c
c       FEM.neu is used by the graphical interface for 
c       post-processing (by using i_form)
c       ---------------------------------------------------------------				   
c
	  character title*70
c
        common /projt/ title
	  common /dimen/ nelem,nelem_x,nelem_y,nelem_xi,nelem_xf,
     1                 npoin,npoin_x,npoin_y,nuvab,ngaus,ngaus_f,
     2                 nnode,m6
	  common /model/ radius,vel_roll,thick_initial,thick_final,
     1                 width,f_tension,b_tension,ale,e2,e11,up,
	2                 w,wp,wf,wt,w1
	  common /mater/ ga,gb,gc,gn,gm,g0
	  common /frict/ fric_m
	  common /compr/ c_factor,beta,wet,force,pc,ps   
        common /dat1/  sgaus(2),sgaus_f(5),wgaus(2),wgaus_f(5)
        common /bem1/  nbem,nodlast(5),nsurf,inpt,imp,nbem1_s
c
	  dimension lnods(4,*),coord(2,*),tangent(*),
     1            stress(*),vnode(2,*),strnrte(3,4,*),
	2            strt_qv(*),strt_eqm(*),
     3            strn_eq(*),s2(4,*),s2_nodes(4,*),
	4            xdie(nbem),ydie(nbem)
c
	  dimension diedf(1,13)
c
	  call prep (coord,xdie,ydie)
c
c	  initialization
c
	  ineu=1
	  dummy=0.0	
        ale1=ale
c
	  n100=0
	  neg=0.0
        wwf=wf/w
        wwp=wp/w
c
c       compute torque
c
        s11=0.
        do 3 ielem=1,nelem
        s11=s11+stress(ielem)
3       continue
c
        s11=s11/nelem
        pcs=pc/s11
        al=sqrt(radius**2-(radius-(thick_initial-thick_final))**2)
        hl=al/(thick_initial+thick_final)
c
        s12=0.
        do 4 ipoin=npoin-npoin_y+1,npoin
        s12=s12+vnode(1,ipoin)
4	  continue
c
        vx=s12/npoin_y/vel_roll
        sv=vx-1.
        dh=thick_initial-thick_final
        hr=dh/thick_initial*100.
	  torque=w/vel_roll*radius	    
c
c	  *** abrir o ficheiro neutro roll.neu 
c
	  thick_initial_2=thick_initial*2.
	  thick_final_2=thick_final*2.
	  dh_2=dh*2.
	  hrdef=radius/(thick_initial*2.)
c
	  open (ineu,file='fem1.neu',form='formatted',status='unknown')
c
        write (ineu,1000) title
        write (ineu,*) npoin
c
	  do 100 ipoin=1,npoin
        vx=vel_roll/sqrt(tangent(ipoin)**2+1.)
        vy=vel_roll*tangent(ipoin)
	  write (ineu,2100) ipoin,
	1                    coord(1,ipoin),coord(2,ipoin),
     2            vnode(1,ipoin)-vx,vnode(2,ipoin)-vy,
	3                    -s2_nodes(1,ipoin),-s2_nodes(2,ipoin)
100	  continue
c
        material=1
        densy=1.0
        fract=1.0
        write(ineu,*) nelem
	  do 200 ielem=1,nelem
	  write (ineu,2200) ielem,material,
	1                    lnods(1,ielem),lnods(2,ielem),
     2                    lnods(3,ielem),lnods(4,ielem),
	3                    dummy,densy,fract
200	  continue
c
	  do 300 ielem=1,nelem
        strnrt1=0.0
        strnrt2=0.0
	  strnrt3=0.0
        do 301 igaus=1,ngaus
	  strnrt1=strnrte(1,igaus,ielem)/ngaus+strnrt1
	  strnrt2=strnrte(2,igaus,ielem)/ngaus+strnrt2
	  strnrt3=strnrte(3,igaus,ielem)/ngaus+strnrt3
301     continue
	  write (ineu,2300) ielem,
	1                    strnrt1,strnrt2,dummy,strnrt3,
	2                    strt_eqm(ielem),strt_qv(ielem)
c
c       reset the matrix 'strnrte' as strain ratio components of
c       element and in the future strain components of element
c
        strnrte(1,1,ielem)=strnrt1
        strnrte(2,1,ielem)=strnrt2
        strnrte(3,1,ielem)=strnrt3
c
        strnrte(1,2,ielem)=0.0
        strnrte(2,2,ielem)=0.0
        strnrte(3,2,ielem)=0.0
c
300	  continue
c
	  do 400 ielem=1,nelem_x
        do 400 jelem=1,nelem_y
        kelem=(ielem-1)*nelem_y+jelem
        ux=0.0
        do 402 inode=1,nnode
        ipoin=lnods(inode,kelem)
        ux=ux+vnode(1,ipoin)*0.25
402     continue
        ipoin=lnods(1,kelem)
        jpoin=ipoin+npoin_y
        dx=(coord(1,jpoin)-coord(1,ipoin))*0.5+(coord(1,jpoin+1)-
     1      coord(1,ipoin+1))*0.5
        if(ielem.le.1)then
        strnrte(1,2,kelem)=strnrte(1,1,kelem)*dx/ux
        strnrte(2,2,kelem)=strnrte(2,1,kelem)*dx/ux
        strnrte(3,2,kelem)=strnrte(3,1,kelem)*dx/ux
        else
        strnrte(1,2,kelem)=strnrte(1,2,kelem-nelem_y)+
     1                     strnrte(1,1,kelem)*dx/ux
        strnrte(2,2,kelem)=strnrte(2,2,kelem-nelem_y)+
     1                     strnrte(2,1,kelem)*dx/ux
        strnrte(3,2,kelem)=strnrte(3,2,kelem-nelem_y)+
     1                     strnrte(3,1,kelem)*dx/ux
        endif
c
        xgash=(strnrte(1,2,kelem)+strnrte(2,2,kelem))*0.5
        xgish=(strnrte(1,2,kelem)-strnrte(2,2,kelem))*0.5
        xgesh=strnrte(3,2,kelem)
        xgosh=sqrt(xgish*xgish+xgesh*xgesh)
        str1=xgash+xgosh
        str3=xgash-xgosh
        if(xgish.eq.0.)xgish=1.0e-8
        ang13=atan(xgesh/xgish)*28.64788976
c
	  write (ineu,2400) kelem,
     1 		            strnrte(1,2,kelem),strnrte(2,2,kelem),
     2                    dummy,strnrte(3,2,kelem),
     3                    strn_eq(kelem),str1,str3,ang13
400	  continue
c
	  do 500 ielem=1,nelem
	  sz=0.5*(s2(1,ielem)+s2(2,ielem))
	  sm=(1./3.)*(s2(1,ielem)+s2(2,ielem)+sz)
c
	  write (ineu,2500) ielem,
     1 		            s2(1,ielem),s2(2,ielem),
     2                    sz,s2(3,ielem),
     3                    stress(ielem),sm
500	  continue
c
        dtemp=0.0
        temp=0.0
	  do 600 ipoin=1,npoin
        write (ineu,2600) jpoin,dtemp,temp
600	  continue
c
c       geometria do rolo
c
        nbem1=nodlast(1)
	  ndies=1
	  mdpoi=nbem1
	  write (ineu,3000) ndies
	  write (ineu,3050) mdpoi
c
	  ndpoi=mdpoi
c
c	  centre of the roll
c
	  diedf(1,1)=ale1
	  diedf(1,2)=thick_final+radius
c
c	  horizontal velocity of the roll
c
	  diedf(1,3)=vel_roll
	  diedf(1,4)=0.0
c
c	  friction factor
c
	  diedf(1,5)=fric_m
c
c	  forces
c
	  diedf(1,6)=0.0
	  diedf(1,7)=0.0
c
        tdie=1.0
        xlength=1.25*ale
        thickmax=radius-sqrt(radius**2-xlength**2)+thick_final
        ntry=0
        do 711 ipto=1,ndpoi
	  x1=xdie(ipto)
	  y1=ydie(ipto)
        if(y1.le.thickmax)goto 711
	  y1=thickmax+0.001
        ydie(ipto)=y1
        ntry=ntry+1
711     continue
c
        ntry1=ntry
        xx=xlength-ale-0.001
        do 713 ipto=1,ndpoi
        y1=ydie(ipto)
        if(y1.le.thickmax)goto 713
        ntry1=ntry1-1
        xdie(ipto)=2.*xx*ntry1/ntry-xx
713     continue
c
	  do 700 idies=1,ndies
	  write (ineu,3100) idies,ndpoi,tdie
	  write (ineu,3200) (diedf(idies,ipos1),ipos1=1,7)
	  npto=ndpoi
c
	  do 800 ipto=1,npto
	  x1=xdie(ipto)
	  y1=ydie(ipto)
	  write (ineu,3300) x1,y1
800	  continue
	  write (ineu,3300) xdie(1),ydie(1)
700	  continue
c
	  close(ineu)
c
c	  output formats 
c
1000    format (1x,a)
1050    format (2x,i5,2x,i5,2x,i5,2x,i5)
1100    format (2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7)
1150    format (2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7)
1200    format (2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7)
1250    format (2x,e18.7,2x,e18.7)
1300    format (2x,e18.7,2x,e18.7)
1350    format (2x,e18.7,2x,e18.7,2x,e18.7,2x,e18.7)
1400    format (2x,e18.7,2x,e18.7,2x,e18.7)
1450    format (2x,e18.7,2x,e18.7,2x,e18.7)
2100	  format (1x,i5,6(1x,e18.7))
2200	  format (1x,i5,5(1x,i5),3(1x,e18.7))
2300	  format (1x,i5,6(1x,e18.7))
2400    format (1x,i5,8(1x,e14.7))
2500	  format (1x,i5,6(1x,e18.7))
2600    format (1x,i5,2x,d18.7,2x,d18.7)
3000	  format (1x,i5)
3050	  format (1x,i5)
3100	  format (1x,i5,4x,i5,1x,d18.7)
3200	  format (1x,d18.7,2x,d18.7,2x,d18.7,2x,d18.7,2x,f8.4,2x,
     1	         d18.7,2x,d18.7)
3300	  format(1x,d18.7,2x,d18.7)
c
        inpt=15
        npbe=(nelem_x+nelem_y)*2+2
        open(inpt,file='fem.dat',form='formatted',status='unknown')
        write(inpt,'(a70)')title
        write(inpt,*)0,1,dummy,dummy
        write(inpt,*)dummy,dummy,dummy,1.0,dummy
        write(inpt,*)dummy,1,dummy,dummy,1,dummy,1
        write(inpt,*)nnode,2,0,1,1,1
        write(inpt,*)1.0
        write(inpt,*)0,0,0
        write(inpt,*)npoin,nelem
        write(inpt,*)npbe,0
        write(inpt,*)(jpoin,(coord(idime,jpoin),idime=1,2),
     1                                         jpoin=1,npoin)
        do 4322 ielem=1,nelem
        write(inpt,*)ielem,material,
     1                  lnods(1,ielem),lnods(2,ielem),
	2                  lnods(3,ielem),lnods(4,ielem)
4322    continue
c
        close(inpt)
c
        open(inpt,file='die.dat',form='formatted',status='unknown')
c
	  write (inpt,3000) ndies
	  write (inpt,3050) mdpoi
c
        tdie=1.0
	  do 702 idies=1,ndies
	  write (inpt,3100) idies,ndpoi,tdie
	  write (inpt,3200) (diedf(idies,ipos1),ipos1=1,7)
	  npto=ndpoi
c
	  do 802 ipto=1,npto
	  x1=xdie(ipto)
	  y1=ydie(ipto)
	  write (inpt,3300) x1,y1
802	  continue
	  write (inpt,3300) xdie(1),ydie(1)
702	  continue
c
        close(inpt)
c
c       output the calculated nodal stresses along the contact arc
c
	  open(15,file='FEMpressure.dat')
	  do 66 ipoin=npoin_x,nelem_xi+1,-1
	  jpoin=ipoin*npoin_y
	  g1=-s2_nodes(2,jpoin)
	  g2=(ale-coord(1,jpoin))/(radius-coord(2,jpoin)+thick_final)
        g2=atan(g2)*180./3.1415926
        write(15,68)jpoin,g1,g2
66	  continue
68      format(2x,i5,2(2x,f8.3))
        close(15)
c
	  return
	  end


