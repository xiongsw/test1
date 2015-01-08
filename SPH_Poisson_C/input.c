/*
 	 	 * input.c
 *
 *  Created on: Oct 9, 2013
 *      Author: Shangwu Xiong
 */
#include <stdio.h>
#include "user_defined.h"
#include <stdlib.h>

void input1(int idat, int ii9)
	{
/*c********************************************************c
c*/
       FILE *input;
		char para_name[256], para_value[256],value1[80],value2[80],value3[80],value4[80],value5[80];
		int num_of_lines=icase=0;
	    int c=0,lc=0;
		double v[100];
//
		printf("Please input from screen.......\n");
		printf(" 1: Poisson problem\n");
		printf(" 2: steady-state journal bearing\n\n");
		if(scanf("%d", &icase) > 0)printf("%d\n", icase);
		if(icase<=1)icase=1;
		if(icase>=2)icase=2;
		nmethod=icase;
//
		if(icase<=1)
		{ if((input=fopen("Poisson.dat","r"))==NULL)
			   {
                 fprintf(stderr, "can not open the input file : Poisson.dat");
                 exit(1);
                }


//     while((c=fgetc(input)) != EOF)
//     {
//        if( c == '\n') num_of_lines ++;
//        lc = c;
//     }
//      if(lc != '\n') num_of_lines++;
//      rewind(input);

  // for(lc=1;lc <= num_of_lines; lc++)
  //  {
  //      fscanf(input, "%s %s", para_name,para_value);

   //     v[lc] = atof(para_value);
    //}
	/**/

//  1st line data
       c=fscanf(input, "%s %s %s %s %s", value1,value2,value3,value4,value5);

       if (c == EOF) {
           if (ferror(input)) {
               perror("fscanf");
           }
           else {
               fprintf(stderr, "Error: fscanf reached the end of file, there are no matching characters\n");
           }
           exit(1);
       }
       else if (c != 5) {
           fprintf(stderr, "Error: fscanf failed to match and only assigned %i input items, 5 expected\n", c);
           exit (1);
       }
	   xmin=atof(value1);
	   xmax=atof(value2);
	   ymin=atof(value3);
	   ymax=atof(value4);
	   alfa=atof(value5);

//  2nd line data
	   lc=fscanf(input, "%s %s",value1,value2);
       if (lc == EOF) {
           if (ferror(input)) {
               perror("fscanf");
           }
           else {
               fprintf(stderr, "Error: fscanf reached the end of file, there are no matching characters\n");
           }
           exit(1);
       }
       else if (lc!= 2) {
           fprintf(stderr, "Error: fscanf failed to match and only assigned %i input items, 2 expected\n", lc);
           exit (1);
       }

	   npoin_x=atoi(value1);
	   npoin_y=atoi(value2);
	   fclose(input);
       if(alfa>1.2)alfa=1.2;
       if(alfa<0.6)alfa=0.6;
       if(npoin_x<=10)npoin_x=10;
       if(npoin_y<=10)npoin_y=10;
       icase=1;
       imod=2;
	   return;}
//
	 else
	 {    if((input=fopen("Bearing.dat","r"))==NULL)
			   {
              fprintf(stderr, "can not open the input file : Bearing.dat");
              exit(1);
             }
// 1st line data
	   c=fscanf(input, "%s %s %s %s %s", value1,value2,value3,value4,value5);
       if (c == EOF) {
           if (ferror(input)) {
               perror("fscanf");
           }
           else {
               fprintf(stderr, "Error: fscanf reached the end of file, there are no matching characters\n");
           }
           exit(1);
       }
       else if (c != 5) {
           fprintf(stderr, "Error: fscanf failed to match and only assigned %i input items, 5 expected\n", c);
           exit (1);
       }

	   rshaft=atof(value1);
	   rshaft1=atof(value2);
	   hlength=atof(value3);
	   clrnce=atof(value4);
	   alfa=atof(value5);


//  2nd line data
	   lc=fscanf(input, "%s %s",value1,value2);
       if (lc == EOF) {
           if (ferror(input)) {
               perror("fscanf");
           }
           else {
               fprintf(stderr, "Error: fscanf reached the end of file, there are no matching characters\n");
           }
           exit(1);
       }
       else if (lc!= 2) {
           fprintf(stderr, "Error: fscanf failed to match and only assigned %i input items, 2 expected\n", lc);
           exit (1);
       }

	   npoin_x=atoi(value1);
	   npoin_y=atoi(value2);

       printf("npoin_x= %d, npoin_y=%d\n",npoin_x,npoin_y );
//
// 3rd line data
	   c=fscanf(input, "%s %s %s %s",value1,value2,value3,value4);
       if (c == EOF) {
           if (ferror(input)) {
               perror("fscanf");
           }
           else {
               fprintf(stderr, "Error: fscanf reached the end of file, there are no matching characters\n");
           }
           exit(1);
       }
       else if (c != 4) {
           fprintf(stderr, "Error: fscanf failed to match and only assigned %i input items, 4 expected\n", c);
           exit (1);
       }

	   omega1=atof(value1);
	   omega2=atof(value2);
	   eta0=atof(value3);
	   sigma=atof(value4);

//c
       xmin=0.0;
       xmax=rshaft*2.*pi;
       ymin=0.0;
       ymax= hlength;
       eccen1[1]=0.0;
       eccen1[2]=-clrnce*0.9;
       eta0=eta0*1.0e-9;
       omega2=0.0;

        if(alfa>1.2)alfa=1.2;
        if(alfa<0.6)alfa=0.6;
        if(npoin_x<=10)npoin_x=10;
        if(npoin_y<=10)npoin_y=10;

       fclose(input);
//c
        icase=1;
        imod=2;
//        return 0;
	 }
}
//*c


