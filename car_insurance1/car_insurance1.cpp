/*
 * car_insurance1.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: Shangwu Xiong
 */
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "MyTroubles.h"
#include "car_insurance1.h"

void car_insurance1::SaveRecords(customer_data *entry1){
	  std::string lines;
	  int number_of_lines = 0;

	  std::fstream myfile;
	myfile.open(inputfile_name,std::ios::in|std::ios::out);
	if(!myfile.is_open()){
		std::cout<<" customers.dat cannot be open! "<<std::endl;
		std::cout<<" Create a new one automatically!!"<<" "<<std::endl;
        myfile.close();
        count=0;
	}
	else
	{ while (std::getline(myfile, lines))  {  ++number_of_lines;		};
	    myfile.close();
	    count=(number_of_lines-1)/6;
	  };
        std::ofstream myfile2;
	    myfile2.open(inputfile_name,std::ios::out|std::ios::app);
      if(myfile2.is_open()){count=count+1; if(count<=1)myfile2<<count<<std::endl;
        myfile2<<count<<"  "<<entry1->FirstName<<" "<<entry1->LastName<<std::endl;
        myfile2<<entry1->DateofBirth.month<<" "<<entry1->DateofBirth.date<<" "<<entry1->DateofBirth.year;
        myfile2<<std::endl;
        myfile2<<entry1->DateOfRecord.month<<" "<<entry1->DateOfRecord.date<<" "<<entry1->DateOfRecord.year<<" , ";
        myfile2<<entry1->DateOfOpening.month<<" "<<entry1->DateOfOpening.date<<" "<<entry1->DateOfOpening.year;
        myfile2<<std::endl;
        myfile2<<entry1->DateOfStart.month<<" "<<entry1->DateOfStart.date<<" "<<entry1->DateOfStart.year<<" , ";
        myfile2<<entry1->DateOfExpiry.month<<" "<<entry1->DateOfExpiry.date<<" "<<entry1->DateOfExpiry.year;
        myfile2<<std::endl;
        myfile2<<entry1->RenewalCost<<" "<<carValue<<" "<<accidents<<" "<<riskSurcharge<<" "<<record_status<<std::endl;
        myfile2<<entry1->car_model.Model<<" "<<entry1->car_model.Manufacturer<<" "<<entry1->car_model.VIN<<" "<<entry1->car_model.Year;
        myfile2<<std::endl;
        myfile2.close();
       };

	head= new customer_data;
	tail= new customer_data;
	head->next = NULL;
	head->prev = NULL;
};

car_insurance1::car_insurance1(int count){
	car_insurance1::inputfile_name="customers.dat";
	std::cout<<"*****constructor is called*****"<<std::endl;
	car_insurance1::carValue=1000.;car_insurance1::accidents=0;car_insurance1::riskSurcharge='0';
	car_insurance1::premium=0.0; car_insurance1::count=0;car_insurance1::record_status='N';
	head= new customer_data;
	head->prev=NULL;
	head->next=NULL;
	tail=head;
};

car_insurance1::car_insurance1(double carValue,double accidents,char riskSurcharge,customer_data *entry1){
	int count=0;
	double premium=0.0;

	head= new customer_data;
    head=entry1;
	head->next = NULL;
	head->prev = NULL;
	tail= head;
	car_insurance1::inputfile_name="customers.dat";
   if(carValue<=0.)	carValue=0.0;
   if(accidents<=0.)accidents=0.0;
   if(riskSurcharge<=0.)riskSurcharge=0.0;

//   std::cout<< "Customer Name:"<<entry1->FirstName<<" "<<entry1->LastName<<std::endl;

//   std::cout << "the value of the car to be insured:"<<carValue<<"$\n"<<std::endl;

	car_insurance1::carValue=carValue;
	car_insurance1::accidents=accidents;
	car_insurance1::riskSurcharge=riskSurcharge;
	car_insurance1::premium=premium;
    car_insurance1::count=count;
    car_insurance1::record_status='N';
};


void car_insurance1::insurance_Base()
{
	car_insurance1::premium = car_insurance1::carValue /base0;
};
void car_insurance1::insurance_RiskSurcharge()
{
	   if(riskSurcharge=='A'||riskSurcharge=='a'){premium=premium+A_a;}
	   else if(riskSurcharge=='B'||riskSurcharge=='b'){premium=premium+B_b;}
	   else if(riskSurcharge=='C'||riskSurcharge=='c'){premium=premium+C_c;}
	   else if(riskSurcharge=='D'||riskSurcharge=='d'){premium=premium+D_d;}
	   else if(riskSurcharge=='E'||riskSurcharge=='E'){premium=premium+E_e;}
	   else {premium=premium+F_f;};
};

void car_insurance1::insurance_accidentsAdjust()
{ Cannot_insured cannot_insured;
	try{
	   if (accidents<=1){premium=premium*(accidents*accident_0)+premium;}
	   else if (accidents==2)
	  	   {premium=premium*(accidents*accident_1)+premium;}
	   else if (accidents<=4)
		   {premium=premium*(accidents*accident_2)+premium;}
	   else if(accidents==5){premium=premium*(accidents*accident_3)+premium;};

	   if (accidents>5)	throw cannot_insured;
     }
     catch (Cannot_insured &e)
     {
    	std::cout<<"Exception:"<<e.what()<<std::endl;
    	std::cout<<"Sorry! Program will exit automatically, bye!"<<std::endl;
        exit(1);
     };
	};

void car_insurance1::printrecord(customer_data *entry1)const{
	std::cout.setf(std::ios::fixed);
	std::cout.setf(std::ios::showpoint);
	std::cout.precision(2);
	   std::cout<<std::endl<< "*** Summary of Car Insurance ***"<<std::endl;
	   std::cout<< "Customer Name:"<<entry1->FirstName<<" "<<entry1->LastName;
       std::cout<<std::endl<<"Insurance will start at:"<<entry1->DateOfStart.month<<"/"<<entry1->DateOfStart.date<<"/"<<entry1->DateOfStart.year;
       std::cout<<std::endl<<"Insurance will expire by:"<<entry1->DateOfExpiry.month<<"/"<<entry1->DateOfExpiry.date<<"/"<<entry1->DateOfExpiry.year;
       std::cout<<std::endl<<"Insured Car:"<<entry1->car_model.Model<<" "<<entry1->car_model.Manufacturer<<","<<"VIN:"<<entry1->car_model.VIN;
	   std::cout<<std::endl<< "the value of the car to be insured:"<<carValue<<"$\n";

    std::cout<<std::endl<< "The premium for this policy holder will be $" << premium <<std:: endl <<std:: endl;
	std::cout<<std::endl<< "Save a huge! Please dial 800-800-800." << std:: endl;
};

void car_insurance1::ReviseRecord (int &choice,customer_data *entry1){
	  std::string lines;
	  int number_of_lines = 0, Count;

	  std::fstream myfile;
	myfile.open(inputfile_name,std::ios::in|std::ios::out);
	if(!myfile.is_open()){
		std::cout<<" customers.dat cannot be open! "<<std::endl;
		std::cout<<"Option is changed to 1: Add/Renew A Insurance Policy"<<" "<<std::endl;
      myfile.close();
      count=0;
      choice=1; return;
	}
	else
	{ while (std::getline(myfile, lines))  {  ++number_of_lines;		};
	    myfile.close();
	    count=(number_of_lines-1)/6;
        if(count<1){choice=1;
                   std::cout<<"Option is changed to 1: Add/Renew A Insurance Policy"<<" "<<std::endl;
                   return;};
        number_of_lines = 0;
        myfile.open(inputfile_name,std::ios::in);
        std::getline(myfile, lines);
        customer_data *tmp;
        for(Count=0; Count<count;Count++){
            tmp = tail;
            tail -> next = new customer_data;
            tail = tail -> next;
            tail -> next = NULL;


        	myfile>>lines>>tail->FirstName>>tail->LastName;
        	myfile>>tail->DateofBirth.month>>tail->DateofBirth.date>>tail->DateofBirth.year;
        	myfile>>tail->DateOfRecord.month>>tail->DateOfRecord.date>>tail->DateOfRecord.year>>lines;
        	myfile>>tail->DateOfOpening.month>>tail->DateOfOpening.date>>tail->DateOfOpening.year;
        	myfile>>tail->DateOfStart.month>>tail->DateOfStart.date>>tail->DateOfStart.year>>lines;
        	myfile>>tail->DateOfExpiry.month>>tail->DateOfExpiry.date>>tail->DateOfExpiry.year;
            myfile>>tail->RenewalCost>>tail->CarValue>>tail->Accidents>>tail->RiskSurcharge>>tail->Record_status;
            myfile>>tail->car_model.Model>>tail->car_model.Manufacturer>>tail->car_model.VIN>>tail->car_model.Year;


        if((entry1->FirstName==tail->FirstName)&&(entry1->LastName==tail->LastName))
           {std::cout<<"Customer's name is found: "<<tail->FirstName<<" "<<tail->LastName <<std::endl;
           number_of_lines =1;
           if(choice==2){tail->Record_status='V';}
           else if(choice==3){tail->Record_status='R';
 // the following policy content are changed
//           (omitted .... )
             };
           };
            tail -> prev = tmp;
        };
        myfile.close();

        if(number_of_lines ==0){std::cout<<"   Customer Name is not found!"<<std::endl;}
        else { myfile.open(inputfile_name,std::ios::out);
               myfile<<Count<<std::endl; Count =0;
        	for(tmp = head -> next; tmp!= NULL; tmp = tmp -> next) {
                 Count=Count+1;
             	myfile<<Count<<" "<<tmp->FirstName<<" "<<tmp->LastName<<std::endl;
             	myfile<<tmp->DateofBirth.month<<" "<<tmp->DateofBirth.date<<" "<<tmp->DateofBirth.year;
             	myfile<<std::endl;
             	myfile<<tmp->DateOfRecord.month<<" "<<tmp->DateOfRecord.date<<" "<<tmp->DateOfRecord.year<<" , ";
             	myfile<<tmp->DateOfOpening.month<<" "<<tmp->DateOfOpening.date<<" "<<tmp->DateOfOpening.year;
             	myfile<<std::endl;
             	myfile<<tmp->DateOfStart.month<<" "<<tmp->DateOfStart.date<<" "<<tmp->DateOfStart.year<<" , ";
             	myfile<<tmp->DateOfExpiry.month<<" "<<tmp->DateOfExpiry.date<<" "<<tmp->DateOfExpiry.year;
                myfile<<std::endl;
                myfile<<tmp->RenewalCost<<" "<<tmp->CarValue<<tmp->Accidents<<" "<<tmp->RiskSurcharge<<" "<<tmp->Record_status<<std::endl;
                myfile<<tmp->car_model.Model<<" "<<tmp->car_model.Manufacturer<<" "<<tmp->car_model.VIN<<" "<<tmp->car_model.Year;
                myfile<<std::endl;
                   };
        	myfile.close();
              };


	  };

};
void car_insurance1::VoidRecord(){

};

car_insurance1::~car_insurance1(){};
