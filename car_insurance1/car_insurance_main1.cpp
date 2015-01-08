/*
 * car_insurance1.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: Shangwu Xiong
 */

// This program is used to calculate 6-month policy for vehicle insurance
// a particular car's value, the number of accidents and residence sub-charges.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include "car_insurance1.h"
using namespace std;

void New_customer_regist(double &carValue, double &accidents, char &riskSurcharge,customer_data &entry);
int check_leapYear(int year);
void check_calendar(int &month,int &date,int &year);

int main()
{   customer_data entry;
	double carValue=0, premium=0;
	double accidents=0;
	char riskSurcharge='A';string LastName, FirstName;
    int choice=0;
	// screen output environment using Eclipse IDE
	    setvbuf(stdout, NULL, _IONBF, 0 );
	    setvbuf(stderr, NULL, _IONBF, 0);

	cout << "Welcome to Use Insurance Software 1.0" << endl;
	cout << "Automotive Policy 6-Month premium Calculator" << endl << endl;
//
    cout<<".:::::SELECT MENU OPTION:::::.\n\n";
       cout<<"1. Add/Renew   A Insurance Policy.\n";
       cout<<"2. Delete/Void A Insurance Policy.\n";
       cout<<"3. Revise/Edit A Insurance Policy.\n";
       cout<<"4. Exit The Program.\n\n\n";
       cout<<"SELECT#";
       cin>>choice;

//
    if(choice==4){cout<<"\nThanks, Bye!";return 0;};
//
    if(choice==2||choice==3){
    	cout<<"Input customer first name insurance will be changed"<<endl;
    	cin>>entry.FirstName;
    	cout<<"Input customer last name insurance will be changed"<<endl;
    	cin>>entry.LastName;

 	     customer_data *entry1=&entry;
    	 car_insurance1 *car_insurance=new car_insurance1(0); // just for fun
    	 car_insurance->ReviseRecord (choice,entry1);
         if(choice==1)cout<<"Choice=1 and will register a new customer!"<<endl;
         delete car_insurance;
        };
//
    if(choice==1){
	              New_customer_regist(carValue, accidents, riskSurcharge,entry);
 	    customer_data *entry1=&entry;
   	    car_insurance1 *car_insurance=new car_insurance1(carValue, accidents,riskSurcharge,entry1);
        car_insurance->insurance_Base();
        car_insurance->insurance_RiskSurcharge();
        car_insurance->insurance_accidentsAdjust();
        premium=car_insurance->getpremium();
        entry.RenewalCost=premium;
         car_insurance->SaveRecords(entry1);
         car_insurance->printrecord(entry1);
 	     delete car_insurance;
        };

	return 0;
};

void New_customer_regist(double &carValue, double &accidents, char &riskSurcharge,customer_data &entry)
		{int Date,Month,Year;char sure='N';
    time_t now = std::time(NULL);
    tm* timePtr = localtime(&now);
    int current_year=(timePtr->tm_year+1900), current_month=timePtr->tm_mon+1, today=timePtr->tm_mday;
    std::cout<<std::endl<<"Today time= ";
    std::cout<<current_month<<"/"<<today<<"/"<<current_year<<std::endl;
    entry.DateOfRecord.date=today;
    entry.DateOfRecord.month=current_month;
    entry.DateOfRecord.year=current_year;
    entry.DateOfOpening.date=today;
    entry.DateOfOpening.month=current_month;
    entry.DateOfOpening.year=current_year;

	       cout<<"Enter customer First Name:"<<endl;
	       cin>>entry.FirstName;
	       cout<<"Enter customer Last Name:"<<endl;
	       cin>>entry.LastName;
	       cout<<"Enter Birthday of customer (month, date,  year):"<<endl;
	       cin>>Month>>Date>>Year;

	       check_calendar(Month,Date,Year);
           entry.DateofBirth.month=Month;
           entry.DateofBirth.date=Date;
           entry.DateofBirth.year=Year;
	       if(current_year-Year<16){std::cout<<"Sorry, you are younger than 16 "; exit(1);}
	       if(current_year-Year>100){std::cout<<"Sorry, you are over 100 "; exit(1);}

       while(sure=='N'||sure=='n')
        {  cout<<"Enter date of Insurance will start (month, date,  year):"<<endl;
           std::cin.clear();
	       cin>>Month>>Date>>Year;
	       check_calendar(Month,Date,Year);
           while((Year<current_year)||(Year>current_year+1)){cout<<"Year should be "<<current_year<<" or "<<current_year+1<<endl;
                                                            cin>>Year;};
           if(Year==current_year){if((Month<current_month)||(Month==current_month&&Date<today))
                                    {cout<<"Insurance has to start from today!";
                                     Month=current_month;Date=today;};
                                };
           cout<<"     Start at "<<Month<<"/"<<Date<<"/"<<Year<<endl<<"    Are you sure?? (Y or y: Yes; N or n: No)"<<endl;
           cin>>sure;};
	       entry.DateOfStart.month=Month;
	       entry.DateOfStart.date=Date;
	       entry.DateOfStart.year=Year;

//	       cout<<"Enter insurance Expire date of customer (month, date,  year):"<<endl;
//	       cin>>Month>>Date>>Year;
//	       check_calendar(Month,Date,Year);
// fixed 6 months
	       Month=Month+6;
	       if(Month>12){Year=Year+1;Month=Month-12;};
	       entry.DateOfExpiry.month=Month;
	       entry.DateOfExpiry.date=Date;
	       entry.DateOfExpiry.year=Year;

//
            cout<<"Enter Car Model: (e.g. Chevrolet)"<<endl;
            cin>>entry.car_model.Model;
            cout<<"Enter Car Manufacturer: (e.g. GM)"<<endl;
            cin>>entry.car_model.Manufacturer;
            cout<<"Enter Car VIN: (e.g. 2G1FB1ED4C9000000)"<<endl;
            cin>>entry.car_model.VIN;
            cout<<"Enter the year when the car was manufactured"<<endl;
            cin>>entry.car_model.Year;
            while(entry.car_model.Year>current_year){cout<<"Year Manufactured should be equal or smaller than "<<current_year<<endl;
                                                     cin>>entry.car_model.Year;};

            cout << "Enter the value of the car to be insured:  $"<<endl;
			cin >> carValue;
            while((current_year-entry.car_model.Year<=5)&&(carValue<5000.)){cout<<"Car is too cheap! and Input value again:";
            	    cin>>carValue;
                     };
			while (carValue<0)
			{
				cout<<"You have entered an invalid number. Please input a positive value\n"
				    << "Enter the value of the car to be insured:  $"<<endl;
				cin >> carValue;
			} ;


			cout<< "How many accidents has the policy holder caused during last three years? "<<endl<<endl;
			cin >> accidents;
			if (accidents<0)
			{
				cout<<"It should be equal or larger than 0. The program set it to zero"<<endl;
				accidents=0;
			}



			cout<<endl << "Enter the geographical risk factor (Class a - f): "<<endl;
			cin >> riskSurcharge;

		};

           void check_calendar(int &month,int &date,int &year){
        	   cout<<month<<"/"<<date<<", "<<year<<endl;
               while((month>12)||(month<1)){cout<<"Month should be from 1 to 12, input again"<<endl;
            	       cin>>month;};
               int date1;string Month;
            switch(month){
              case 1: {date1=Jan; Month="January"; break;}
              case 2: {date1=Feb+check_leapYear(year);Month="February";break;}
              case 3: {date1=Mar; Month="March";break;}
              case 4: {date1=Apr; Month="April";break;}
              case 5: {date1=May; Month="May";break;}
              case 6: {date1=Jun; Month="June";break;}
              case 7: {date1=Jul; Month="July";break;}
              case 8: {date1=Aug; Month="August";break;}
              case 9: {date1=Sep; Month="September";break;}
              case 10:{date1=Oct; Month="October";break;}
              case 11:{date1=Nov; Month="November";break;}
              case 12:{date1=Dec; Month="December";break;}
              default: cout<<"wrong format!"<<endl;
              };
            while((date>date1)||(date<1)){cout<<"Date of "<<Month<<" should be from 1 to "<<date1<<", input again"<<endl;
         	       cin>>date;};
            while(year<1900){cout<<"Year should be greater than 1900, input again"<<endl;
         	       cin>>year;};

           };

           int check_leapYear(int year){ //checks whether the year is a leap year??
               if(year % 400 == 0 || (year % 100!=0 && year % 4 ==0))
                  return 1;
               return 0;
           };



