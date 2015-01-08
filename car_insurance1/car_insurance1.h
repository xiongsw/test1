/*
 * car_insurance1.h
 *
 *  Created on: Dec 2, 2014
 *      Author: Shangwu Xiong
 */


#ifndef CAR_INSURANCE1_H_
#define CAR_INSURANCE1_H_
#include <iostream>

enum Calendar_date{Jan=31,Feb=28,Mar=31,Apr=30,May=31,Jun=30,Jul=31,Aug=31,Sep=30,Oct=31,Nov=30,Dec=31};

  struct car{
	  std::string Model;
	  std::string Manufacturer;
      std::string VIN;
      int Year;
      };
  struct date{
		int date, month, year;
	    };
///* a customer_data of the doubly linked list */
	struct customer_data {
		std::string FirstName; //!: First name of the account holder
		std::string LastName; //!: Last name of the account holder
		date DateofBirth ; //: Date of birth of the account holder
        date DateOfRecord; // Date the customer's insurance record is updated (today)
        date DateOfOpening; //: Date the account was opened
        date DateOfStart; //Date the insurance will start
		date DateOfExpiry ;//: Date after which the insurance will expire
		double RenewalCost; //: cost paid for the renewal of the insurance account
		double CarValue;
		double Accidents;
		char RiskSurcharge,Record_status;
        car car_model; //: Model,Manufacturer,VIN, Year	of a car
	    struct customer_data *next;
	    struct customer_data *prev;
	};

class car_insurance1{
private:
	int count;
	double carValue,accidents,premium;
// record_status='N=new, V=void, C=cancelled, R=revised'
	char riskSurcharge, record_status;
	const char *inputfile_name;
	customer_data *head;
	customer_data *tail;

public:
    enum region_risk_factor{A_a=20,B_b=30,C_c=40,D_d=55,E_e=70,F_f=100};
    const static double base0=30.0, kids=2.0;
    const static double accident_0=0., accident_1=0.01, accident_2=0.02,accident_3=0.03;

    car_insurance1(int count);
    car_insurance1(double carValue=1000.,double accidents=0,char riskSurcharge='0', customer_data *entry1=NULL);
//    car_insurance1(double carValue,double accidents,char riskSurcharge, customer_data *entry1);
    void insurance_Base();
    void insurance_RiskSurcharge();
    void insurance_accidentsAdjust();
    void ReviseRecord(int &choice,customer_data *entry1);
    void VoidRecord();
    void SaveRecords(customer_data *entry1);
    void printrecord(customer_data *entry1)const;

    double getcarValue(){return carValue;};
    double getaccidents(){return accidents;};
    double getpremium(){return premium;};


	virtual ~car_insurance1();


};



#endif /* CAR_INSURANCE1_H_ */
