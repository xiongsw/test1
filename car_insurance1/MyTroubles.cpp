/*
 * MyTroubles.cpp
 *
 *  Created on: Dec 3, 2014
 *      Author: Shangwu Xiong
 */
#include "MyTroubles.h"

// constructor for trouble
Trouble::Trouble(const char* pStr):pMessage(pStr){};

Trouble::~Trouble(){};

//Return the message
const char* Trouble:: what() const
{
	return pMessage;
	};

//constructor for Cannot_insured
Cannot_insured::Cannot_insured(const char* pStr):Trouble(pStr){};
