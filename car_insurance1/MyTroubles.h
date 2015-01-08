/*
 * MyTroubles.h
 *
 *  Created on: Dec 3, 2014
 *      Author: Shangwu Xiong
 */
// reference: Ivor Horton (Beginning C++: The complete Language, 1998)
#ifndef MYTROUBLES_H_
#define MYTROUBLES_H_

// base exception class
class Trouble
{
public:
	Trouble(const char* pStr="There is a problem");
	virtual ~Trouble();
	virtual const char* what() const;

private:
	const char* pMessage;
	};
// derived exception class
class Cannot_insured: public Trouble
{
public:
	Cannot_insured(const char* pStr="too many accidents to insure your vehicle");
};



#endif /* MYTROUBLES_H_ */
