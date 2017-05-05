/*
 * CLRClass.h
 *
 *  Created on: 2017年4月27日
 *      Author: clover
 */

#ifndef CLRCLASS_H_
#define CLRCLASS_H_

class CLRClass {
	int m_a;
public:
	CLRClass(int _a);
	virtual ~CLRClass();
	//void operator << ();
	void print_para ();
};

#endif /* CLRCLASS_H_ */
