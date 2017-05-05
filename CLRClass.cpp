/*
 * CLRClass.cpp
 *
 *  Created on: 2017年4月27日
 *      Author: clover
 */

#include "CLRClass.h"
#include <iostream>

CLRClass::CLRClass(int _a) {
	// TODO Auto-generated constructor stub
	m_a = _a;

	std::cout<<"Hello World " << m_a << std::endl;
}

CLRClass::~CLRClass() {
	// TODO Auto-generated destructor stub
	std::cout<<"Bye World " << m_a << std::endl;
}

void CLRClass::print_para  (){
	std::cout << m_a <<std::endl;

}
