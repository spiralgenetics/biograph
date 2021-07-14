#ifndef IO__DOUBLE_PAIR_H__
#define IO__DOUBLE_PAIR_H__

#include "modules/io/transfer_object.h"


struct double_pair
{
public:
	double_pair() : first(0.0), second(0.0) {}
	double_pair(double the_first, double the_second) : first(the_first), second(the_second) {}
	
	double_pair& operator+=(const double_pair &rhs)
	{
		first += rhs.first;
		second += rhs.second;
		
		return *this;
	}

	TRANSFER_OBJECT { 
		VERSION(0);
		FIELD(first); 
		FIELD(second); 
	}
	
	double first;
	double second;
};

#endif // IO__DOUBLE_PAIR_H__
