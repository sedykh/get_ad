// stack.h: interface for the stack class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(STACK_CLASS)
#define STACK_CLASS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "typedefs.h"			//basic types

class stack  
{
public:
	stack();
	stack(const stack &);
	virtual ~stack();

	UNSIGNED_4B_TYPE Pop();
	void Push(UNSIGNED_4B_TYPE);
	bool IsEmpty();
	void Dump();

private:
	apvector<UNSIGNED_4B_TYPE> S;
};

#endif // !defined(STACK_CLASS)
