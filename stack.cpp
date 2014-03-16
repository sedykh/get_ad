// stack.cpp: implementation of the stack class.
//
//////////////////////////////////////////////////////////////////////
#include "stack.h"


//-------   memory leaks catcher for the current source-file  --------
#ifdef ADV_LEAK_CATCHER
#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif
//--------------------------------------------------------------------


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

stack::stack()
{
}

stack::stack(const stack &St)
: S(St.S)
{
}

stack::~stack()
{
}


UNSIGNED_4B_TYPE stack::Pop()
{
	UNSIGNED_4B_TYPE L = S.length();
	UNSIGNED_4B_TYPE C = ZERO;
	
	if (L > ZERO)
	{
		C = S[L-1];
		S.resize(L-1);
	}

	return (C);
}

void stack::Push(UNSIGNED_4B_TYPE C)
{
	UNSIGNED_4B_TYPE L = S.length();
	S.resize(L+1);
	S[L] = C;
}

bool stack::IsEmpty()
{
	return (S.length() == ZERO);
}

void stack::Dump()
{
	S.resize(ZERO);
}
// ------------------ end of stack subroutines -------------------//



