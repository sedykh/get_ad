// set.cpp: implementation of the set class.
// Implementaion of the set as a mathematical abstraction :) 
//////////////////////////////////////////////////////////////////////


#include "set.h"


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

set::set()
{
	S.resize(ZERO);
}

set::set(SIGNED_4B_TYPE i)
{
	PutInSet(i);
}

set::~set()
{
}

set::set(const set &initSet)
: S(initSet.S)
{
	//	Copying the contents	
}

//some special constructors
set::set(apvector <UNSIGNED_2B_TYPE> &Filler)
{
	for (UNSIGNED_2B_TYPE di = ZERO; di < Filler.length(); di++)
		PutInSet(Filler[di]);
}

set::set(apvector <SIGNED_4B_TYPE> &Filler)
{
	for (SIGNED_4B_TYPE di = ZERO; di < Filler.length(); di++)
		PutInSet(Filler[di]);
}

set::set(SIGNED_4B_TYPE Min, SIGNED_4B_TYPE Max)
{
	for (SIGNED_4B_TYPE di = Min; di < Max; di++)
		PutInSet(di);
}

//-------------------- structure 'set' subroutines ------------------------- //
void set::deflate()
{//service subroutine to automaticaly downsize the storage
	UNSIGNED_2B_TYPE sz = S.length();
	while (sz) 
	if (S[--sz]) 
	{
		S.resize(++sz);
		return;
	}

	S.resize(sz);
}

bool set::IsInSet(set &SI)
{
	SIGNED_4B_TYPE Byte, L = SI.S.length();

	if (L == ZERO)
		return true;

	if (L > S.length())
		return false;

	for (Byte = ZERO; Byte < L; Byte++)
	if ((S[Byte] & (SI.S[Byte])) == SI.S[Byte])
		continue;
	else
		return false;

	return true;
}

bool set::IsInSet(SIGNED_4B_TYPE index)
{
	SIGNED_4B_TYPE Byte, position, Bit;

	Byte = index >> 3;	//index / 8
	if ((Byte >= S.length()) || (index < ZERO))
		return false;

	
	position	= index - (Byte << 3);
	Bit			= 1 << position;

	if (S[Byte] & Bit)
		return true;

	return false;
}

void set::RemoveFromSet(SIGNED_4B_TYPE index)
{	
	SIGNED_4B_TYPE Byte, position, Bit;

	if (index < 0)
		return;

	Byte = index >> 3;

	//if set-variable is too small then no need to delete
	position = S.length();
	if (Byte < position)		
	{	
		position  = index - (Byte << 3);
		Bit =  1 << position;

		if ((S[Byte] & Bit) == Bit)
		{
			S[Byte] -= Bit;
			deflate();
		}
	}
}

void set::PutInSet(SIGNED_4B_TYPE index)
{	
	SIGNED_4B_TYPE i, Byte, position, Bit;

	if (index < 0)
		return;

	Byte = index >> 3;		//index / 8

	//if set-variable is too small then expand!
	position = S.length();
	if (Byte >= position)
	{
		S.resize(Byte+1);
		for (i= position; i<=Byte; i++)
			S[i] = ZERO;
	};
	
	position  = index - (Byte << 3);
	Bit =  1 << position; //2^position

	S[Byte] = (S[Byte] | Bit);
}

const set & set::operator = (const set &S1)
{
	S = S1.S;

	return (*this);
}

void set::BitOperation(const set &S1, const set &S2, UNSIGNED_1B_TYPE TypeOperation)
//
//NOTE:			this routine should work even if S1 or/and S2 is equal to *this
//
{
	UNSIGNED_1B_TYPE B1, B2;	
	SIGNED_4B_TYPE s1l = S1.S.length(), s2l = S2.S.length();
	S.resize(max(s1l, s2l));
	
	for (SIGNED_4B_TYPE g = 0; g<S.length(); g++)	
	{
		if (g < s1l)	B1 = S1.S[g];	else	B1 = ZERO;
		if (g < s2l)	B2 = S2.S[g];	else	B2 = ZERO;
			
		switch (TypeOperation)
		{
			case OR:
				S[g] = (B1 | B2);
				break;
			case AND:
				S[g] = (B1 & B2);
				break;
			case XOR:
				S[g] = (B1 ^ B2);
				break;
			default:
				break;
		};//switch

	};//for

	deflate();
}

set	set::operator | (const set &S1)
{//UNITE two sets
	set R;	
	R.BitOperation((*this), S1, OR);	
	return (R);
}

set set::operator & (const set & S1)
{
	set R;	
	R.BitOperation((*this), S1, AND);	
	return (R);
}

set set::operator ^ (const set & S1)
{
	set R;	
	R.BitOperation((*this), S1, XOR);	
	return (R);
}

set set::operator - (const set & S1)
{
	set R, S2;	
	S2.BitOperation((*this), S1, AND);
	R.BitOperation ((*this), S2, XOR);
	return (R);
}



const set & set::operator |= (set &S1)
{
	BitOperation((*this), S1, OR);
	return (*this);
}

const set & set::operator &= (set &S1)
{
	BitOperation((*this), S1, AND);
	return (*this);
}

const set & set::operator ^= (set &S1)
{
	BitOperation((*this), S1, XOR);
	return (*this);
}

const set & set::operator -= (set &S1)
{
	set X;
	X.BitOperation((*this), S1, AND);
	BitOperation((*this), X, XOR);
	return (*this);
}

bool set::IsEmpty()
{	//it's faster than usage of set::Size!

	for (SIGNED_4B_TYPE g = 0; g<S.length(); g++)	
		if (S[g] != ZERO)
			return false;

	return true;
}

SIGNED_4B_TYPE set::Size()
{
	SIGNED_4B_TYPE C = ZERO;
	
	for (SIGNED_4B_TYPE g = 0; g<S.length(); g++)			
		C +=	(S[g]& 1)		  + 
				((S[g]& 2)	>> 1) + 
				((S[g]& 4)	>> 2) + 
				((S[g]& 8)	>> 3) +
				((S[g]& 16)	>> 4) + 
				((S[g]& 32)	>> 5) + 
				((S[g]& 64)	>> 6) + 
				((S[g]& 128)>> 7);
	
	return (C);
}

void set::Dump()
{
	S.resize(ZERO);
}

void set::GetList(apvector<SIGNED_4B_TYPE> &L)
{
	SIGNED_4B_TYPE C = ZERO, p, g, k;
	L.resize(C);

	for (g = 0; g<S.length(); g++)	
		for (k = 0; k<8; k++)
		{
			p = 1 << k;
			if (S[g] & p)
			{
				C++;
				L.resize(C);
				L[C-1] = g*8+k;
			};
		};	
}

void set::GetList(apvector<UNSIGNED_2B_TYPE> &L)
{
	UNSIGNED_2B_TYPE C = ZERO, p, g, k;
	L.resize(C);

	for (g = 0; g<S.length(); g++)	
		for (k = 0; k<8; k++)
		{
			p = 1 << k;
			if (S[g] & p)
			{
				C++;
				L.resize(C);
				L[C-1] = g*8+k;
			};
		};	
}

bool set::GetElement(SIGNED_4B_TYPE &C)
{
	SIGNED_4B_TYPE p, g, k;	

	for (g = 0; g<S.length(); g++)	
		for (k = 0; k<8; k++)
		{
			p = 1 << k;
			if (S[g] & p)
			{
				C = g*8+k;
				return true;
			}
		};

	return false;
}

bool set::GetElement(UNSIGNED_4B_TYPE &C)
{
	SIGNED_4B_TYPE p, g, k;	

	for (g = 0; g<S.length(); g++)	
		for (k = 0; k<8; k++)
		{
			p = 1 << k;
			if (S[g] & p)
			{
				C = g*8+k;
				return true;
			}
		};

	return false;
}

UNSIGNED_4B_TYPE set::SaveSet(FILE * wbf)
//description:		saves set into the binary file wbf
//
//precondition:		wbf must be open in binary writing mode
//					set should not have more than 16k elements
//
//postcondition:	# of bytes written
{
	UNSIGNED_1B_TYPE Mod = ZERO; //default mode to store the set bitwise
	UNSIGNED_2B_TYPE sz, L = S.length();
	apvector<SIGNED_4B_TYPE> A;
	
	GetList(A);
	sz = A.length();
	sz <<= 2; //x4
	if (L > sz)
	{//store by elements;
		Mod = 1;
		L = A.length();
	}

	fwrite(&Mod, 1, 1, wbf); //way of storage
	fwrite(&L, 2, 1, wbf);

	if (Mod == ZERO)
	{
		for (sz = ZERO; sz < L; sz++)
			fwrite(&(S[sz]), 1, 1, wbf);	
	}
	else
	{
		for (sz = ZERO; sz < L; sz++)
			fwrite(&(A[sz]), 4, 1, wbf);

		L <<= 2;
	}

	return (L + 3);
}


UNSIGNED_4B_TYPE set::LoadSet(FILE * wbf)
//description:		loads set from the binary file wbf
//
//precondition:		wbf must be open in binary reading mode
//
//postcondition:	returns # of bytes read
{
	UNSIGNED_1B_TYPE Mod = ZERO;
	UNSIGNED_2B_TYPE sz = ZERO, x;
	SIGNED_4B_TYPE El;	

	fread(&Mod, 1, 1, wbf);
	fread(&sz, 2, 1, wbf);
	
	if (Mod == ZERO)
	{
		S.resize(sz);
		for (x = ZERO; x < sz; x++)
			fread(&S[x], 1, 1, wbf);
	}
	else
	{
		S.resize(ZERO);
		for (x = ZERO; x < sz; x++)
		{
			fread(&El, 4, 1, wbf);
			PutInSet(El);
		}

		sz <<= 2;
	}
	
	return (sz + 3);
}

// ------------------------- end of set subroutines ----------------------------- //
