// set.h: interface for the set class.
// Here is a header file for the implementation of the set
//////////////////////////////////////////////////////////////////////

#if !defined(SET_CLASSES)
#define SET_CLASSES

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "typedefs.h"			//use basic types

//type of operation for the BitOperation()
#define		OR			0
#define		AND			1
#define		XOR			2

class set  
{
public:
	set();
	set(const set &);
	set(SIGNED_4B_TYPE);
	
	//some special constructors
	set(apvector <UNSIGNED_2B_TYPE> &);
	set(apvector <SIGNED_4B_TYPE> &);
	set(SIGNED_4B_TYPE, SIGNED_4B_TYPE);

	virtual ~set();

public:
	const set & operator = (const set &);
	set operator | (const set &);		//union
	set operator & (const set &);		//intersection
	set operator ^ (const set &);		//addition (union of unique parts)
	set operator - (const set &);		//subtraction (unique part)

	const set & operator |= (set &);
	const set & operator &= (set &);
	const set & operator ^= (set &);
	const set & operator -= (set &);
	

	bool IsInSet(set &);					//checks if the given set is a subset of the set
	bool IsInSet(SIGNED_4B_TYPE);			//checks whether the element is in the set

	void PutInSet(SIGNED_4B_TYPE);			//adds element to the set
	void RemoveFromSet(SIGNED_4B_TYPE);		//removes element from the set
	
	SIGNED_4B_TYPE Size();					//returns number of elements in the set
	bool IsEmpty();							//checks whether the set is empty
	void Dump();							//clears the set
	void GetList(apvector<SIGNED_4B_TYPE> &);//returns all elements of the set in the array
	void GetList(apvector<UNSIGNED_2B_TYPE> &);
	bool GetElement(SIGNED_4B_TYPE &);
	bool GetElement(UNSIGNED_4B_TYPE &);

	UNSIGNED_4B_TYPE SaveSet(FILE *);
	UNSIGNED_4B_TYPE LoadSet(FILE *);

private:
	void deflate();
	void BitOperation(const set &, const set &, UNSIGNED_1B_TYPE);		//invisible function that handles all bitwise operations
	
	apvector<UNSIGNED_1B_TYPE> S;
};

#endif // !defined(SET_CLASSES)
