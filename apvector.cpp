/*******************************************************************
6/16/2009 - rand_shuffle()

APCS vector class  IMPLEMENTATION
see vector.h for complete documentation of functions
vector class consistent with a subset of the standard C++ vector class
as defined in the draft ANSI standard (part of standard template library)
*******************************************************************/


#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include "apvector.h"

template <class itemType>
apvector<itemType>::apvector()
//postcondition: vector has a capacity of 0 items, and therefore it will
//               need to be resized
    : mySize(0),
      myList(0)
{

}

template <class itemType>
apvector<itemType>::apvector(int size)
// precondition: size >= 0
// postcondition: vector has a capacity of size items
   : mySize(size),
	 myList(0)
{
	if (size > 0)
		myList = new itemType[size];
}

template <class itemType>
apvector<itemType>::apvector(int size, const itemType & fillValue)
// precondition: size >= 0
// postcondition: vector has a capacity of size items, all of which are set
//                by assignment to fillValue after default construction
    : mySize(size),
      myList(0)      
{
	if (size == 0)
		return;

	myList = new itemType[size];    
    for(int k = 0; k < size; k++)
    {
        myList[k] = fillValue;
    }
}

template <class itemType>
apvector<itemType>::apvector(const apvector<itemType> & vec)
// postcondition: vector is a copy of vec
    : mySize(vec.length()),
      myList(0)
{
	if (mySize == 0)
		return;

	myList = new itemType[mySize];
    
        // copy elements
    for(int k = 0; k < mySize; k++)
	{
        myList[k] = vec.myList[k];
    }
}

template <class itemType>
apvector<itemType>::~apvector ()
// postcondition: vector is destroyed
{
	if (mySize > 0)
		delete [] myList;
}

template <class itemType>
const apvector<itemType> &
apvector<itemType>::operator = (const apvector<itemType> & rhs)
// postcondition: normal assignment via copying has been performed;
//                if vector and rhs were different sizes, vector
//                has been resized to  match the size of rhs
{
    if (this != &rhs)                           // don't assign to self!
    {
		if (mySize > 0)
			delete [] myList;                       // get rid of old storage

        mySize = rhs.length();
        if (mySize > 0)
			myList = new itemType [mySize];         // allocate new storage

         // copy rhs        
        for(int k=0; k < mySize; k++)
        {
            myList[k] = rhs.myList[k];
        }
    }
    return *this;                               // permit a = b = c = d
}

template <class itemType>
int apvector<itemType>::length() const
// postcondition: returns vector's size (number of memory cells
//                allocated for vector)
{
    return mySize;
}

template <class itemType>
itemType & apvector<itemType>::operator [] (int k)
// description: range-checked indexing, returning kth item
// precondition: 0 <= k < length()
// postcondition: returns the kth item
{

    if (k < 0 || mySize <= k)
    {
        cerr << "Illegal vector index: " << k << " max index = ";
        cerr << (mySize-1) << endl;
        exit(1);
    }
    return myList[k];
}

template <class itemType>
const itemType & apvector<itemType>::operator [] (int k) const
// safe indexing, returning const reference to avoid modification
// precondition: 0 <= index < length
// postcondition: return index-th item
// exception: exits if index is out-of-bounds
{
    if (k < 0 || mySize <= k)
    {
        cerr << "Illegal vector index: " << k << " max index = ";
        cerr << (mySize-1) << endl;
        exit(1);
    }
    return myList[k];
}

template <class itemType>
void apvector<itemType>::resize(int newSize)
// description:  resizes the vector to newSize elements
// precondition: the current capacity of vector is length(); newSize >= 0
// postcondition: the current capacity of vector is newSize; for each k
//                such that 0 <= k <= min(length, newSize), vector[k]
//                is a copy of the original; other elements of vector are
//                initialized using the 0-argument itemType constructor
//                Note: if newSize < length, elements may be lost
{
    int k;
    int numToCopy = newSize < mySize ? newSize : mySize;

    // allocate new storage and copy element into new storage

    itemType * newList = NULL;
	if (newSize > 0)
	{
		newList = new itemType[newSize];	
		for(k=0; k < numToCopy; k++)    
			newList[k] = myList[k];
	}    

	if (mySize > 0)
		delete [] myList;                      // de-allocate old storage

    mySize = newSize;                      // assign new storage/size
    myList = newList;
}

template <class itemType>
void apvector<itemType>::rand_shuffle()
{//NB: if (mySize > RAND_MAX) shuffling may not be effective
	int c, l, k = (mySize >> 1);
	itemType T;
	while (k--)
	for (c = 0; c < mySize; c++) 
	{//swap c with the random element l
		l = (rand() % mySize);
		T			 = myList[c];
		myList[c]	 = myList[l];
		myList[l]	 = T;
	}
}
