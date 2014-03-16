// hash.cpp: implementation of the one-dimentional hash table.
//			 When more dimentions are required, the Hash function should
//			 be adjusted so that it should unfold and project all x-D range
//			 of values into 1D.
//
//////////////////////////////////////////////////////////////////////
//KEYS[] contain sorted hash values.
//INDX[i] stores the position in HASH[] of the KEYS[i] hash value.
//
//thus, HASH[INDX[i]] stores a list of values b, for each of which
//      KEYS[i] is equal to hashFunction(b).
//////////////////////////////////////////////////////////////////////

#include "hash.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template <class HashValueType>
hash<HashValueType>::hash()
{	
	nHASH = ZERO;
	FixedExpand = DEF_FIXED_EXPAND;
	FlexbExpand = DEF_FLEXB_EXPAND;
}

template <class HashValueType>
hash<HashValueType>::hash(UNSIGNED_4B_TYPE fix, REALNUM_TYPE flex)
{
	nHASH = ZERO;
	FixedExpand = fix;
	FlexbExpand = flex;
}

template <class HashValueType>
hash<HashValueType>::hash(const hash & Hh)
: HASH(Hh.HASH), KEYS(Hh.KEYS), INX(Hh.INX)
{	
	nHASH		= Hh.nHASH;
	FixedExpand = Hh.FixedExpand;
	FlexbExpand = Hh.FlexbExpand;
}

template <class HashValueType>
hash<HashValueType>::~hash()
{	
}

template <class HashValueType>
UNSIGNED_4B_TYPE hash<HashValueType>::FindHashKeyPoz(HashValueType V)
//description:		returns position in KEYS where hashkey V should be
//
//perconditon:		nHASH must be > ZERO!
//
//postcondition:		none
{

	UNSIGNED_4B_TYPE a = ZERO, e = nHASH, m;

	while (e > a + 1)
	{
		m = (a + e) >> 1;
		if ( KEYS[m] < V)
			a = m;
		else
			e = m;
	};

	if (KEYS[a] >= V)			//correction for a deadlock
		return a;
	
	return e;
}

template <class HashValueType>
UNSIGNED_4B_TYPE hash<HashValueType>::FindHashKey(HashValueType V)
//returns position in HASH corresponding to hashkey V
// (via KEYS -> INX -> HASH)
{
	UNSIGNED_4B_TYPE a = ZERO;

	if (nHASH == ZERO)	
		return ZERO;

	a = FindHashKeyPoz(V);

	if (nHASH == a)
		return nHASH;
	
	if (KEYS[a] != V)
		return nHASH;

	return INX[a];
}

template <class HashValueType>
void hash<HashValueType>::AddRecord(UNSIGNED_4B_TYPE Index)
//description:		adds Index as a new hash-element into a proper devision,
//					for that hashFunction is called
{	
	HashValueType Rv = hashFunction(Index);
	UNSIGNED_4B_TYPE i, Poz, Key = FindHashKey(Rv);
	
	if (Key == nHASH)
	{		
		if ( nHASH == UNSIGNED_4B_TYPE(HASH.length()) )
		{
			UNSIGNED_4B_TYPE E = (UNSIGNED_4B_TYPE)FlexbExpand*nHASH;
			HASH.resize(nHASH + FixedExpand + E);
			KEYS.resize(nHASH + FixedExpand + E);
			INX.resize(nHASH + FixedExpand + E);
		}

		//Key insertion
		if (nHASH == ZERO)
			Poz = ZERO;
		else
			Poz = FindHashKeyPoz(Rv);

		for (i = Key; i > Poz; i--)
		{
			KEYS[i] = KEYS[i-1];
			INX[i] = INX[i-1];
		}
		KEYS[Poz] = Rv;
		INX[Poz] = Key;
		nHASH++;
	}

	HASH[Key].Insert(Index);	
}

template <class HashValueType>
UNSIGNED_4B_TYPE hash<HashValueType>::RetrieveRecord(UNSIGNED_4B_TYPE KeyIndex, UNSIGNED_4B_TYPE &Index)
//description:		operates on the local hash-list for the hash-entry KeyIndex
//					i.e. a current hash-element is retrieved and returned
//precondition:		
//postcondition:	returns size of the list, and the current element in Index
{
	UNSIGNED_4B_TYPE N = nHASH, Sz = ZERO;
	
	if (KeyIndex < N)
	{
		Sz = HASH[KeyIndex].Size();

		if(Sz != ZERO)
			Index = HASH[KeyIndex].Next();
	};

	return (Sz);
}

template <class HashValueType>
void hash<HashValueType>::Wipe()
{
	nHASH = ZERO;
	HASH.resize(ZERO);
	KEYS.resize(ZERO);
	INX.resize(ZERO);
}

template <class HashValueType>
void hash<HashValueType>::SaveHash(FILE * hf)
//description:		saves hash table into a binary file
//precondition:		hf must be open for writing in binary mode
{
	UNSIGNED_4B_TYPE x, hi, 
		bsz = sizeof(HashValueType), hsz = sizeof(UNSIGNED_4B_TYPE);
	
	fwrite(&bsz, hsz, 1, hf);							//datablock size
	fwrite(&FixedExpand, hsz, 1, hf);					//FixedExpand
	fwrite(&FlexbExpand, sizeof(REALNUM_TYPE), 1, hf);	//FlexbExpand;
	 
	fwrite(&nHASH, hsz, 1, hf);//# hash entries
	for (x = ZERO; x < nHASH; x++)
	{
		fwrite(&(KEYS[x]), bsz, 1, hf);
		fwrite(&(INX[x]), hsz, 1, hf);
	}

	for (x = ZERO; x < nHASH; x++)
	{
		bsz = HASH[x].Size();
		fwrite(&bsz, hsz, 1, hf);
		while (bsz)
		{
			hi = HASH[x].Next();
			fwrite(&hi, hsz, 1, hf);
			bsz--;
		}
	}
}

template <class HashValueType>
bool hash<HashValueType>::LoadHash(FILE * hf)
//description:		loads hash table from a binary file
//precondition:		hf must be open for reeading in binary mode
{
	UNSIGNED_4B_TYPE x, hi, bsz, hsz = sizeof(UNSIGNED_4B_TYPE);

	Wipe();
	fread(&bsz, hsz, 1, hf);							//datablock size
	if (bsz != sizeof(HashValueType))
		return false;

	fread(&FixedExpand, hsz, 1, hf);					//FixedExpand
	fread(&FlexbExpand, sizeof(REALNUM_TYPE), 1, hf);	//FlexbExpand;
	 
	fread(&nHASH, hsz, 1, hf);//# hash entries
	KEYS.resize(nHASH);
	INX.resize(nHASH);
	HASH.resize(nHASH);
	for (x = ZERO; x < nHASH; x++)
	{
		fread(&KEYS[x], bsz, 1, hf);
		fread(&INX[x], hsz, 1, hf);
	}

	for (x = ZERO; x < nHASH; x++)
	{		
		fread(&bsz, hsz, 1, hf);
		while (bsz)
		{			
			fread(&hi, hsz, 1, hf);
			HASH[x].Insert(hi);
			bsz--;
		}
	}

	return true;
}
