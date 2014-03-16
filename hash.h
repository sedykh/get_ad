// hash.h: interface for the hash class.
//
//NB: technically, any object can serve as an input-parameter to hash-function, but
//then the same object-type will have to be stored in a hashlist, which would inflate its size,
//hence currently the internal record is restricted to UNSIGNED_4B_TYPE. i.e., only (UNSIGNED_4B_TYPE)objects are hashed.
//since 4 byte pointer is compatible, it is possible to hash references to objects, instead.
//////////////////////////////////////////////////////////////////////

#if !defined(HASH_)
#define HASH_

#include "typedefs.h"
#include "list.h"

#define	DEF_FIXED_EXPAND		1024
#define	DEF_FLEXB_EXPAND		0.1

typedef slist<UNSIGNED_4B_TYPE>	hashlist;

template <class HashValueType>

class hash  
{
public:
	hash();	
	hash(UNSIGNED_4B_TYPE, REALNUM_TYPE);
	hash(const hash &);	
	virtual ~hash();

	void SaveHash(FILE *);
	bool LoadHash(FILE *);

	void Wipe();
	void AddRecord(UNSIGNED_4B_TYPE);	
	UNSIGNED_4B_TYPE RetrieveRecord(UNSIGNED_4B_TYPE, UNSIGNED_4B_TYPE &);
	UNSIGNED_4B_TYPE FindHashKey(HashValueType);
	UNSIGNED_4B_TYPE FindHashKeyPoz(HashValueType);
	
	virtual HashValueType hashFunction(UNSIGNED_4B_TYPE) = 0;
	
	UNSIGNED_4B_TYPE nHASH;
	apvector<hashlist> HASH;
	apvector<HashValueType> KEYS;			//stores hash keys sorted in ascenting order
	apvector<UNSIGNED_4B_TYPE> INX;			//stores pointers to HASH for respective hash values stored in KEYS

	UNSIGNED_4B_TYPE FixedExpand;
	REALNUM_TYPE FlexbExpand;
};

#include "hash.cpp"

#endif // !defined(HASH_)
