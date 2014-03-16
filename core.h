//core.h contains basic types definitions and some global constants 
//which are used throughout all files of core-project
//A.Sedykh, February 28, 2001-2010

#if !defined(CORE_)
#define CORE_

#include "typedefs.h"			//basic types and libraries

#include "bonds.h"
#include "list.h"
#include "set.h"
#include "stack.h"
#include "hash.h"

//---------------------------------------------------------------------
//---------------------------------------------------------------------
#ifdef WIN_PORT2
#define CARBON_COLOR		RGB(140, 240, 240)
#define NITROGEN_COLOR		RGB(0, 0, 255)
#define OXYGEN_COLOR		RGB(255, 0, 0)
#define SULFUR_COLOR		RGB(240, 240, 200)
#define PHOSPHORUS_COLOR	RGB(160, 250, 100)
#define CHLORINE_COLOR		RGB(230, 170, 20)
#define BROMINE_COLOR		RGB(230, 140, 120)
#define FLOURINE_COLOR		RGB(220, 190, 180)
#define DEFAULT_COLOR		RGB(255, 255, 255)
#define BACKGROUND_COLOR	RGB(0, 0, 0)

#else

#ifndef COLORREF
	typedef UNSIGNED_4B_TYPE	COLORREF;
#endif

#ifndef LPVOID
	typedef GENERIC_POINTER		LPVOID;
#endif

#define CARBON_COLOR		0x8cf0f0
#define NITROGEN_COLOR		0x0000ff
#define OXYGEN_COLOR		0xff0000
#define SULFUR_COLOR		0xf0f0c8
#define PHOSPHORUS_COLOR	0xa0fa64
#define CHLORINE_COLOR		0xe6aa14
#define BROMINE_COLOR		0xe68c78
#define FLOURINE_COLOR		0xdcbeb4
#define DEFAULT_COLOR		0xffffff
#define BACKGROUND_COLOR	0x000000

#endif
//---------------------------------------------------------------------
//---------------------------------------------------------------------


typedef apvector<set> apvector_set_type;

struct dat
{
	UNSIGNED_1B_TYPE	CONFG;		//configuration
	REALNUM_TYPE		COV_R;		//covalent radius
	REALNUM_TYPE		E_NGTV;		//electro-negativity
};

struct state 
{
	UNSIGNED_1B_TYPE VAL;
	apvector<dat> CFG;
};

struct Elements
{
	STRING_TYPE			SHORTNAME;
	STRING_TYPE			FULLNAME;
	UNSIGNED_2B_TYPE	ATOMIC_NUMBER;
	REALNUM_TYPE		ATOMIC_WEIGHT;
	REALNUM_TYPE		ATOMIC_RADIUS;
	REALNUM_TYPE		ATOMIC_ENGTV;				//electronegativity
	REALNUM_TYPE		ATOMIC_IP1, ATOMIC_IP2;		//ionization potentials
	apvector<state>		STATES;
	COLORREF			COLOR;
};


class AtomTypeHash : public hash<UNSIGNED_4B_TYPE>
{
public:	
	UNSIGNED_4B_TYPE hashFunction(UNSIGNED_4B_TYPE);

};


// ---------- Ring structures  ---------- //


// ---- for the array of rings
struct Ring
{
	UNSIGNED_1B_TYPE			done;
	set							setRing;
	apvector<SIGNED_4B_TYPE>	pRing;
	list<SIGNED_4B_TYPE>		lRing;	
	
	//flags
	bool						flat;					
	bool						pi_s;			//	a connected system of pi-orbitals
	bool						arom;
	bool						dlcl;			//	when it's bonds can be set to AROMATIC-BOND without problems.
};

struct CEdge			//a structure for storing ring edges
{UNSIGNED_4B_TYPE A1, A2;};	



struct atomic_in4cycl
{
	apvector<SIGNED_4B_TYPE> cycles;
};
// -------- End of ring structures -------- /


//--------    local structure that is needed in graph similarity search    -----------//
//			  used in MaximalCommonSubstructure(), CheckFragment() of the molAtom module
struct Cpair {	UNSIGNED_4B_TYPE A, B, P;};		//a pair of corresponding atoms and the degree of their likeness
struct Chain { set setA, setB, setCpair;};		//Chain of corresponding atoms


//------------------    Functions' predeclarations   ----------------------//

//ring handling structures and subroutines
void CleanRing(Ring &);
void InitializeRing(Ring &);
void InitializeRings(apvector<Ring> &);
void FindRingSystems(apvector<Ring> &, apvector<SIGNED_4B_TYPE> &, UNSIGNED_1B_TYPE = 2);
void MinimizeRings(apvector<Ring> &, apvector<SIGNED_4B_TYPE> &);

void FinalizeRings(apvector<Ring> &, apvector<SIGNED_4B_TYPE> &, SIGNED_4B_TYPE);

void SetUpColors();
bool LoadElements();
char CheckResidue(STRING_TYPE);
void GeneratePrimeNumbers(apvector<UNSIGNED_4B_TYPE> &, UNSIGNED_4B_TYPE);

void PutInLogFile(STRING_TYPE); //saves data into the default logfile, added 2008
void GetTimeStamp(STRING_TYPE &);

//---------------------------------------------------------
//Windows based
#ifdef WIN_PORT2
bool Message(STRING_TYPE, bool = false);
bool GetFile(STRING_TYPE &,	//default name of file and result
			 STRING_TYPE ,	//title of operation
			 STRING_TYPE ,	//name of filter - for example - datafiles
			 STRING_TYPE , 
			 HWND =NULL  );	//default extention and also standard filter

bool PutFile(STRING_TYPE &,	//default name of file and result
			 STRING_TYPE ,	//title of operation
			 STRING_TYPE ,	//name of filter - for example - datafiles
			 STRING_TYPE ,
			 HWND = NULL ,	//default extention and also standard filter
			bool = false);	//mode for loading multiple files
#endif //#ifdef WIN_PORT2
//---------------------------------------------------------

bool CheckStrEnding(STRING_TYPE &, STRING_TYPE);
void CutStrEnding(STRING_TYPE &);

REALNUM_TYPE String2Number(STRING_TYPE &);
void SplitString(STRING_TYPE &, STRING_TYPE, apvector<STRING_TYPE> &);
SIGNED_4B_TYPE SaveSetAsText(FILETYPE_OUT &, set &, SIGNED_4B_TYPE = 1);
SIGNED_4B_TYPE LoadSetAsText(FILETYPE_IN &, set &, SIGNED_4B_TYPE = 1);

REALNUM_TYPE getMetricDistance(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &, REALNUM_TYPE = 2.0, UNSIGNED_1B_TYPE = 0);

REALNUM_TYPE RRound(REALNUM_TYPE);
SIGNED_4B_TYPE Round(REALNUM_TYPE);

UNSIGNED_4B_TYPE GetRandomNumber(UNSIGNED_4B_TYPE);
void SortRandomly (apvector<REALNUM_TYPE> &);
UNSIGNED_4B_TYPE GetCRC32 (const GENERIC_POINTER, UNSIGNED_4B_TYPE);

//finds a value in a sorted array
SIGNED_4B_TYPE FindArrPoz(apvector<SIGNED_4B_TYPE> &, SIGNED_4B_TYPE);

//returns a subset of elements, every call gives next combination
bool GetCombination(apvector<SIGNED_4B_TYPE> &, apvector<SIGNED_4B_TYPE> &, UNSIGNED_2B_TYPE = 0);
bool GetCombination(set &, set &, UNSIGNED_2B_TYPE = 0);

void BubbleSort(apvector<REALNUM_TYPE> &, apvector<SIGNED_4B_TYPE> &);

//--- Qsort interface
int QSortCompareGreater(const void *, const void *);
int QSortCompareLess(const void *, const void *);
extern REALNUM_TYPE *QSortScore;


//------------------    SHARED VARIABLES    ------------------ //
extern apvector<Elements> ATABLE;	//table of atomic elements

extern set OrganicElements, AlkaliMetals, Metals, Halogens, Hbonders,
AromaticityElements, Group5Elements, Group6Elements, 
Loc, Deloc, PiBonds, CovBonds,
HphobResidues, HphilResidues, AmbiResidues, AromResidues;

extern bool SupressMessages;

extern STRING_TYPE ModulePath;		//a full path to the place where executable and configuration files are
extern STRING_TYPE StartDirName;	//current directory
extern STRING_TYPE LOG_FILENAME;	//current log file name

//ring variables
extern apvector <Ring>	  rings, backup_rings;	//structure to store rings, used by molecule-class
extern apvector<SIGNED_4B_TYPE>  ringSys;		//structure to store ring-systems (it points to rings[])
extern apvector <Ring>	 delocs, contours;		//stucture to store pi-contours
extern apvector <atomic_in4cycl> CMAtoms;		//to store cyclic information

//structures for finding and storing hits\subgraphs, etc.
extern apvector<Cpair> A2Corrs;
extern apvector<Chain> Corrs;

#endif		//#define CORE_
