/*
Class-library for the methods for separating a dataset into a test and training subsets.

Coordinators: Aleksandr Sedykh, Alexander Golbraikh, Christopher Grulke

sfexcl_split() - Spherical Exclusion algorithm
ref:			Golbraikhm A, Tropsha A et al, J Comp-Aid Mol Design 2003, 17: 241-253
Below is the description of typical behavioral patterns of selecting new spheres:
Selection of new spheres is done either randomly or by a two-step selection process:
step1 (inner cycle), step2(outer cycle)
if SFEXCL_NEXTSF_SPHERES is on, then in the inner cycle spheres are scanned, 
and unassigned points - in the outer cycle, if SFEXCL_NEXTSF_SPHERES then it is the reverse

1) SFEXCL_NEXTSF_SPHERES is ON and STEP1/STEP2 is:
MIN/MIN		snake-like creeping addition of new spheres
MIN/MAX		lattice-like even equidistant coverage of space, empty spots of the lattice are covered the last
MAX/MAX		Covers corners, creeps along edges then closes in on the center
MAX/MIN		tumor-like growth of new spheres
SUM/MIN		tumor-like growth, like MAX/MIN but more dense
SUM/MAX		similar to MAX/MAX but more regular (usually quickly captures 4 corners first)
2) SFEXCL_NEXTSF_SPHERES is OFF and STEP1/STEP2 is:
MIN/MIN		same as in 1)
MIN/MAX		tumor-like, same as in 1) MAX/MIN
MAX/MAX		same as in 1)
MAX/MIN		growth is controlled by the existing sphere that is closest to the center of the descriptors space 
(because such point always has a minimal max-distance to an unassaigned datapoint). 
If initial starting point is in the corner, then the behavior is like that of MAX/MAX
SUM/MIN		same as MAX/MIN
SUM/MAX		similar to MAX/MAX
*/

#if !defined(DATASET_CLASSES)
#define DATASET_CLASSES

#include "qsar.h"

// -------------------------------------------------------------------------------------
//variables to store information about descriptors (labels, normalisation coefficients)
class descr_hash: public hash<UNSIGNED_4B_TYPE>
{//Hash to store descriptors
public:
	descr_hash(){};
	descr_hash(UNSIGNED_4B_TYPE R){ FixedExpand = R; };
	UNSIGNED_4B_TYPE hashFunction(UNSIGNED_4B_TYPE ix){ if (ix < UNSIGNED_4B_TYPE(L.length())) return GetCRC32(GENERIC_POINTER(L[ix].c_str()), L[ix].length()); return ZERO; };
	void Wipe() { L.resize(ZERO); A.resize(ZERO); B.resize(ZERO); hash<UNSIGNED_4B_TYPE>::Wipe(); };
	//nested data, stored here for convenience of handling
	apvector<STRING_TYPE> L;			//labels for descriptors -> {A, B}
	apvector<REALNUM_TYPE> A, B;		//rescaling coefficients for descriptors in hit-matrix (Pars)
};
// -------------------------------------------------------------------------------------

class id_hash: public hash<UNSIGNED_4B_TYPE>
{//Hash to store ids
public:
	id_hash(){};
	id_hash(UNSIGNED_4B_TYPE R){ FixedExpand = R; };
	UNSIGNED_4B_TYPE hashFunction(UNSIGNED_4B_TYPE ix){ if (ix < UNSIGNED_4B_TYPE(ids.length())) return GetCRC32(GENERIC_POINTER(ids[ix].c_str()), ids[ix].length()); return ZERO; };
	void Wipe() { ids.resize(ZERO); hash<UNSIGNED_4B_TYPE>::Wipe(); };
	//nested data, stored here for convenience of handling
	apvector<STRING_TYPE> ids;			//labels for descriptors -> {A, B}
};

//--------------------------------------------------------
//work mode options for sfexcl_split, combine them using |
//_R_ - algorithm for calculating sphere's radius
//DEFAULT, sphere's R is defined through descriptor space volume
#define	SFEXCL_R_BYDIST				1	//sphere's R is defined through min and max distances in the dataset
#define	SFEXCL_R_BYUSER				2	//R is directly supplied by the user, lowest priority

//NEXTSF - algorithms for selecting next sphere:
#define	SFEXCL_NEXTSF_RAND			4	//next sphere is picked randomly
#define	SFEXCL_NEXTSF_STEP1_MIN		8	//minimal distance is taken on step1, default is max d
#define	SFEXCL_NEXTSF_STEP2_MIN		16	//minimal distance is taken on step2, default is max d
#define	SFEXCL_NEXTSF_STEP1_SUMDIST	32	//sum of distances is taken on step1, applicable only with SFEXCL_NEXTSF_SPHERES
#define SFEXCL_NEXTSF_SPHERES		64	//spheres are analyzed (in step1) relative to candidate points (step2), default is the reverse
#define SFEXCL_SPLIT_FIRST2TRN		128	//in assigning points from a sphere: first batch should go training set (next one - to test set).

//modes for seeding initial spheres
#define	SFEXCL_SEED_BYACTS		512		//seeds points by activity
#define	SFEXCL_SEED_MINACT		1024	//adds min and max activity to the seeding set of points
#define	SFEXCL_SEED_MAXACT		256		//adds min and max activity to the seeding set of points

//default is Euclidean\Minkowski metric
#define SFEXCL_METRIC_COSINE	2048	//Cosine-based 1 - ()^metrikV
#define SFEXCL_METRIC_CORR		4096	//Correlation, similar to Cosine but mean-centered
#define SFEXCL_METRIC_TANIMOTO	8192	//Tanimoto, does not satisfy triangular inequiality
//--------------------------------------------------------

#define DUMMY_VAL	-999		//used for default value for distances, should be < 0

typedef apvector_set_type lneib;			//structure type for local neighbours
typedef apvector<SIGNED_4B_TYPE> dparr;		//structure type for datapoints array

class dataset 
{//dataset class implementation
public:
	dataset();
	dataset(const dataset &data_in);
	virtual ~dataset();
	void dump();
	//algorithms of splitting data onto test and training sets

	//random pick of nFraction points or rtFraction part as a test set; activity distribution can be preserved
	void rand_split(UNSIGNED_4B_TYPE nFraction, UNSIGNED_1B_TYPE nActClasses = 1, bool eqlSizeActBin = true, set *pExcl = NULL);
	void rand_split(REALNUM_TYPE rtFraction = 0.1, UNSIGNED_1B_TYPE nActClasses = 1, bool eqlSizeActBin = true, set *pExcl = NULL);
	
	REALNUM_TYPE sfexcl_split(set &seed, UNSIGNED_2B_TYPE Mode = 0, REALNUM_TYPE f = 0.2, bool ifForceDistCalc = false, UNSIGNED_1B_TYPE toTest = 1, UNSIGNED_1B_TYPE toTrain = 1, UNSIGNED_1B_TYPE nSeed = 0, REALNUM_TYPE Kmetric = 2.0);
	REALNUM_TYPE sfexcl_split(UNSIGNED_2B_TYPE Mode = 0, REALNUM_TYPE f = 0.2, bool ifForceDistCalc = false, UNSIGNED_1B_TYPE toTest = 1, UNSIGNED_1B_TYPE toTrain = 1, UNSIGNED_1B_TYPE nSeed = 0, REALNUM_TYPE Kmetric = 2.0);

	void lgo_split(lneib &sLGO, UNSIGNED_2B_TYPE nLGO, UNSIGNED_1B_TYPE nActClasses = 1, bool eqlSizeActBin = true);
	void lgo_split(lneib &sLGO, REALNUM_TYPE rtLGO = 0.1, UNSIGNED_1B_TYPE nActClasses = 1, bool eqlSizeActBin = true);

	void kohonen_split();								//kohonen maps
	void sca_split();									//stochastic clustering algorithm

//retrieval
   //friend istream& operator >> (istream& i, dataset& m);
   //friend ostream& operator << (ostream& o, dataset& m);

	bool load(istream& i, istream *i1 = NULL, UNSIGNED_1B_TYPE mtxType = 0);
	bool save(ostream& o, ostream *o1 = NULL, UNSIGNED_1B_TYPE mtxType = 0);

	void reduce_dimensions(set &);
	void scale_dimensions(UNSIGNED_1B_TYPE = 0);
	bool expandby(dataset &, bool = false);
	void normalizeby(dataset &);
	
	dataset subset(set &);
	dataset get_training_set();
	dataset get_test_set();

	SIGNED_4B_TYPE get_Ndatapoints();
	SIGNED_4B_TYPE get_Ndimensions();
	
	STRING_TYPE get_sid(SIGNED_4B_TYPE posN);		//get datapoint's id-name from its sequential number
	STRING_TYPE get_dscr(SIGNED_4B_TYPE posN);		//get descriptor's id-name from its sequential number
	SIGNED_4B_TYPE get_sid_pos(STRING_TYPE stID);	//get datapoint's sequential number from its id-name
	SIGNED_4B_TYPE get_dscr_pos(STRING_TYPE stID);	//get descriptor's sequential number from its id-name
	
	SIGNED_4B_TYPE ExpandDescriptors(dataset &);	//Aug, 2010
	void get_DimRowValues(SIGNED_4B_TYPE X, apvector<REALNUM_TYPE> &Dt);	//Feb, 2010
	void get_DimValues(SIGNED_4B_TYPE D, apvector<REALNUM_TYPE> &Dim);

	void get_ActValues(apvector<REALNUM_TYPE> &Acts);
	void set_ActValues(apvector<REALNUM_TYPE> &Acts);

	REALNUM_TYPE get_Act(SIGNED_4B_TYPE dp1);
	void set_Act(SIGNED_4B_TYPE dp1, REALNUM_TYPE val);

	REALNUM_TYPE get_MaxAct();
	REALNUM_TYPE get_MinAct();
	REALNUM_TYPE get_AverageAct();
	set get_ActPoints(REALNUM_TYPE lowAct, REALNUM_TYPE hiAct);
	
	void get_distMatr(matrix<REALNUM_TYPE> &DST);
	bool set_distMatr(matrix<REALNUM_TYPE> &DST);

	REALNUM_TYPE get_indDistance(SIGNED_4B_TYPE dp1, SIGNED_4B_TYPE dp2, set * pDims = NULL, REALNUM_TYPE metricV = 2.0, UNSIGNED_1B_TYPE metricKind = 0);
	REALNUM_TYPE get_Distance(SIGNED_4B_TYPE dp1, SIGNED_4B_TYPE dp2);
	
	REALNUM_TYPE get_MinNonZeroDistance(SIGNED_4B_TYPE dp = INVALID);
	REALNUM_TYPE get_MinDistance(SIGNED_4B_TYPE dp = INVALID);
	REALNUM_TYPE get_MaxDistance(SIGNED_4B_TYPE dp = INVALID);
	REALNUM_TYPE get_AverageDistance(SIGNED_4B_TYPE dp = INVALID, REALNUM_TYPE cutoff = DUMMY_VAL);
	
	//NB: this set of functions may help to restore respective parameters quickly,
	//be it is too cumbersome and calc_dist_pars() can be used instead
	void set_MinNonZeroDistance(REALNUM_TYPE rtD);
	void set_MinDistance(REALNUM_TYPE rtD);
	void set_MaxDistance(REALNUM_TYPE rtD);
	void set_AverageDistance(REALNUM_TYPE rtD);
	
	void calc_dist_pars();	//calculates min, max, aver, etc.
	void calc_dist(REALNUM_TYPE rtDef = 0, REALNUM_TYPE metricV = 2.0, UNSIGNED_1B_TYPE metricKind = 0);

	void add_dp(apvector<REALNUM_TYPE> &dscr, REALNUM_TYPE actdp = 0, STRING_TYPE sidp = "", STRING_TYPE pidp = "");
	void remove_dp(SIGNED_4B_TYPE dp);

	//returns nearest neighbours for a single datapoint datap
	void get_NearNeib(SIGNED_4B_TYPE datap, apvector<SIGNED_4B_TYPE> &NNeib, UNSIGNED_2B_TYPE kNeibours = 1, REALNUM_TYPE cutoffR = 0);

	//calculates stats for kNeibours throughout the dataset, calls get_NearNeib()
	void get_NearNeibDistances(apvector<REALNUM_TYPE> &Stats, UNSIGNED_2B_TYPE kNeibours = 1, REALNUM_TYPE cutoffR = 0, UNSIGNED_1B_TYPE kmode = 1);

	//special function for kNN mainly, calculates specified number of neighbours for each datapoint's each dimension
	//void calc_allDneibs(UNSIGNED_2B_TYPE nln, apvector<lneib> &dln, apvector<lneib> &dmln);

	void randomizeY(UNSIGNED_1B_TYPE nActClasses = 1, bool eqlSizeActBin = true);
	REALNUM_TYPE get_OccupVol( REALNUM_TYPE Probe, UNSIGNED_4B_TYPE EstimationMethod = 0, bool ifForceDistCalc = false, REALNUM_TYPE metricV = 2.0, UNSIGNED_1B_TYPE metricKind = 0);

	void sort_act();
	void get_act_bins(apvector<UNSIGNED_4B_TYPE> &abins, UNSIGNED_4B_TYPE nbins, bool ifbsame = true);

	void RemoveLowVarDims(REALNUM_TYPE = 0.01, UNSIGNED_2B_TYPE = 3);
	void RemoveHiCorrDims(REALNUM_TYPE = 0.99);

	bool isScaled();	//returns true if scaling coefficients exist.
public:
	set test, train;				//places to store test and training subsets

private:
	descr_hash vars;				//descriptor variables (labels, etc)

	apvector<REALNUM_TYPE> act;		//activities of data points
	matrix<REALNUM_TYPE> patt;		//pattern matrix of descriptors vs data points

	apvector<STRING_TYPE> pid;		//datapoint seq# or ids that don't have to be unique
	id_hash sid;					//hash of datapoints ids
	apvector<SIGNED_4B_TYPE> sact;	//datapoints sorted in increasing order of activities

	//service data which should be recalculated as rarely as possible! (due to heavy load)
	matrix<REALNUM_TYPE> dist;				//matrix of cross-distances for the datapoints
	REALNUM_TYPE distn0Min, distMin, distMax, distAv;	//min, max and average distances in dist() matrix
};//

#endif // !defined(DATASET_CLASSES)
