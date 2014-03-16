/*
TODO: 
(Add here problems and suggestions to be worked upon)


1)

2) implement SCA stochastic cluster analysis, but what to use as similarity?


DONE:
Aug 11 2011
			isScaled() added for high-level check for the presence of scaling in the dataset
May 26 2011	
			minor fix in load() to prevent silent handling of over-the-limit entries
			RemoveLowVarDims() expanded to perform filtering by min.occurrences

Nov 21 2010	
			ExpandDescriptors() fixed.
Oct 4 2010
			minor fix in rand_split() to make act-independent splits if nActBins = 1
Aug 20 2010
			added member-functions():
			void set_ActValues(apvector<REALNUM_TYPE> &Acts);
			void set_Act(SIGNED_4B_TYPE dp1, REALNUM_TYPE val);
			SIGNED_4B_TYPE ExpandDescriptors(dataset &X);
			
June 20 2010
			string buffer overflow in save(): Bfr[20] changed to Bfr[100]

May 16 2010
			R-calculation fix in sfexcl_split()

May 12 2010 
			svm load/save critical fix (dscr indexation was from 0, not from 1)!!

May 5 2010
			get_act_bins() fixed, wrong range-binning of continuous data!!
April 27 2010
			 minor fix in load(), empty lines in the file end caused wrong initialization of .vars.B() if #descriptors = 1
April 14 2010
			 critical fix in scale_dimensions(): vars.A[u] = INVALID; ( was "vars.A = INVALID" );
April 3 2010
			 adjusted group size added to lgo_split() to make splits more even.
March 30 2010 
			 minor fix: RRound() added into lgo_split(), 11 splits were generated instead of 10 for 127 compounds of binary activity

February 2010 
			  get_DimRowValues() added to retrieve values of a descriptor column

Feb 28 2010	  added:
			  normalizeby(), scale_dimensions(), reduce_dimensions()
			  RemoveLowVarDims(), RemoveHiCorrDims()
			  introduced dependency on qsar-class
				
December 2009 Minor fix for descriptor checking in load() 
			  Fix in get_dscr_pos()
			  altering load(): ability to load x-file alone, if mtxType=10
							 skipping empty lines in inside the matrix

November 2009 SVM format support, minor fix in get_NearNeibDistances()
August 8 2009 fixing Load()
June 16 2009 fixing Y-randomization
June 14 2009 minor fix for get_dscr() and get_sid() to return "" instead of INVALID for error.
March 19 2009 get_DimValues() fixed
March 18 2009 set_distMin/Max/Min0/Av() added to manage these vars + calc_dist_par() to recalc them.
March 5 2009 set_distMatr() and get_distMatr() implemented for saving/loading to/from external source
Feb 25 2009 add_dp() and remove_dp() implemented to manage individual datapoints
Jan 28 2009	all getdistance-f() modified to work with single datapoints, nonzero-mindist added.
Jun 16 2008	"==DUMMY_VAL" comparison was replaced by "< 0" check
May 25 2008	load() rewritten to tighten security when loading data
*/

//-------   memory leaks catcher for the current source-file  --------
#ifdef ADV_LEAK_CATCHER
#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif
//--------------------------------------------------------------------
#include "dataset.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
dataset::~dataset()
{
}

dataset::dataset()
{
	test.Dump();
	train.Dump();
	vars.Wipe();
	act.resize(0);
	patt.SetSize(0,0);
	pid.resize(0);
	sid.Wipe();
	sact.resize(0);
	dist.SetSize(0,0);
	distAv = distMax = distMin = distn0Min = DUMMY_VAL;
}

dataset::dataset(const dataset &dsInit)
: test(dsInit.test), train(dsInit.train), vars(dsInit.vars), act(dsInit.act), patt(dsInit.patt), pid(dsInit.pid), sid(dsInit.sid), sact(dsInit.sact), dist(dsInit.dist)
{//Copying the contents
	distAv	 = dsInit.distAv;
	distMax	 = dsInit.distMax;
	distMin	 = dsInit.distMin;
	distn0Min= dsInit.distn0Min;
}

void dataset::dump()
{
	dataset dtst;
	(* this) = dtst;
}

void dataset::sort_act()
{
	SIGNED_4B_TYPE N = act.length();
	SIGNED_4B_TYPE *AA = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N), ic = ZERO;
	REALNUM_TYPE *FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, N);
	QSortScore = FF;
	for (; ic < N; ic++)
	{
		AA[ic] = ic;	
		FF[ic] = act[ic]; 
	}
	qsort(AA, (size_t)N, sizeof(SIGNED_4B_TYPE), QSortCompareGreater);
	DROP_MEM_BLOCKS(FF);
	QSortScore = NULL;

	sact.resize(N);
	for (ic = 0; ic < N; ic++)	sact[ic] = AA[ic];
	DROP_MEM_BLOCKS(AA);
}


void dataset::set_MinNonZeroDistance(REALNUM_TYPE rtD)
{
	distn0Min = rtD;
}

void dataset::set_MinDistance(REALNUM_TYPE rtD)
{
	distMin = rtD;
}

void dataset::set_MaxDistance(REALNUM_TYPE rtD)
{
	distMax = rtD;
}

void dataset::set_AverageDistance(REALNUM_TYPE rtD)
{
	distAv = rtD;
}

void dataset::calc_dist_pars()
{
	UNSIGNED_4B_TYPE i, i1, N = patt.RowNO();
	distAv = distMax = distMin = distn0Min = DUMMY_VAL;
	
	if ((dist.RowNO() != N) || (!dist.IsSquare())) return;
	REALNUM_TYPE X;
	for (i = 0; i < N; i++)
	{		
		for (i1 = i + 1; i1 < N; i1++)
		{
			X = dist(i, i1);
			if (distMax < 0)
			{//checks if vars are not initialized, DUMMY_VAL const should be < 0, of course!
				distAv = distMax = distMin = distn0Min = X;
			}
			else
			{
				if (distMax < X) distMax = X;
				if (distMin > X) 
				{					
					distMin = X;
					if (X > 0)	distn0Min = X;
				}
				distAv += X;
			}
		}//for i1
	}//for i

	distAv *= 2;
	distAv /= N;
	distAv /= --N;
}

void dataset::calc_dist(REALNUM_TYPE rtDef, REALNUM_TYPE metricV, UNSIGNED_1B_TYPE metricKind)
//metricKind: 0 - Euclidean, 1 - Cosine-based, 2 - Correlation coff, 3 - Tanimoto
//all coefficients are turned into distances in a way of  = 1 - (coff)^metricV
{
	UNSIGNED_4B_TYPE i, i1, j, N = patt.RowNO(), D = patt.ColNO();
	dist.SetSize(N, N);
	apvector<REALNUM_TYPE> V1(D), V2(D);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < D; j++)	V1[j] = patt(i, j);
		dist(i,i) = rtDef;
		for (i1 = i + 1; i1 < N; i1++)
		{
			for (j = 0; j < D; j++)	V2[j] = patt(i1, j);
			dist(i, i1) = dist(i1, i) = getMetricDistance(V1, V2, metricV, metricKind);
		}//for i1
	}//for i

	calc_dist_pars();
}

void dataset::get_act_bins(apvector<UNSIGNED_4B_TYPE> &abins, UNSIGNED_4B_TYPE nbins, bool ifbsame)
{
	//abins.resize(0);
	if (nbins < 2) return;

	abins.resize(nbins);
	UNSIGNED_4B_TYPE i, N = act.length(), r = N / nbins;
	if (ifbsame)
	{
		for (i = 0; i < nbins; i++) abins[i] = r;
	}
	else
	{
		REALNUM_TYPE aBin, aRange; 
		aRange = get_MaxAct() - get_MinAct();
		aBin = aRange / nbins;
		REALNUM_TYPE aPrev = get_MinAct() - aRange, aNext = get_MinAct() + aBin;
		for (i = 0; i < nbins; i++)
		{
			aNext = (i == nbins-1) ? (get_MaxAct() + aRange) : (get_MinAct() + aBin*(i+1));
			abins[i] = get_ActPoints(aPrev, aNext).Size();
			aPrev = aNext;
		}		
	}

	//check abins
	for (r = i = 0; i < nbins; i++) r += abins[i];
	for (i = 0; i < nbins; i++)
	{
		if (r == N) break;
		if (r < N) { abins[i]++; r++; };
		if (r > N) { abins[i]--; r--; };
	}
}

//--------------------------------------------------------------
//algorithms of splitting data onto test and training sets
//--------------------------------------------------------------
void dataset::rand_split(REALNUM_TYPE rtFraction, UNSIGNED_1B_TYPE nActClasses, bool eqlSizeActBin, set *pExcl)
{
	rand_split( UNSIGNED_4B_TYPE (act.length() * rtFraction), nActClasses, eqlSizeActBin, pExcl);	
}

void dataset::rand_split(UNSIGNED_4B_TYPE nFraction, UNSIGNED_1B_TYPE nActClasses, bool eqlSizeActBin, set *pExcl)
{//random pick of nFraction; can keep activity distribution
	UNSIGNED_1B_TYPE nBins = max(nActClasses, 1);
	UNSIGNED_4B_TYPE Nmls = act.length();
	if (Nmls < 2) return;
	
	set setSeed(0, Nmls);
	apvector<SIGNED_4B_TYPE> spnts(sact);
	UNSIGNED_4B_TYPE i, r, added, currNmls;
	apvector<UNSIGNED_4B_TYPE> actNmls(nBins, Nmls / nBins), 
	testNmls(nBins, UNSIGNED_4B_TYPE(ceil(REALNUM_TYPE(nFraction) / nBins)) );
	if (nBins > 1) 
	{
		get_act_bins(actNmls, nBins, eqlSizeActBin);
		for (i = 0; i < nBins; i++)	testNmls[i] = (actNmls[i] * nFraction) / Nmls;		
		//check testNmls
		for (r = i = 0; i < nBins; i++) r += testNmls[i];
		for (i = 0; i < nBins; i++)
		{
			if (r == nFraction) break;
			if (r < nFraction) { testNmls[i]++; r++; };
			if (r > nFraction) { testNmls[i]--; r--; };
		}
	}
	else 
		setSeed.GetList(spnts); //no need to use ids sorted by activity, will be activity-independent.

	test.Dump();	
	for (currNmls = added = i = 0; i < nBins; i++)	
	if (testNmls[i])
	{
		Nmls = actNmls[i];
		set tNmls;
		while ( (actNmls[i] < Nmls + testNmls[i]) && (added < nFraction) )
		{//get randomly some portion
			r = currNmls + GetRandomNumber(Nmls);
			//check against the subset of points to be excluded
			if (pExcl)
			while (pExcl->IsInSet(spnts[r]) && (UNSIGNED_4B_TYPE(tNmls.Size()) < Nmls) )
			{
				tNmls.PutInSet(spnts[r]);
				r = currNmls + GetRandomNumber(Nmls);
			}

			if (UNSIGNED_4B_TYPE(tNmls.Size()) == Nmls) break; //to prevent a case when a bin is totally covered by pExcl

			test.PutInSet(spnts[r]);
			added++;
			spnts[r] = spnts[--Nmls + currNmls];
		}
		currNmls += actNmls[i];

		if (added == nFraction) break;
	}

	train = setSeed - test;
}

REALNUM_TYPE dataset::sfexcl_split(UNSIGNED_2B_TYPE Mode, REALNUM_TYPE f, bool ifForceDistCalc, UNSIGNED_1B_TYPE toTest, UNSIGNED_1B_TYPE toTrain, UNSIGNED_1B_TYPE nSeed, REALNUM_TYPE Kmetric)
{
	set nouserseed;
	return (sfexcl_split(nouserseed, Mode, f, ifForceDistCalc, nSeed));
}

REALNUM_TYPE dataset::sfexcl_split(set &seed, UNSIGNED_2B_TYPE Mode, REALNUM_TYPE f, bool ifForceDistCalc, UNSIGNED_1B_TYPE toTest, UNSIGNED_1B_TYPE toTrain, UNSIGNED_1B_TYPE nSeed, REALNUM_TYPE Kmetric)
//sphere exclusion implementation with some extra features added to the basic published version
//ref: Golbraikhm A, Tropsha A et al, J Comp-Aid Mol Design 2003, 17: 241-253
//There are two ways of calculating sphere's radius R:
//1) R =  f*(V/N)^1/k, where f is interpreted as Dissimilarity [0.2 .. 5.2]; 
//2) R =  mind.dist + f(max.dist- min.dist), where f should be varied [0 .. 0.25]
//
//Parameters:
//Mode				see the header-file
//seed				is the set of user-supplied points for seeding, nSeed - number of points to be seeded randomly in addition to that
//toTest, toTrain	#points from the sphere that are placed, respectively, into test and training set, 
//					which is done alternatingly (toTrain #points to training set, then toTest #points to test set, and repeated)
{
	UNSIGNED_4B_TYPE C, N = patt.RowNO(), D = patt.ColNO();
	SIGNED_4B_TYPE Z, Z1, el, el1;
	REALNUM_TYPE R, mnD, rtD;

	if ((N < 2) || (D < 2)) return 0;
	if ( (N != dist.ColNO()) || (!dist.IsSquare()) || ifForceDistCalc )
	{//!may be affected by changed metric-function
		if (Mode & SFEXCL_METRIC_COSINE)
			calc_dist(0, Kmetric, 1); //cosine-metric
		else
			if (Mode & SFEXCL_METRIC_CORR)
				calc_dist(0, Kmetric, 2); //similar to cosine-metric, but mean-centered
			else
				if (Mode & SFEXCL_METRIC_TANIMOTO)
					calc_dist(0, Kmetric, 3); //Tanimoto, though strange to apply it to non-discrete data : )
				else
					calc_dist(0, Kmetric, 0); //Euclidean-like
	}

	//----------------------
	//prepare seeding set
	if (nSeed)
	{
		if ((Mode & SFEXCL_SEED_BYACTS) == SFEXCL_SEED_BYACTS)	
			rand_split((UNSIGNED_4B_TYPE)nSeed, nSeed); 
		else 
			rand_split((UNSIGNED_4B_TYPE)nSeed, 1);
		seed |= test;
	}

	if ((Mode & SFEXCL_SEED_MINACT) == SFEXCL_SEED_MINACT)
	{
		C = 0;
		while (++C < N) if (act[sact[C]] > act[sact[0]]) break;
		Z = GetRandomNumber(C);
		seed.PutInSet(sact[Z]);
	}

	if ((Mode & SFEXCL_SEED_MAXACT) == SFEXCL_SEED_MAXACT)
	{
		C = 1;
		while (++C < N) if (act[sact[N-C]] < act[sact[N-1]]) break;
		Z = GetRandomNumber(--C) + 1;
		seed.PutInSet(sact[N - Z]);
	}

	if ((Mode & SFEXCL_R_BYDIST) == SFEXCL_R_BYDIST)
		R = f * (distAv - distMin) + distMin; //distMax is too big to use effectively!
	else
	{
		if ((Mode & SFEXCL_R_BYUSER) == SFEXCL_R_BYUSER) 
			R = f;
		else
		{//default
			REALNUM_TYPE l, h, V = 1;
			for (C = 0; C < D; C++)
			{//calculate volume of the descriptors' space
				patt.GetColScale(C, C, l, h);
				V *= pow(h - l, 1.0/D);
			}
			V /= pow(N, 1.0/D);
			R = f * V;
		}
	}

	if (R < distMin) R = distMin;
	
	//----------------------
	//prepare R-neiborhood-subsets to speed up splitting
	lneib Rneibs(N);
	for (C = 0; C < N - 1; C++)
	for (Z = C + 1; Z < SIGNED_4B_TYPE(N); Z++)
		if (dist(C, Z) < R)
		{
			Rneibs[C].PutInSet(Z);
			Rneibs[Z].PutInSet(C);
		}
	//----------------------
	test.Dump();
	train = seed;
	C = train.Size();
	SIGNED_4B_TYPE nspnts = 0;
	apvector<SIGNED_4B_TYPE> pnts_i(N), pnts(N - C), spnts(N), vecRneibs;
	for (Z1 = Z = 0; Z < SIGNED_4B_TYPE(N); Z++)
	{
		if (seed.IsInSet(Z)) continue;
		pnts[Z1++] = Z;
	}
	pnts.rand_shuffle(); //randomize datapoint positions
	for (Z = 0; Z < SIGNED_4B_TYPE(N - C); Z++)	pnts_i[pnts[Z]] = Z; //store randomized positions
	
	while (C < N)
	{//go on as long as there are points to exhaust

		if (seed.IsEmpty())
		{
			if ( ((Mode & SFEXCL_NEXTSF_RAND) == SFEXCL_NEXTSF_RAND) || (nspnts == 0) )
			//get random seeding point from the rest of the data
				Z = 0; //GetRandomNumber( N - C ); //array has been randomized already anyway
			else
			{//el would store -> pnts[]
				mnD = 0;
				if ((Mode & SFEXCL_NEXTSF_STEP2_MIN) == SFEXCL_NEXTSF_STEP2_MIN) mnD = distMax * N;
				if ((Mode & SFEXCL_NEXTSF_SPHERES) == SFEXCL_NEXTSF_SPHERES)
					for (el = Z = 0; Z + C < N; Z++)
					{
						rtD = 0;
						if ( ((Mode & SFEXCL_NEXTSF_STEP1_SUMDIST) != SFEXCL_NEXTSF_STEP1_SUMDIST) &&
							((Mode & SFEXCL_NEXTSF_STEP1_MIN) == SFEXCL_NEXTSF_STEP1_MIN) )
							rtD = distMax * N;

						for (Z1 = 0; Z1 < nspnts; Z1++)
						{
							if ((Mode & SFEXCL_NEXTSF_STEP1_SUMDIST) == SFEXCL_NEXTSF_STEP1_SUMDIST)
								rtD += dist(pnts[Z], spnts[Z1]);
							else
								if ( ((Mode & SFEXCL_NEXTSF_STEP1_MIN) == SFEXCL_NEXTSF_STEP1_MIN) ^ (dist(pnts[Z], spnts[Z1]) > rtD) )	
									rtD = dist(pnts[Z], spnts[Z1]);
						}

						if ( ((Mode & SFEXCL_NEXTSF_STEP2_MIN) == SFEXCL_NEXTSF_STEP2_MIN) ^ ( mnD < rtD ) )
						{
							mnD = rtD;
							el = Z;
						}
					}
				else //if ((Mode & SFEXCL_NEXTSF_SPHERES) == SFEXCL_NEXTSF_SPHERES)
					for (Z = 0; Z < nspnts; Z++)
					{
						el1 = 0;
						for (Z1 = 1; Z1 + C < N; Z1++)
						if ( ((Mode & SFEXCL_NEXTSF_STEP1_MIN) == SFEXCL_NEXTSF_STEP1_MIN) ^ 
							(dist(spnts[Z], pnts[Z1]) > dist(spnts[Z], pnts[el1])) )
							el1 = Z1;

						if ( ((Mode & SFEXCL_NEXTSF_STEP2_MIN) == SFEXCL_NEXTSF_STEP2_MIN) ^ ( mnD < dist(spnts[Z], pnts[el1]) ) )
						{
							mnD = dist(spnts[Z], pnts[el1]);
							el = el1;
						}
					}
				Z = el;	
			}//if ((Mode & SFEXCL_NEXTSF_RAND) == SFEXCL_NEXTSF_RAND) .. else

			el = pnts[Z];
			train.PutInSet(el);
			C++;
			pnts[Z] = pnts[N - C];
			pnts_i[el] = N;
			pnts_i[pnts[Z]] = Z;
		}//if (seed.IsEmpty())
		else
		{
			seed.GetElement(el);
			seed.RemoveFromSet(el);
		}

		spnts[nspnts++] = el;

		//now let's get the points within R-distance from el and distribute them between train and test
		Rneibs[el].GetList(vecRneibs);
		vecRneibs.rand_shuffle(); //to randomize the order of points!
		for (Z1 = el1 = 0; el1 < vecRneibs.length(); el1++)
		{
			Z = vecRneibs[el1];
			if (test.IsInSet(Z) || train.IsInSet(Z)) continue;
			C++;
			pnts[pnts_i[Z]] = pnts[N - C];
			pnts_i[pnts[N - C]] = pnts_i[Z];
			pnts_i[Z] = N;
			if ((Mode & SFEXCL_SPLIT_FIRST2TRN) == SFEXCL_SPLIT_FIRST2TRN)
			{
				if (Z1++ < toTrain)	train.PutInSet(Z); else	test.PutInSet(Z);
			}
			else
				if (Z1++ < toTest)	test.PutInSet(Z);	else	train.PutInSet(Z);
			if (Z1 == toTrain + toTest) Z1 = 0;
		}
	}//global-loop while (C < N)
	return (R);
}

void dataset::lgo_split(lneib &sLGO, REALNUM_TYPE rtLGO, UNSIGNED_1B_TYPE nActClasses, bool eqlSizeActBin)
{
	lgo_split(sLGO, (UNSIGNED_2B_TYPE)RRound(act.length() * rtLGO), nActClasses, eqlSizeActBin); //RRound added 03.30.2010
}

void dataset::lgo_split(lneib &sLGO, UNSIGNED_2B_TYPE nLGO, UNSIGNED_1B_TYPE nActClasses, bool eqlSizeActBin)
{	
	sLGO.resize(0);
	SIGNED_4B_TYPE Z, rep = Round(REALNUM_TYPE(act.length())/nLGO), REM;
	sLGO.resize(rep);
	REM = rep*nLGO - act.length();
	set setExclude;
	test.Dump();
	for (Z = 1; Z < rep; Z++)
	{
		setExclude |= test;
		UNSIGNED_4B_TYPE adj_nLGO = (UNSIGNED_4B_TYPE)nLGO;
		if (REM > 1) { adj_nLGO--; REM--; };
		if (REM < -1) { adj_nLGO++; REM++; };
		rand_split(adj_nLGO, nActClasses, eqlSizeActBin, &setExclude);
		sLGO[Z] = test;
	}
	sLGO[0] = train - setExclude;
}

void dataset::kohonen_split()
{//kohonen maps
}

void dataset::sca_split()
{//stochastic clustering algorithm
}


//retrieval
bool dataset::load(istream& i, istream *i1, UNSIGNED_1B_TYPE mtxType)
//mtxType = 0 (.xa file); = 1 (.x and .a files); = 2 (.svm file); = 10 (.x file only, no activities)
//NB: descriptor labels should not contain blanks or tabs!
//for extra-security, all reading in is done via string variable first
{
	if ((mtxType == 1) && (i1 == NULL)) return false; // if activity file required, generate error

	STRING_TYPE xch = " ", axch = "\t", sx = "  ", line;
	apvector<STRING_TYPE> lspl;
	char Bfr[20];
	STRING_TYPE stBfr;
	SIGNED_4B_TYPE N, D, n, d, q, im0, im2;

	line.getlinewithtabs(i);
	line.parse_string();
	line.replace(axch, xch);
	while (line.replace(sx, xch)); //treat consecutive delimiters as one

	SplitString(line, xch, lspl);
	if (lspl.length() < 2) return false;


	//---------------------------------------------------------------------------
	if (mtxType == 2)
	{//SVM format, easier to handle separately
		N = 256;
		D = lspl.length() - 1;
		patt.SetSize(N, D);
		act.resize(N);
		n = 0;
		while (!i.eof())
		{
			if (n == N)
			{
				N = n << 1;
				patt.SetSize(N, D);
				act.resize(N);
			}

			act[n] = atof(lspl[0].c_str());

			for (d = 1; d < lspl.length(); d++)
			{
				im0 = lspl[d].find(':');
				im2 = atoi(lspl[d].substr(0, im0).c_str());
				stBfr = lspl[d].substr(im0+1, lspl[d].length());
				if (im2 > D) 
				{ 
					D = im2;
					patt.SetSize(N, D);
				}
				patt(n, im2 - 1) = atof(stBfr.c_str());
			}
			n++;
			line.getlinewithtabs(i);
			line.parse_string();
			line.replace(axch, xch);
			while (line.replace(sx, xch)); //treat consequtive delimitors as one
			SplitString(line, xch, lspl);
			if (lspl.length() < 2) break;
		};//while (!i.eof())

		//only now true N and D are known - can generate SVM descriptor labels
		N = n;
		patt.SetSize(N, D);
		act.resize(N);
		sid.ids.resize(N);
		pid.resize(N);

		vars.Wipe();
		vars.L.resize(D);
		for (n = 0; n < N; n++)
		{
			sprintf(Bfr, "SVMDATA%d", n+1);
			pid[n] = Bfr;
			pid[n].parse_string();
			sid.ids[n] = pid[n];
			sid.AddRecord(n);
		}
		
		for (d = 0; d < D; d++)
		{
			sprintf (Bfr, "SVMDSCR%d", d+1);
			vars.L[d] = Bfr;
			vars.L[d].parse_string();
			vars.AddRecord(d);
		}

		sort_act();
		dist.SetSize(0, 0);
		return true;
	}//if (mtxType == 2)
	//---------------------------------------------------------------------------


	N = atoi(lspl[0].c_str());
	D = atoi(lspl[1].c_str());
	if ( (N < 1) || (D < 1) ) return false;
	dump();

	//descriptors
	vars.L.resize(D);
	line.getlinewithtabs(i); 
	line.parse_string();
	line.replace(axch, xch);
	while (line.replace(sx, xch));	
	
	SplitString(line, xch, lspl);
	if (lspl.length() < D)
	{//allows to read first D data-fields only, fixed Dec 3 2009
		dump();
		return false;
	}
	if (lspl.length() > D) //May 2011
		cout << "Warning, too many descriptor labels: " << lspl.length() << endl;

	for (d = 0; d < D; d++)
	{
		vars.L[d] = lspl[d];
		vars.L[d].parse_string();
		vars.AddRecord(d);
	}

	patt.SetSize(N, D);
	sid.ids.resize(N);
	pid.resize(N);
	//3 extra dimensions are seq#, ID and activity (only in a new format)
	act.resize(N);
	if (mtxType == 0) q = 3; else q = 2;
	for (n = 0; (n < N) && !i.eof(); n++) 
	{
		do 
		{
			line.getlinewithtabs(i); 
			line.parse_string();
		} while (line.length() == 0); //to skip empty lines within the matrix file

		line.replace(axch, xch);
		while (line.replace(sx, xch));
		SplitString(line, xch, lspl);		
		if (lspl.length() < (D + q)) break;	

		if (lspl.length() > (D + q)) //May 2011
			cout << "Line#" << n+1 << " has too many records: " << lspl[0] << " " << lspl[1] << " " << lspl[2] << " ..." << endl;

		pid[n] = lspl[0];
		pid[n].parse_string();
		sid.ids[n] = lspl[1];
		if (mtxType == 0)
			act[n] = atof(lspl[2].c_str());
		else
		{
			if ((i1 == NULL) || (mtxType == 10)) // if x-file only, then make up fake activities
				act[n] = 0;
			else
			{//can check for match with any of id columns in x-file
				(*i1) >> sid.ids[n]; 
				(*i1) >> act[n];
			}
		}

		sid.ids[n].parse_string();
		sid.AddRecord(n);
		for (d = 0; d < D; d++)	patt(n, d) = atof(lspl[d + q].c_str());
	}

	if (n < N) 
	{
		dump();
		return false;
	}

	if (!i.eof())
	{//read normalising data directly
		vars.B.resize(D);
		for (d = 0; (d < D) && !i.eof(); d++) i >> vars.B[d];
		if (d == D)
		{
			vars.A.resize(D);
			for (d = 0; (d < D) && !i.eof(); d++) { i >> vars.A[d]; vars.A[d] -= vars.B[d]; };
			if (d < D) 
			{
				vars.A.resize(0);
				vars.B.resize(0);
			}
		}
		else
			vars.B.resize(0);
	}

	//initialize
	sort_act();
	dist.SetSize(0, 0);

	return true;
}

bool dataset::save(ostream& o, ostream *o1, UNSIGNED_1B_TYPE mtxType)
{
	if ((mtxType ==1) && (o1 == NULL)) return false;

	SIGNED_4B_TYPE N = patt.RowNO(), D = patt.ColNO(), n, d;
	
	char Bfr[100], xch = ' ';
	STRING_TYPE stBfr;

	if (mtxType ==2)
	{//SVM format, quite simple (no descriptor or data labels, no scaling info. NB: only 3 decimal digits are saved!
		for (n = 0; n < N; n++)
		{
			sprintf(Bfr, "%12.3f", act[n]);
			stBfr = Bfr;
			stBfr.parse_string();
			o << stBfr;

			for (d = 0; d < D; d++)
			{
				if (fabs(patt(n, d)) < 0.001) continue; //3 digits! see below sprintf(.. "%12.3f" .. )

				sprintf(Bfr, "%12.3f", patt(n, d));
				stBfr = Bfr; stBfr.parse_string();
				sprintf(Bfr, "%d:%s", (d + 1), stBfr.c_str());
				stBfr = Bfr; stBfr.parse_string();
				o << xch << stBfr;
			}//for d
			o << endl; // each line should have line end in LIBSVM format!
		}//for n

		return true;
	}

	o << N << xch << D << endl;

	o << vars.L[0];
	for (d = 1; d < D; d++) o << xch << vars.L[d];
	o << endl;

	for (n = 0; n < N; n++) 
	{
		o << pid[n] << xch << sid.ids[n]; //seq#, id

		sprintf(Bfr, "%12.6f", act[n]);
		stBfr = Bfr;
		stBfr.parse_string();

		if (mtxType == 0) 
			o << xch << stBfr; //act
		else
			(*o1) << sid.ids[n] << xch << stBfr << endl; //act

		for (d = 0; d < D; d++)  
		{
			sprintf(Bfr, "%12.6f", patt(n, d));
			stBfr = Bfr;
			stBfr.parse_string();

			o << xch << stBfr;
		}
		o << endl;
	}

	if (vars.B.length() != D) return true;

	for (d = 0; d < D; d++) 
	{
		sprintf(Bfr, "%12.6f", vars.B[d]);
		stBfr = Bfr;
		stBfr.parse_string();
		if (d + 1 < D) stBfr += xch;
		o << stBfr;
	}
	o << endl;

	for (d = 0; d < D; d++) 
	{
		sprintf(Bfr, "%12.6f", vars.A[d] + vars.B[d]);
		stBfr = Bfr;
		stBfr.parse_string();
		if (d + 1 < D) stBfr += xch;
		o << stBfr;
	}
	o << endl;

	return true;
}

dataset dataset::subset(set &Subset)
{
	set Sub = Subset; //otherwise Subset is changed
	dataset dnew;
	SIGNED_4B_TYPE n = 0, N = Sub.Size(), d, D = patt.ColNO(), e;

	dnew.vars = vars;
	dnew.act.resize(N);
	dnew.patt.SetSize(N, D);
	
	dnew.sid.ids.resize(N);
	dnew.pid.resize(N);
	while (n < N)
	{
		Sub.GetElement(e);
		Sub.RemoveFromSet(e);
		dnew.act[n]	 = act[e];

		dnew.pid[n] = pid[e];
		dnew.sid.ids[n]	 = sid.ids[e];
		dnew.sid.AddRecord(n);

		for (d = 0; d < D; d++)	dnew.patt(n, d) = patt(e, d);
		n++;
	}

	dnew.sort_act();
	dnew.dist.SetSize(0, 0);
	dnew.distAv = dnew.distMax = dnew.distMin = DUMMY_VAL;

	return dnew;
}

dataset dataset::get_training_set()
{
	return subset(train);
}

dataset dataset::get_test_set()
{
	return subset(test);
}


STRING_TYPE dataset::get_sid(SIGNED_4B_TYPE posN)
//get datapoint's id-name from its sequential number
{
	if ( (posN < 0) || posN > sid.ids.length() ) return "";
	return sid.ids[posN];
}

STRING_TYPE dataset::get_dscr(SIGNED_4B_TYPE posN)
//get descriptor's id-name from its sequential number
{
	if ( (posN < 0) || posN > vars.L.length() ) return "";
	return vars.L[posN];
}

SIGNED_4B_TYPE dataset::get_sid_pos(STRING_TYPE stID)
//get datapoint's sequential number from its id-name
{
	UNSIGNED_4B_TYPE n, e = UNSIGNED_4B_TYPE(sid.ids.length()), a, m = 0;
	sid.ids.resize(e + 1);
	sid.ids[e] = stID;
	a = sid.FindHashKey( sid.hashFunction(e) );
	n = sid.RetrieveRecord(a, m);
	while (n--)
	{
		if (m <= e) //in case of an obsolete hash-record for a deleted entry
		if (sid.ids[m].length() == stID.length())
		if (sid.ids[m] == stID)
		{
			sid.ids.resize(e);
			return m;
		}
		sid.RetrieveRecord(a, m);
	}
	sid.ids.resize(e);
	return INVALID;
}

SIGNED_4B_TYPE dataset::get_dscr_pos(STRING_TYPE stID)
//get descriptor's sequential number from its unique name
{
	UNSIGNED_4B_TYPE n, e = UNSIGNED_4B_TYPE(vars.L.length()), a, m = 0;
	vars.L.resize(e + 1);
	vars.L[e] = stID;
	a = vars.FindHashKey( vars.hashFunction(e) );
	n = vars.RetrieveRecord(a, m);
	while (n--)
	{
		if (vars.L[m].length() == stID.length())
		if (vars.L[m] == stID)
		{
			vars.L.resize(e);
			return m;
		}
		vars.RetrieveRecord(a, m);
	}

	vars.L.resize(e);
	return INVALID;
}

REALNUM_TYPE dataset::get_MinNonZeroDistance(SIGNED_4B_TYPE dp)
//returns adjusted minimum distance in the distance-matrix dist()
//or DUMMY_VAL if dist-matrix was incorrect
//(overlapping points with 0-distance are excluded from calculation)
{
	if (dp == INVALID)	return distn0Min;

	SIGNED_4B_TYPE i, N = act.length();
	if ( (SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()) ) return DUMMY_VAL;
	REALNUM_TYPE rtV = distMax;
	for (i = 0; i < N; i++)
	{
		if (i == dp) continue;
		if (dist(dp, i) < rtV) rtV = dist(dp, i);
	}
	return rtV;
}

REALNUM_TYPE dataset::get_MinDistance(SIGNED_4B_TYPE dp)
//returns minimum distance in the distance-matrix dist() 
//or DUMMY_VAL if dist-matrix was incorrect
//NB: min.dist will be 0 if there are overlapping points!
{
	if (dp == INVALID)	return distMin;

	SIGNED_4B_TYPE i, N = act.length();
	if ( (SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()) ) return DUMMY_VAL;
	REALNUM_TYPE rtV = distMax;
	for (i = 0; i < N; i++)
	{
		if (i == dp) continue;
		if (dist(dp, i) < rtV) rtV = dist(dp, i);
	}
	return rtV;
}

REALNUM_TYPE dataset::get_MaxDistance(SIGNED_4B_TYPE dp)
//returns maximum distance in the distance-matrix dist()
//or DUMMY_VAL if dist-matrix was incorrect
{
	if (dp == INVALID)	return distMax;
	
	SIGNED_4B_TYPE i, N = act.length();
	if ( (SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()) ) return DUMMY_VAL;
	REALNUM_TYPE rtV = 0;
	for (i = 0; i < N; i++)
	{
		if (i == dp) continue;
		if (dist(dp, i) > rtV) rtV = dist(dp, i);
	}
	return rtV;
}

REALNUM_TYPE dataset::get_AverageDistance(SIGNED_4B_TYPE dp, REALNUM_TYPE cutoff)
//returns average distance in the distance-matrix dist() for a given cutoff (max.dist to use)
//or DUMMY_VAL if dist-matrix was incorrect
{
	SIGNED_4B_TYPE i, j, N = act.length();
	if ( (SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()) ) return DUMMY_VAL;
	bool noCutoff = ((cutoff < 0) || (cutoff > distMax) || (cutoff < distMin));
	UNSIGNED_4B_TYPE n = 0;
	REALNUM_TYPE rtV = 0;

	if (dp == INVALID)
	{//scan whole matrix;
		if (noCutoff) return distAv;
		
		for (i = 0; i < N - 1; i++)
		for (j = i + 1; j < N; j++)
		{
			if (dist(i, j) > cutoff) continue;
			rtV += dist(i, j);
			n++;
		}
	}
	else
	{
		for (i = 0; i < N; i++)
		{
			if (i == dp) continue;
			if (!noCutoff)	if (dist(i, dp) > cutoff) continue;
			rtV += dist(i, dp);
			n++;
		}
	}

	if (n)	rtV /= n;
	return rtV;
}

REALNUM_TYPE dataset::get_Distance(SIGNED_4B_TYPE dp1, SIGNED_4B_TYPE dp2)
{
	if ( (dp1 < 0) || (dp2 < 0) ) return DUMMY_VAL;
	if ( SIGNED_4B_TYPE(min(dist.RowNO(), dist.ColNO())) <  max(dp1, dp2) ) return DUMMY_VAL;
	
	return dist(dp1, dp2);
}

REALNUM_TYPE dataset::get_indDistance(SIGNED_4B_TYPE dp1, SIGNED_4B_TYPE dp2, set * pDims, REALNUM_TYPE metricV, UNSIGNED_1B_TYPE metricKind)
//calculates distance between datapoints dp1 and dp2, based on pDims descriptors
{
	//if ( (!dist.IsSquare()) || (dist.RowNO() < patt.RowNO()) )		return DUMMY_VAL;

	SIGNED_4B_TYPE d, D = patt.ColNO();
	apvector<SIGNED_4B_TYPE> lDims;
	if (pDims == NULL) 
	{
		lDims.resize(D);
		for (d = 0; d < D; d++) lDims[d] = d;
	}
	else
	{
		pDims->GetList(lDims);
		D = lDims.length();
	}

	apvector<REALNUM_TYPE> V1(D), V2(D);

	for (d = 0; d < D; d++)
	{
		V1[d] = patt(dp1, lDims[d]);
		V2[d] = patt(dp2, lDims[d]);
	}

	return getMetricDistance(V1, V2, metricV, metricKind);
}

void dataset::get_NearNeibDistances(apvector<REALNUM_TYPE> &Stats, UNSIGNED_2B_TYPE kNeibours, REALNUM_TYPE cutoffR, UNSIGNED_1B_TYPE kmode)
//calculates av, max, stdev for kNeibours-number throughout the dataset
//NB: if ifSqrDist is true, then average and the stdev of the squared distance is returned!
{
	bool ifSqrDist = ((kmode & 1) == 1);
	bool ifKstrict = ((kmode & 2) == 2);

	SIGNED_4B_TYPE i, N = act.length();

	Stats.resize(0);
	if ((SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()))
		return;

	REALNUM_TYPE v = 0, mxkd = distMin, Spread;
	apvector<REALNUM_TYPE> veckNNdists(N * max(kNeibours, 1), 0);
	apvector<SIGNED_4B_TYPE> NNlist;

	UNSIGNED_4B_TYPE j, ln, n = 0;
	for (i = 0; i < N; i++)
	{
		get_NearNeib(i, NNlist, kNeibours, cutoffR);		
		
		if (ifKstrict) //so as not to copy more than kNeibours-points
			ln = min(kNeibours, NNlist.length());
		else
			ln = NNlist.length();

		j = (UNSIGNED_4B_TYPE) veckNNdists.length();
		if (j < n + ln)
			veckNNdists.resize(j + ln + 1);

		for (j = 0; j < ln; j++)
		{//sum up distances of only k neighbours
			veckNNdists[j + n] = dist(i, NNlist[j]);
			if (ifSqrDist)	veckNNdists[j + n]	*= dist(i, NNlist[j]);
			v += veckNNdists[j + n];
		}
		n += ln;

		if ( mxkd < dist(i, NNlist[ln-1]) ) mxkd = dist(i, NNlist[ln-1]);
	}//for i < N

	v /= n; //needed average
	
	//now calculate the spread
	Spread = 0;
	for (j = 0; j < n; j++)	Spread += sqr(veckNNdists[j] - v);
	Spread /= n-1;
	Spread = sqrt(Spread);

	Stats.resize(3);
	Stats[0] = v;
	Stats[1] = mxkd;
	Stats[2] = Spread;
}

void dataset::get_NearNeib(SIGNED_4B_TYPE datap, apvector<SIGNED_4B_TYPE> &NNeibs, UNSIGNED_2B_TYPE kNeibours, REALNUM_TYPE cutoffR)
//returns array of nearest-neighbours 
//sorted by their distance to datap from smallest to biggest
//precondition:	dist() should be filled, distMax should be valid!
{
	SIGNED_4B_TYPE j, N = act.length();

	NNeibs.resize(0);
	if ( (SIGNED_4B_TYPE(dist.ColNO()) != N) || (!dist.IsSquare()) )
		return;

	SIGNED_4B_TYPE *AA = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N);
	REALNUM_TYPE *FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, N);
	QSortScore = FF;

	for (j = 0; j < N; j++)
	{
		if (j == datap) //can not just continue, indexing will be broken!
			FF[j] = 1 + distMax;
		else
			FF[j] = dist(datap, j);
		AA[j] = j;
	}

	qsort(AA, (size_t)N, sizeof(SIGNED_4B_TYPE), QSortCompareGreater);
	
	REALNUM_TYPE Rc = cutoffR;
	if (kNeibours > 0)
	{//use cut-off from last k-th neighbour
		j = kNeibours;
		Rc = FF[ AA[j - 1] ];
	}
	else //use a radial cut-off supplied by the user
		j = 0;
	
	for (; j < N; j++)
	if (FF[ AA[j] ] > Rc ) //would recheck one neighbour, the last one
		break;

	NNeibs.resize(j);
	while (j-- > 0) //save nearest neighbours
		NNeibs[j] = AA[j];

	QSortScore = NULL;
	DROP_MEM_BLOCKS(FF);
	DROP_MEM_BLOCKS(AA);
}

void dataset::randomizeY(UNSIGNED_1B_TYPE nActClasses, bool eqlSizeActBin)
{
	UNSIGNED_4B_TYPE N = act.length(), nB =  max(nActClasses, 1);

	if (N < nActClasses) return;

	apvector<UNSIGNED_4B_TYPE> abins(nB, N/nB);
	if (nB > 1) get_act_bins(abins, nB, eqlSizeActBin);
	
	UNSIGNED_4B_TYPE r, i, ca;
	for (ca = i = 0; i < nB; i++)	
	{
		apvector<REALNUM_TYPE> ract(abins[i]);
		for (r = 0; r < abins[i]; r++)	ract[r] = act[sact[r + ca]];
		//shuffle randomly
		SortRandomly(ract); //more accurate than ract.rand_shuffle();
		for (r = 0; r < abins[i]; r++)	act[sact[r + ca]] = ract[r];
		ca += abins[i];
	}
	sort_act(); //update sact-array
}

void dataset::reduce_dimensions(set &Dims)
{//retains only descriptors that are in Dims
	if (Dims.IsEmpty())  return;
	
	apvector<SIGNED_4B_TYPE> vecDims;
	Dims.GetList(vecDims);
	SIGNED_4B_TYPE d, n, D = vecDims.length(), N = patt.RowNO();
	if ((SIGNED_4B_TYPE)patt.ColNO() == D) return;

	dataset datasetSub(*this);
	
	datasetSub.dist.SetSize(0,0);
	datasetSub.distMax = distMin = distAv = DUMMY_VAL;
	datasetSub.vars.Wipe();
	datasetSub.patt.SetSize(N, D);
	if (vars.B.length() > D)
	{
		datasetSub.vars.A.resize(D);
		datasetSub.vars.B.resize(D);
	}

	for (d = 0; d < D; d++)
	{
		datasetSub.vars.L.resize(D);
		datasetSub.vars.L[d] = vars.L[vecDims[d]];
		datasetSub.vars.AddRecord(d);
		if (vars.B.length() > D)
		{
			datasetSub.vars.A[d] = vars.A[vecDims[d]];
			datasetSub.vars.B[d] = vars.B[vecDims[d]];
		}

		for (n = 0; n < N; n++)	datasetSub.patt(n, d) = patt(n, vecDims[d]);
	}

	this->dump();
	(*this) = datasetSub;
}

REALNUM_TYPE dataset::get_OccupVol(REALNUM_TYPE Probe, UNSIGNED_4B_TYPE EstimationMethod, bool ifForceDistCalc, REALNUM_TYPE metricV, UNSIGNED_1B_TYPE metricKind)
{
	SIGNED_4B_TYPE d, n, D = patt.ColNO(), N = patt.RowNO();
	if ( (N != SIGNED_4B_TYPE(dist.ColNO())) || (!dist.IsSquare()) || ifForceDistCalc )	calc_dist(0, metricV, metricKind);

	REALNUM_TYPE v = 0, Probe2 = 2*Probe;
	UNSIGNED_4B_TYPE nRefpnts = 0, D2 = D << 1;


	//prepare R-neiborhood-subsets (only the low-id spheres will know about their high-id neighbors)
	apvector<SIGNED_4B_TYPE> vecRneib;
	lneib Rneibs(N);
	SIGNED_4B_TYPE z, d1;
	for (n = 0; n < N - 1; n++)
	for (z = n + 1; z < N; z++)
	if (dist(n, z) < Probe2)	
	{
		Rneibs[n].PutInSet(z);
		if (EstimationMethod) Rneibs[z].PutInSet(n);
	}
	//----------------------

	if (EstimationMethod == 1)
	{//the number of points without neighbours for the specified R
		for (n = 0; n < N; n++)	if (Rneibs[n].IsEmpty()) nRefpnts++;
		v = nRefpnts;
	}

	if (EstimationMethod == 2)
	{//the average number of neighbours for the specified R
		for (n = 0; n < N; n++)	nRefpnts += Rneibs[n].Size();
		v = nRefpnts;
		v /= N;
	}

	if (EstimationMethod == 0)
	{//by spheres and ref-points

		apvector<REALNUM_TYPE> refpnt(D), neibpnt(D);
		UNSIGNED_1B_TYPE ref2;
		for (n = 0; n < N; n++)
		{
			nRefpnts += D2;
			if (Rneibs[n].IsEmpty())	continue;
			for (d = 0; d < D; d++)	refpnt[d] = patt(n, d);

			Rneibs[n].GetList(vecRneib);
			for (d1 = 0; d1 < D; d1++)
			{
				ref2 = 2;
				refpnt[d1] += Probe;
				do
				{
					for (z = 0; z < vecRneib.length(); z++)
					{
						for (d = 0; d < D; d++)	neibpnt[d] = patt(vecRneib[z], d);
						if (Probe > getMetricDistance(refpnt, neibpnt, metricV, metricKind))	
							break;
					}
					if (z < vecRneib.length())	
						nRefpnts--;
					refpnt[d1] -= Probe2;
				} while (--ref2);
				refpnt[d1] += Probe2 + Probe;
			}//d1
		}//n

		//calculate volume as Probe's spherical volume  * #ref-points/ (2 * Dimensions);
		v = nRefpnts;

		/*//report volume relative to n-hypersphere's volume. (otherwise absolute volume is too mall for high D)!

		//Recursively calculate the volume of the D-hypersphere: Vn = Vn-2 * 2piR^2/n; V0 = 1; V1 = 2R; V2 = piR^2;
		D2 = D >> 1;
		ref2 = 0;
		if ((D2 << 1) == D )  v *= 0.5; else { v = Probe; ref2 = 1; };
		Probe2 *= pi * Probe; //2piR^2		
		for (d = 1; d <= D2; d++)	
		{
			v *= Probe2;
			v /= REALNUM_TYPE((d << 1) + ref2);
		}//*/

		v /= D2 ;//D;
	}

	return v;
}


void dataset::get_DimValues(SIGNED_4B_TYPE D, apvector<REALNUM_TYPE> &Dim)
//returns array of descriptors' values for a datapoint D
{
	SIGNED_4B_TYPE PN = patt.RowNO(), PD = patt.ColNO();
	if ( (D >= PN) || (D < 0) ) return;
	Dim.resize(PD);
	for (SIGNED_4B_TYPE i = 0; i < PD; i++)	Dim[i] = patt(D, i);
}

void dataset::get_DimRowValues(SIGNED_4B_TYPE X, apvector<REALNUM_TYPE> &Dt)
//returns array of values for a descriptor X
{
	SIGNED_4B_TYPE PN = patt.RowNO(), PD = patt.ColNO();
	if ( (X >= PD) || (X < 0) ) return;
	Dt.resize(PN);
	for (SIGNED_4B_TYPE i = 0; i < PN; i++)	Dt[i] = patt(i, X);
}

void dataset::get_ActValues(apvector<REALNUM_TYPE> &Acts)
{
	Acts = act;
}

void dataset::set_ActValues(apvector<REALNUM_TYPE> &Acts)
{
	if ( Acts.length() == act.length() )
	{
		act = Acts;
		sort_act();
	}
}

void dataset::set_Act(SIGNED_4B_TYPE dp1, REALNUM_TYPE val)
{
	if ( (act.length() > dp1) && (dp1 >= 0) )
	{
		act[dp1] = val;
		sort_act();
	}
}

REALNUM_TYPE dataset::get_Act(SIGNED_4B_TYPE dp1)
{
	if ( (act.length() > dp1) && (dp1 >= 0) )
		return (act[dp1]);
	return DUMMY_VAL;
}

REALNUM_TYPE dataset::get_MaxAct()
{
	if (sact.length())
		return (act[sact[sact.length() - 1]]);
	return DUMMY_VAL;
}

REALNUM_TYPE dataset::get_MinAct()
{
	if (sact.length())
		return (act[sact[0]]);
	return DUMMY_VAL;
}

REALNUM_TYPE dataset::get_AverageAct()
{
	SIGNED_4B_TYPE i = 0, N = act.length();
	REALNUM_TYPE v = 0;
	for (; i<N; i++)	v += act[i];
	v /= N;
	return (v);
}

set dataset::get_ActPoints(REALNUM_TYPE lowAct, REALNUM_TYPE hiAct)
//returns all points with lowAct < Act <= hiAct
{
	set setP;
	SIGNED_4B_TYPE i = 0, N = sact.length();

	for (; i<N; i++) 
	{
		if (lowAct > act[sact[i]]) continue;
		if (hiAct < act[sact[i]]) break;
		setP.PutInSet(sact[i]);
	}	
	return (setP);
}

SIGNED_4B_TYPE dataset::get_Ndatapoints()
{
	return SIGNED_4B_TYPE(patt.RowNO());
}

SIGNED_4B_TYPE dataset::get_Ndimensions()
{
	return SIGNED_4B_TYPE(patt.ColNO());
}

bool dataset::expandby(dataset &dtstA, bool UseAsTest)
/*	description:	expands current dataset by data from dtstA;
	precondition:	dimensions should match
	postcondition:	dist-matrix is resized to 0; 
					NB: datapoints are copied without check of duplicate ids!
						but dimensions are checked for matching ids!

					returns false if dimensions do not match
*/				
{
	if (dtstA.get_Ndimensions() != get_Ndimensions()) return false;
	
	SIGNED_4B_TYPE i, j, nA = dtstA.get_Ndatapoints(), n = get_Ndatapoints(), d = get_Ndimensions();

	apvector<SIGNED_4B_TYPE> di(d);
	for (j = 0; j < d; j++)
	{
		di[j] = get_dscr_pos(dtstA.vars.L[j]);
		if (di[j] == INVALID) return false;
	}

	patt.SetSize(n + nA, d);
	sid.ids.resize(n + nA);
	pid.resize(n + nA);	
	act.resize(n + nA);
	set newTest;
	
	for (i = 0; i < nA; i++)
	{		
		act[n + i] = dtstA.act[i];
		pid[n + i] = dtstA.pid[i];
		sid.ids[n + i] = dtstA.sid.ids[i];
		sid.AddRecord(n + i);
		for (j = 0; j < d; j++)
		{
			patt(n + i, di[j]) = dtstA.patt(i, j);
			if (UseAsTest) newTest.PutInSet(n + i);
		}
	}

	if (UseAsTest)	test = newTest;

	sort_act();
	dist.SetSize(0, 0);

	return true;
}

void dataset::add_dp(apvector<REALNUM_TYPE> &dscr, REALNUM_TYPE actdp, STRING_TYPE sidp, STRING_TYPE pidp)
{
	SIGNED_4B_TYPE dp = get_Ndatapoints(), dmn = get_Ndimensions();
	if (dmn != dscr.length()) return;
	dp++;
	patt.SetSize(dp, dmn);
	act.resize(dp);
	sact.resize(dp);
	pid.resize(dp);
	sid.ids.resize(dp);
	dp--;

	SIGNED_4B_TYPE i = 0;
	for (; i < dmn; i++) patt(dp, i) = dscr[i];
	
	act[dp] = actdp;
	sact[dp]= dp;
	
	SIGNED_4B_TYPE nSwap = 0;
	for (i = 0; i < dp; i++) 
	if (act[sact[i]] > act[sact[dp]]) Swap(sact[dp], sact[i], nSwap);

	char bb[20];
	sprintf(bb, "#%d#", dp);
	if (pidp.length())	pid[dp]	= pidp;	else pid[dp] = bb;
	if (sidp.length())	sid.ids[dp]	= sidp;	else sid.ids[dp] = bb;
	sid.AddRecord(dp);
}

void dataset::remove_dp(SIGNED_4B_TYPE dp)
{
	SIGNED_4B_TYPE i, dpn = get_Ndatapoints(), dmn = get_Ndimensions();
	
	dpn--;
	for (i = 0; i < dpn; i++)	if (sact[i] == dp) break;
	while (i < dpn) { sact[i] = sact[i+1]; i++; };	
	if (dp < dpn)
	{
		for (i = 0; i < dmn; i++)	patt(dp, i) = patt(dpn, i);
		act[dp] = act[dpn];
		pid[dp]	= pid[dpn];
		sid.ids[dp]	= sid.ids[dpn];
		sid.AddRecord(dp); //add updated hash-key
	}

	patt.SetSize(dpn, dmn);
	act.resize(dpn);
	sact.resize(dpn);
	pid.resize(dpn);
	sid.ids.resize(dpn);
}


void dataset::get_distMatr(matrix<REALNUM_TYPE> &DST)
{
	DST = dist;
}

bool dataset::set_distMatr(matrix<REALNUM_TYPE> &DST)
{
	if (DST.IsSquare())
		dist = DST;
	else return false;
	return true;
}



void dataset::normalizeby(dataset &RF)
{
	UNSIGNED_4B_TYPE dz, dnAll = get_Ndimensions(), dn = RF.get_Ndimensions(), dm = get_Ndatapoints();
	UNSIGNED_4B_TYPE im0, im1, im2;	

	if ( isScaled() ) scale_dimensions(3);	//denormalize back to the original values, if scaling exists! Aug 11 2011

	RF.vars.L.resize(dn + 1);
	apvector<SIGNED_4B_TYPE> dorder(dnAll, INVALID);
	for (dz = 0; dz < dnAll; dz++)
	{
		RF.vars.L[dn] = vars.L[dz];
		im1 = RF.vars.FindHashKey( RF.vars.hashFunction(dn) );
		im0 = RF.vars.RetrieveRecord(im1, im2);
		while (im0 > 0)
		{
			if ( vars.L[dz].length() == RF.vars.L[im2].length() )
			if ( vars.L[dz] == RF.vars.L[im2] )
			{
				dorder[dz] = im2;
				break;
			}			
			RF.vars.RetrieveRecord(im1, im2);
			im0--;
		}
	}//for dz

	RF.vars.L.resize(dn); //size back
	vars	= RF.vars;

	set	setPM; //dimensions to keep
	//resort the matrix to comply with specified order in the file!
	matrix<REALNUM_TYPE>  newpatt(dm, dn);
	
	for (dz = 0; dz < dnAll; dz++)
	{
		if (dorder[dz] == INVALID) continue;
		setPM.PutInSet(dorder[dz]);
		for (im0 = 0; im0 < dm; im0++) newpatt(im0, dorder[dz]) =  patt(im0, dz);
	}
	patt = newpatt;
	reduce_dimensions(setPM);
	if ( isScaled() )	scale_dimensions(2);
}

void dataset::scale_dimensions(UNSIGNED_1B_TYPE scl_mode)
/* 		scl_mode: 
(ONLY IF there are NO normalizing coefficients present)
		0 - range-scaling
		1 - autoscaling
(ONLY IF there ARE normalizing coefficients)
		2 - renormalize using existing coefficients; 
		3 - denormalize back to original values
*/
{
	UNSIGNED_4B_TYPE n = get_Ndimensions(), PN = get_Ndatapoints();
	bool NoCoffs = true;

	UNSIGNED_4B_TYPE f, u;
	apvector<REALNUM_TYPE> apParCol(PN);

	if (vars.A.length() == (SIGNED_4B_TYPE)n)
	{
		NoCoffs = false;
		if (scl_mode < 2) return; //do nothing
	}
	else
	{
		if (scl_mode > 1) return; //do nothing

		vars.A.resize(n);
		vars.B.resize(n);
		for (u = 0; u < n; u++)
		{ //clean up
			vars.A[u] = 1;
			vars.B[u] = 0;
		};
	}

	QSAR qq;
	//rescale parameters
	for (u = 0; u < n; u++)
	{
		if (NoCoffs)
		{
			for (f = 0; f < PN; f++)	apParCol[f] = patt(f, u);
			if (scl_mode == 1)
			{//Autoscale			
				vars.B[u] = qq.meanV(apParCol);
				vars.A[u] = qq.stdev(apParCol);
			}
			else
			{//range-scale
				vars.B[u]	 = qq.minV(apParCol);
				vars.A[u]	 = qq.maxV(apParCol) - vars.B[u];			
			}
		}

		if (vars.A[u] < SMALL_NUMBER)	//warning, very small coefficient
			vars.A[u] = INVALID;

		//rescale
		for (f = 0; f < PN; f++)
		if (scl_mode == 3)
		{//denormalize back to original values
			patt(f, u)	*= vars.A[u];
			patt(f, u)	+= vars.B[u];
		}
		else
		{//normalize
			patt(f, u)	-= vars.B[u];
			patt(f, u)	/= vars.A[u];
		}
	}//for u

	if (scl_mode == 3)
	{//remove coefficients, if matrix was restored to original values
		vars.A.resize(0);
		vars.B.resize(0);
	}
}

void dataset::RemoveHiCorrDims(REALNUM_TYPE MaxCorr)
{//removes dimensions whose squared correlation coefficient is higher than MaxCorr

	UNSIGNED_4B_TYPE f, u, u1, np = get_Ndatapoints(), nd = get_Ndimensions();
	apvector<REALNUM_TYPE> dv1(np), dv2(np);
	QSAR qq;
	REALNUM_TYPE r;
	set Dims(0, nd);
	for (u = 0; u < nd - 1; u++)
	{
		for (f = 0; f < np; f++)	dv1[f] = patt(f, u);
		for (u1 = u + 1; u1 < nd; u1++)
		{
			for (f = 0; f < np; f++)	dv2[f] = patt(f, u1);
			
			r = qq.correl(dv1, dv2);
			if (r * r > MaxCorr)
			{//highly correlated descriptors were found
				Dims.RemoveFromSet(u);
				break;
			}
		}//for u1
	}//for u
	reduce_dimensions(Dims);
}

void dataset::RemoveLowVarDims(REALNUM_TYPE MinVariance, UNSIGNED_2B_TYPE MinOccurrence)
//Removes dimensions, whose st.dev is less than MinVariance ( in total scale units)
//or whose #of non-zero signals is less than MinOccurrence
{
	UNSIGNED_4B_TYPE u, f, np = get_Ndatapoints(), nd = get_Ndimensions();
	apvector<REALNUM_TYPE> dv( np );
	QSAR qq;
	REALNUM_TYPE V, mxV, mnV;
	UNSIGNED_2B_TYPE occ0, occ1;
	set Dims(0, nd);

	for (u = 0; u < nd; u++)
	{
		for (f = 0; f < np; f++) dv[f] = patt(f, u);
		mxV = qq.maxV(dv); mnV = qq.minV(dv);

		for (occ0 = occ1 = f = 0; f < np; f++) 
		{//min.occurrence check
			if (dv[f] > mnV) occ1++;
			if (dv[f] < mxV) occ0++;
		}

		V = mxV - mnV;
		if (V > 0)	V = qq.stdev(dv) / V;
		if (V < MinVariance)	Dims.RemoveFromSet(u);
		if ((occ0 < MinOccurrence)||(occ1 < MinOccurrence))	Dims.RemoveFromSet(u);
	}//for u
	reduce_dimensions(Dims);
}

SIGNED_4B_TYPE dataset::ExpandDescriptors(dataset &X)
/*	Addes descriptor-matrix from X into the current dataset.
	NB: this may introduce redundant descriptors, 
	as well as descriptors that are scaled differently.
*/
{
	SIGNED_4B_TYPE r, c, dp = get_Ndatapoints(), dmn = get_Ndimensions(), Xdmn = X.get_Ndimensions();

	if ( dp != X.get_Ndatapoints() ) return INVALID;		//should be the same set of data points		
	if ( isScaled() != X.isScaled() ) 	return INVALID;		//scaling should match!	

	/* //as alternative, in case of scaling mismatch, we can undo the scaling:
		if ( (!isScaled()) && X.isScaled() ) X.scale_dimensions(3);
		if ( (!X.isScaled()) && isScaled() ) scale_dimensions(3);
	*/
	
	//fill up the expanded matrix
	patt.SetSize(dp, dmn + Xdmn);
	for (r = 0; r < dp; r++)	
	for (c = 0; c < Xdmn; c++)
		patt(r, dmn + c) = X.patt(r, c);
	
	//merge scaling info
	if (isScaled())
	{
		vars.A.resize(dmn + Xdmn);
		vars.B.resize(dmn + Xdmn);
		for (c = 0; c < Xdmn; c++)
		{
			vars.A[dmn + c] = X.vars.A[c];
			vars.B[dmn + c] = X.vars.B[c];
		}
	}

	//merge descriptor labels
	UNSIGNED_4B_TYPE a, m, n;
	vars.L.resize(dmn + Xdmn);
	for (c = 0; c < Xdmn; c++)
	{
		vars.L[dmn + c] = X.vars.L[c];		
		do 
		{
			a = vars.FindHashKey( vars.hashFunction(dmn + c) );
			n = vars.RetrieveRecord(a, m);
			while (n > 0)			
			{
				n--;
				if ( SIGNED_4B_TYPE(m + 1) > (dmn + c) ) continue;
				if (vars.L[m].length() == X.vars.L[c].length())
				if (vars.L[m] == X.vars.L[c])
				{
					vars.L[dmn + c] += "(1)";
					n = 1; //recheck the new label
					break;
				}
				vars.RetrieveRecord(a, m);
			}
		} while (n > 0);
		vars.AddRecord(dmn + c);
	}

	//invalidate distances
	dist.SetSize(0,0);
	distAv = distMax = distMin = distn0Min = DUMMY_VAL;
	return 0;
}

bool dataset::isScaled()
//returns true if scaling coefficients exist
{
	return (vars.A.length() > 0);
}
