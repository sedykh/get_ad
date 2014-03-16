// QSAR.h: interface for the QSAR class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(QSAR_CLASS)
#define QSAR_CLASS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "core.h"

//constants for PLS and PCA
#define DEFAULT_NCOMPONENTS			10		//PLS
#define DEFAULT_MIN_EIGENVALUE		0.1		//PCA

class QSAR  
{
public:
	QSAR();
	virtual ~QSAR();
	
	REALNUM_TYPE Fvalue(REALNUM_TYPE, REALNUM_TYPE, UNSIGNED_4B_TYPE, UNSIGNED_4B_TYPE);
	REALNUM_TYPE rankcorrel(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);	
	REALNUM_TYPE correl(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);
	REALNUM_TYPE q2etc(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &, UNSIGNED_1B_TYPE = 0);
	REALNUM_TYPE q2F13(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &, REALNUM_TYPE REF = 0, UNSIGNED_1B_TYPE = 0);

	REALNUM_TYPE trendline(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &, REALNUM_TYPE &);
	REALNUM_TYPE trendline0(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);

	REALNUM_TYPE sqrR0(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);

	REALNUM_TYPE middleV(apvector<REALNUM_TYPE> &);	//fixed 06.05.2013
	REALNUM_TYPE sumV(apvector<REALNUM_TYPE> &);	//added 06.05.2013
	REALNUM_TYPE maxV(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE minV(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE meanV(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE absdiffV(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE varianceV(apvector<REALNUM_TYPE> &, bool = false);
	REALNUM_TYPE varianceVext(apvector<REALNUM_TYPE> &, REALNUM_TYPE, bool = false);
	REALNUM_TYPE stdev(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE SS(apvector<REALNUM_TYPE> &);
	REALNUM_TYPE RSS(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &,  bool = false);
	REALNUM_TYPE MSE(apvector<REALNUM_TYPE> &, apvector<REALNUM_TYPE> &);
	void centralize(apvector<REALNUM_TYPE> &, REALNUM_TYPE&, REALNUM_TYPE&);
	
	REALNUM_TYPE erf(REALNUM_TYPE);	//fast and lossy approximation of the error function
	REALNUM_TYPE norm_cdf(REALNUM_TYPE, REALNUM_TYPE = 0, REALNUM_TYPE = 1);	//phi(x, mu, sigma), gaussian cumulative distribution function based on the above erf()

	//matrix based algorithms
	void QR(matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &);
	UNSIGNED_4B_TYPE GetEigenVectors(matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, REALNUM_TYPE = 0.0001, UNSIGNED_4B_TYPE = 100);
	UNSIGNED_4B_TYPE GetEigenVectorsEff(matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, REALNUM_TYPE = 0.0001, UNSIGNED_4B_TYPE = 100);
	UNSIGNED_4B_TYPE NIPALS(matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, REALNUM_TYPE = 0.0001, UNSIGNED_4B_TYPE = 100);	
	UNSIGNED_4B_TYPE PLSAlgorithm(UNSIGNED_4B_TYPE , matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, 
		matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, REALNUM_TYPE = 0.0001, UNSIGNED_4B_TYPE = 100);

	void getCrossMatrix(matrix<REALNUM_TYPE> &, matrix<REALNUM_TYPE> &, REALNUM_TYPE = 2.0, UNSIGNED_1B_TYPE = 0);

	//SAR
	void HierCluster(matrix<REALNUM_TYPE> &, apvector<apvector_set_type>&, UNSIGNED_4B_TYPE);
	void get_conf_mtx(apvector<SIGNED_4B_TYPE> &, apvector<SIGNED_4B_TYPE> &, matrix<SIGNED_4B_TYPE> &);
	UNSIGNED_4B_TYPE get_groupN(matrix<SIGNED_4B_TYPE> &, SIGNED_4B_TYPE = -1, bool = true);
	REALNUM_TYPE get_ccr(matrix<SIGNED_4B_TYPE> &, apvector<REALNUM_TYPE> &, REALNUM_TYPE = 0.0, UNSIGNED_1B_TYPE = 0);
	REALNUM_TYPE get_ccri(matrix<SIGNED_4B_TYPE> &, SIGNED_4B_TYPE = -1, bool = false);

	REALNUM_TYPE get_mcc(matrix<SIGNED_4B_TYPE> &);
	REALNUM_TYPE get_max_density_cutoff(apvector<REALNUM_TYPE> &, REALNUM_TYPE = 0);

	//service
	void remove_unpredicted(apvector<REALNUM_TYPE> &e, apvector<REALNUM_TYPE> &p, REALNUM_TYPE MaxV);
};

#endif // !defined(QSAR_CLASS)
