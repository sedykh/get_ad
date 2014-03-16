// QSAR.cpp: implementation of various SAR and QSAR methods
//////////////////////////////////////////////////////////////////////

#include "qsar.h"
//HierCluster() in future should return a tree-structure explicitely, not apvector<apvector<set>>

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

QSAR::QSAR()
{

}

QSAR::~QSAR()
{

}

REALNUM_TYPE QSAR::Fvalue(REALNUM_TYPE TSS, REALNUM_TYPE RSS, UNSIGNED_4B_TYPE DFM, UNSIGNED_4B_TYPE DFE)
/*
DF - degrees of freedom; DFM - ... of a model; DFE - of errors
Example: 
1-degree of freedom model (e.g. linear regression)
y = f(x); e.g. y = ax + b;

RSS = Sum{ (yi - f(xi))^2 }; Residual Sum of Squares
ESS = Sum{ (<y> - f(xi))^2 }; Explained variance
TSS = Sum{ (yi - <y>)^2 }; Total Sum of Squares; = variance
TSS = ESS + RSS
R2 = ESS/TSS = 1 - RSS/TSS;

F = {(TSS - RSS) /DFM} / {RSS/DFE}, 
DFM = 1; DFE = #datapoints - 2;
*/
{
	REALNUM_TYPE FV = TSS - RSS;
	FV /= DFM * RSS;
	FV *= DFE;
	return FV;
}

REALNUM_TYPE QSAR::rankcorrel(apvector<REALNUM_TYPE> &V1, apvector<REALNUM_TYPE> &V2)
{//Spearman rank correlation, less sensitive to outliers at the tail of the distribution
	SIGNED_4B_TYPE i,j,i0, N1 = V1.length(), N2=V2.length();
	if (min(N1,N2) < 2) return INVALID;
	
	REALNUM_TYPE *F1	= GRAB_MEM_BLOCKS(REALNUM_TYPE, N1), *F2	= GRAB_MEM_BLOCKS(REALNUM_TYPE, N2), x;
	SIGNED_4B_TYPE *A1 = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N1), *A2 = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N2);
	for (i = 0; i < N1; i++){ A1[i] = i;	F1[i] = V1[i]; }
	for (i = 0; i < N2; i++){ A2[i] = i;	F2[i] = V2[i]; }
	QSortScore = F1;	qsort(A1, (size_t)N1, sizeof(UNSIGNED_4B_TYPE), QSortCompareGreater);
	QSortScore = F2;	qsort(A2, (size_t)N2, sizeof(UNSIGNED_4B_TYPE), QSortCompareGreater);
	
	//assign ranks
	apvector<REALNUM_TYPE> RV1(N1), RV2(N2);
	x = F1[A1[0]] - 1;	i0 = 1;
	for (i = 0; i < N1; i++) 
	{
		RV1[ A1[i] ] = i + 1;
		if (F1[ A1[i] ] > x)
		{
			if (i > i0)
			{
				x  = i + i0;	x /= 2; 
				for (j=i0-1; j<i; j++) RV1[ A1[j] ] = x;		
			}
			x = F1[ A1[i] ]; i0=i + 1;
		}
	}
	if (i > i0)
	{
		x  = i + i0;	x /= 2; 
		for (j=i0-1; j<i; j++) RV1[ A1[j] ] = x;		
	}


	x = F2[A2[0]] - 1;	i0 = 1;
	for (i = 0; i < N2; i++) 
	{
		RV2[ A2[i] ] = i + 1;
		if (F2[ A2[i] ] > x) 
		{
			if (i > i0)
			{
				x  = i + i0;	x /= 2; 
				for (j=i0-1; j<i; j++) RV2[ A2[j] ] = x;		
			}
			x = F2[ A2[i] ]; i0=i + 1;
		}
	}
	if (i > i0)
	{
		x  = i + i0;	x /= 2; 
		for (j=i0-1; j<i; j++) RV2[ A2[j] ] = x;		
	}

	DROP_MEM_BLOCKS(F1);
	DROP_MEM_BLOCKS(A1);
	DROP_MEM_BLOCKS(F2);
	DROP_MEM_BLOCKS(A2);
	QSortScore = NULL;

	return correl(RV1, RV2);
}

REALNUM_TYPE QSAR::correl(apvector<REALNUM_TYPE> &V1, apvector<REALNUM_TYPE> &V2)
{//Pearson correlation
	UNSIGNED_4B_TYPE N = min(V1.length(), V2.length()), Vx = ZERO;
	if (N < 2) return ZERO;
	if (V1.length() != V2.length()) return ZERO;

	REALNUM_TYPE vl = ZERO, vM1 = meanV(V1), vM2 = meanV(V2), vD = stdev(V1)*stdev(V2);
	if (vD > 0)
	{//correlate only if V1 and V2 have some variation!
		for (; Vx < N; Vx++)	vl += (V1[Vx] - vM1)*(V2[Vx] - vM2);
		vl /= (N - 1)*vD;
	}
	return (vl);
}

REALNUM_TYPE QSAR::q2etc(apvector<REALNUM_TYPE> &V_OBS, apvector<REALNUM_TYPE> &V_PRED, UNSIGNED_1B_TYPE ST_MODE)
/*	description:	calculates statistical factors, q2, q2', MAEq
	V_OBS[] - observed activities
	V_PRED[]- predicted activities
	ST_MODE - 0 = q2; 1 = q2'; 2 = MAEq; 3 = MAEq'
*/
{
	UNSIGNED_4B_TYPE N = min(V_OBS.length(), V_PRED.length());
	if (N < 2) return 0;

	bool ifAlt = ((ST_MODE & 1) == 1), ifErr = ((ST_MODE & 2) == 2);

	REALNUM_TYPE vl = RSS(V_OBS, V_PRED, ifErr);	

	if (ifErr)
	{
		if (ifAlt)
			vl	/= absdiffV(V_PRED);
		else
			vl	/= absdiffV(V_OBS);
	}
	else
	{
		vl /= N;
		if (ifAlt)
			vl /= varianceV(V_PRED, true);
		else
			vl /= varianceV(V_OBS, true);
	}
	
	return (1.0 - vl);
}

REALNUM_TYPE QSAR::q2F13(apvector<REALNUM_TYPE> &EXT_OBS, apvector<REALNUM_TYPE> &EXT_PRED, REALNUM_TYPE REF, UNSIGNED_1B_TYPE ST_MODE)
/*Calculates "F1" and "F3" versions of Q2 for external sets (Consonni et al, J Chem Inf Model 2009, 49, pp1669-1678)
ST_MODE:
		= 0 for F3 (normalization by training set variance),
		=1 for F1 (normalization by external set variance but using training set mean value)
REF should contain:
		RSS of training set (ST_MODE = 0), or 
		mean observed activity value of the training set (ST_MODE=1)

NB:		Both F1,F3 are training-set dependendent, i.e. same prediction results
		F1-version is non-additive, and too optimistic if EXT_OBS are very different from the training set mean
		F3-version is additive, can be reliably calculated even for 1 data point (if REF is not 0)
		F3-version can be viewed as RMSE normalized by training set variance
*/
{	
	REALNUM_TYPE vl = MSE(EXT_OBS, EXT_PRED);	
	if (ST_MODE == 0)
	{
		if (REF > 0) vl	/=	REF; else vl = 1.0;
	}
	else
		vl /= varianceVext(EXT_OBS, REF, true);
	
	return (1.0 - vl);
}

REALNUM_TYPE QSAR::meanV(apvector<REALNUM_TYPE> &V)
{	
	UNSIGNED_4B_TYPE N = V.length();
	if (N == ZERO) return ZERO;

	REALNUM_TYPE vl = V[ZERO];
	for (UNSIGNED_4B_TYPE Vx = 1; Vx < N; Vx++)	vl += V[Vx];
	return (vl / N);
}

REALNUM_TYPE QSAR::minV(apvector<REALNUM_TYPE> &V)
{
	UNSIGNED_4B_TYPE N = V.length();
	if (N == ZERO) return ZERO;

	REALNUM_TYPE vl = V[ZERO];
	for (UNSIGNED_4B_TYPE Vx = 1; Vx < N; Vx++)	
		if (vl > V[Vx]) vl = V[Vx];
	return (vl);
}
REALNUM_TYPE QSAR::maxV(apvector<REALNUM_TYPE> &V)
{
	UNSIGNED_4B_TYPE N = V.length();
	if (N == ZERO) return ZERO;

	REALNUM_TYPE vl = V[ZERO];
	for (UNSIGNED_4B_TYPE Vx = 1; Vx < N; Vx++)	
		if (vl < V[Vx]) vl = V[Vx];
	return (vl);
}

REALNUM_TYPE QSAR::sumV(apvector<REALNUM_TYPE> &V)
{//added 06.05.2013
	UNSIGNED_4B_TYPE N = V.length();
	if (N == ZERO) return ZERO;

	REALNUM_TYPE vl = V[ZERO];
	for (UNSIGNED_4B_TYPE Vx = 1; Vx < N; Vx++)	vl += V[Vx];
	return (vl);
}

REALNUM_TYPE QSAR::middleV(apvector<REALNUM_TYPE> &V)
{//calculates median, fixed 06.05.2013
	UNSIGNED_4B_TYPE i = 0, N = V.length();
	if (N == ZERO) return ZERO;
	if (N == 1) return V[0];
	if (N == 2) return ((V[0]+V[1])/2);

	REALNUM_TYPE *F	= GRAB_MEM_BLOCKS(REALNUM_TYPE, N), X;
	UNSIGNED_4B_TYPE *A = GRAB_MEM_BLOCKS(UNSIGNED_4B_TYPE, N);
	
	for (i = 0; i < N; i++){ A[i] = i;	F[i] = V[i]; };	
	QSortScore = F;	
	qsort(A, (size_t)N, sizeof(UNSIGNED_4B_TYPE), QSortCompareGreater);
	i = N >> 1; 
	if ((i << 1) == N) X = (V[A[i]] + V[A[i+1]])/2; else X = V[A[i]];
	DROP_MEM_BLOCKS(F);
	DROP_MEM_BLOCKS(A);
	QSortScore = NULL;
	return (X);
}

REALNUM_TYPE QSAR::absdiffV(apvector<REALNUM_TYPE> &V)
//calculates sum of absolute differences from average, similar to varianceV
{
	UNSIGNED_4B_TYPE N = V.length(), Vx = ZERO;
	if (N < 2) return ZERO;

	REALNUM_TYPE vl = ZERO, vMean = meanV(V);
	for (; Vx < N; Vx++)	vl += fabs(V[Vx] - vMean);
	return (vl);
}

REALNUM_TYPE QSAR::varianceVext(apvector<REALNUM_TYPE> &V, REALNUM_TYPE MV, bool Biased)
//calculates variance of V from supplied average value MV
//if V represents complete finite population then Biased mode can be used
{
	UNSIGNED_4B_TYPE N = V.length(), Vx = ZERO;	
	REALNUM_TYPE vl = 0;
	for (; Vx < N; Vx++)	vl += sqr(V[Vx] - MV);
	if (Biased) 
		vl /= max(1, N); 
	else 
		if (N > 1) vl /= (N - 1);
	return (vl);
}

REALNUM_TYPE QSAR::varianceV(apvector<REALNUM_TYPE> &V, bool Biased)
//calculates variance of V, 
//if V is complete finite population then Biased mode can be used
{
	return varianceVext(V, meanV(V), Biased);
}

REALNUM_TYPE QSAR::stdev(apvector<REALNUM_TYPE> &V)
{
	REALNUM_TYPE vse = varianceV(V);
	if (vse > 0) return sqrt(vse);
	return 0;
}

REALNUM_TYPE QSAR::RSS(apvector<REALNUM_TYPE> &V1, apvector<REALNUM_TYPE> &V2, bool AbsDiff)
{
	SIGNED_4B_TYPE N = min(V1.length(), V2.length()), Vx = 0;
	if (V1.length() != V2.length()) N = 0;
	REALNUM_TYPE vl = 0;
	for (; Vx < N; Vx++)	
	if (AbsDiff)
		vl += fabs(V1[Vx] - V2[Vx]);
	else
		vl += sqr(V1[Vx] - V2[Vx]);
	return (vl);
}

REALNUM_TYPE QSAR::MSE(apvector<REALNUM_TYPE> &V1, apvector<REALNUM_TYPE> &V2)
{//Mean Squared Error = RSS/n
	return RSS(V1, V2)/ max(V1.length(), 1);	
}

void QSAR::centralize(apvector<REALNUM_TYPE> &V, REALNUM_TYPE &COFF1, REALNUM_TYPE &COFF2)
{
	UNSIGNED_4B_TYPE N = V.length(), Vx = 0;
	COFF1 = COFF2 = ZERO;
	if (N < 2) return;
		
	COFF2 = meanV(V);
	COFF1 = maxV(V) - minV(V);
	
	if (COFF1 == ZERO) return;

	for (; Vx < N; Vx++)	
	{
		V[Vx] -= COFF2;
		V[Vx] /= COFF1;
	}
}


void QSAR::QR(matrix<REALNUM_TYPE> &A, matrix<REALNUM_TYPE> &Q, matrix<REALNUM_TYPE> &R)
//QR decomposition, where Q is orthogonal, R is upper triangular
{//Gram-Schmidt method
	UNSIGNED_4B_TYPE m = A.RowNO(), n = A.ColNO(), t, s;
	if (m < n) return;

	matrix<REALNUM_TYPE> e, u, a;
	REALNUM_TYPE fx;

	Q.Null(m, n);
	R.Null(n, n);
		
	for (t = ZERO; t < n; t++)
	{		
		a = u = A.GetCol(t);
		for (s = ZERO; s < t; s++)
		{
			e = Q.GetCol(s);
			u -= (e * ((~e) * a));
		}
		fx = u.Norm();
		if (fx > ZERO) 
			for (s = ZERO; s < m; s++)
				Q(s, t) = u(s, 0) / fx;		
	}
	//now Q is orthogonal
	R = (~Q) * A;
}


UNSIGNED_4B_TYPE QSAR::GetEigenVectors(matrix<REALNUM_TYPE> &A, matrix<REALNUM_TYPE> &EV, matrix<REALNUM_TYPE> &EL, REALNUM_TYPE nThreshold, UNSIGNED_4B_TYPE nMaxItr)
	//QR algorithm, Gershgorin circle is used for the end condition	
	//if A is symmetrical, EL will converge to diagonal matrix of eigenvalues, 
	//EV - to orthogonal matrix of right eigenvectors.
	//NOTE: convergence conditions are made proportional to the matrix dimension (n)
	//NOTE: this algorithm is slow for large matrices due to many expensive QR decompositions.
{
	EV.Null(); 
	EL.Null();
	if (A.ColNO() != A.RowNO()) return ZERO;
	
	matrix<REALNUM_TYPE> Ai, Qi, Ri;	
	UNSIGNED_4B_TYPE n = A.ColNO(), itr = ZERO;
	
	EL = A;
	EV.Unit(n);	
	do
	{
		Ai = EL;
		QR(EL, Qi, Ri);
		EL  = Ri * Qi;
		EV *= Qi;
	} while ((EL.modNorm(1) > nThreshold*n) && (++itr < nMaxItr*n));

	return (itr);
}
UNSIGNED_4B_TYPE QSAR::GetEigenVectorsEff(matrix<REALNUM_TYPE> &A, matrix<REALNUM_TYPE> &EV, matrix<REALNUM_TYPE> &EL, REALNUM_TYPE threshold, UNSIGNED_4B_TYPE MaxItr)
{	//Power iterations algorithm, may be used for large matrices.
	//similar to NIPALS algorithm of PCA
	//
	//if A is symmetrical, EL will be a diagonal matrix of eigenvalues, 
	//EV - to orthogonal matrix of right eigenvectors.
	//
	//NOTE: this algorithm gives eigenvectors that may differ in sign,
	//but are colinear with eigenvectors of GetEigenVectors();
	
	EV.Null(); EL.Null();
	if (A.ColNO() != A.RowNO()) return ZERO;	
	if (A.ColNO() < 3)	return GetEigenVectors(A, EV, EL, threshold);

	matrix<REALNUM_TYPE> Ai, Qi, Ri, rEigVec, lEigVec;
	REALNUM_TYPE Norm;
	UNSIGNED_4B_TYPE tot_itr = ZERO, itr, i, k, n = A.ColNO();
	
	Ai = A;
	EL.Unit(n);
	EV.Null(n, n);
	
	for (k = ZERO; k < n; k++)
	{
		rEigVec = Ai.GetCol(0);
		rEigVec /= rEigVec.Norm();
		lEigVec = Ai.GetRow(0);
		lEigVec /= lEigVec.Norm();

		Qi = Ai * Ai; //to do iterations by the power of 2
		itr = ZERO;
		do 
		{
			Ri = rEigVec;
			rEigVec  = Qi * rEigVec;
			rEigVec /= rEigVec.Norm();
			Ri		-= rEigVec;
			Norm	 = Ri.Norm();

			Ri = lEigVec;
			lEigVec *= Qi;
			lEigVec /= lEigVec.Norm();
			Ri		-= lEigVec;
			Norm	+= Ri.Norm();
			Norm	/= 2;
		} while ((Norm > threshold) && (++itr < MaxItr));

		tot_itr += itr;		
		Qi = (~rEigVec) * Ai * rEigVec;
		//copy eigenvalue and eigenvector
		EL(k,k) = Qi(0, 0);
		for (i = 0; i < n; i++)		
			EV(i, k) = rEigVec(i, 0);
		
		//reduce the matrix
		Ri = rEigVec * lEigVec * EL(k,k);	//same as Ri = Ai* rEigVec * lEigVec;	Ri = rEigVec * lEigVec * Ai;
		Ai -= Ri;		
	}//for k	
	
	return (tot_itr);
}

UNSIGNED_4B_TYPE QSAR::NIPALS(matrix<REALNUM_TYPE> &X, matrix<REALNUM_TYPE> &Score, matrix<REALNUM_TYPE> &Loading, REALNUM_TYPE threshold, UNSIGNED_4B_TYPE MaxItr)
//ref: Geladi, P., Kowalski, Analytica Chimica Acta, 185(1986)1-17
//NIPALS - Nonlinear iterative partial least squares, 
//calculates principal components for PCA and PCR
//for matrix X{n x m} => Score{n x n'}*Loading{n' x m}, n - number of rows, m - number of columns
//Score() columns and Loading'() rows are collinear with eigenvectors of XX' and X'X matrices
//
//to emulate NIPALS use: X{n x m}, n > m;
//GetEigenVectorsEff(X'X, Loading', ...)
//Score = X*Loading';

{	
	matrix<REALNUM_TYPE> Ei, Ri, t, p_;
	UNSIGNED_4B_TYPE tot_itr = ZERO, itr, i, k, m = X.ColNO(), n = X.RowNO();
	
	Ei = X;
	Score.Null(n, m);
	Loading.Null(m, m);
	
	for (k = ZERO; k < m; k++)
	{
		t	 = Ei.GetCol(0);
		itr = ZERO;
		do 
		{
			Ri = t;
			p_	 = (~t)*Ei;
			p_	/= p_.Norm();
			t	 = Ei*(~p_);
			Ri -= t;
		} while ((Ri.Norm() > threshold) && (++itr < MaxItr));
		tot_itr += itr;		
		
		//copy obtained components		
		for (i = 0; i < n; i++)	Score(i, k) = t(i, 0);
		for (i = 0; i < m; i++)	Loading(k, i) = p_(0, i);
		
		//reduce the matrix
		Ri = t * p_;	//same as Ri = Ai* rEigVec * lEigVec;	Ri = rEigVec * lEigVec * Ai;
		Ei -= Ri;
		if (Ei.Norm() < threshold) break;
	}//for k
	
	return (tot_itr);
}
UNSIGNED_4B_TYPE QSAR::PLSAlgorithm
(UNSIGNED_4B_TYPE NComponents, 
 matrix<REALNUM_TYPE> &E0,matrix<REALNUM_TYPE> &F0,//input matrices
 matrix<REALNUM_TYPE> &P,matrix<REALNUM_TYPE> &W,//results 
 matrix<REALNUM_TYPE> &Q, matrix<REALNUM_TYPE> &B,
 REALNUM_TYPE threshold, UNSIGNED_4B_TYPE MaxItr)
//ref: Geladi, P., Kowalski, Analytica Chimica Acta, 185(1986)1-17
//assuming that X and Y are mean centered and scaled
//X = TP' + residual; Y = UQ' + residual {similar to PCA, but with inner relationship U(i) = B(i)T(i)}
//Y = TBQ' + F; where ||F|| is to be  minimized.
//in order to obtain orthogonal X scores as in PCA it is necessary to introduce weights (w[])
//T = XW
{
	REALNUM_TYPE ff;
	UNSIGNED_4B_TYPE i, Itr, iret = ZERO;
	matrix<REALNUM_TYPE> u, wt, qt, t, pt, prevt, b, _m;	
	matrix<REALNUM_TYPE> X = E0, Y = F0;

	
	//B.Null(NComponents, 1);
	B.Null(NComponents, NComponents); //diagonal, March 27,2008

	//untransposed matrices:
	P.Null(X.ColNO(), NComponents);
	W.Null(X.ColNO(), NComponents);
	Q.Null(Y.ColNO(), NComponents);
	
	
	while (iret < NComponents)
	{
		t.Null();
		Itr = ZERO;
		
		u = Y.GetCol(0);						//(1) u-start = some y(j)
		if (u.Norm() < SMALL_NUMBER) break;
		
		do 
		{
			prevt = t;
			wt	 = (~u) * X;					//(2) w' = u'X/u'u
			wt	/= wt.Norm();					//(3) w' normalization
			
			_m	 = ~wt;							//(4) t = Xw/w'w
			t	 = X * _m;
			_m	 = wt * _m;
			t	/= _m.Norm();
			
			if (Y.ColNO() == 1)
			{//if only one column no need for iterations, q = 1, u is the original column
				qt.SetSize(1, 1); 
				qt(0, 0) = 1;
				break;
			}
		
			qt	 = (~t) * Y;					//(5) q' = t'Y/t't
			qt	/= qt.Norm();					//(6) q' normalization

			_m = ~qt;							//(7) u = Yq/q'q
			u	 = Y * _m;
			_m	 = qt * _m;
			u	/= _m.Norm();

			_m = t - prevt;						//(8) check convergence
		} while ((_m.Norm() > threshold) || (++Itr < MaxItr));

		_m = ~t;								//(9) p' = t'X/t't
		pt	 = _m * X;
		_m	*= t;
		pt	/= _m.Norm();

		ff = pt.Norm();							//normalizations
		pt	/= ff;								//(10) p
		t	/= ff;								//(11) t
		wt	/= ff;								//(12) w'
		

		_m	 = (~t) * t;						//(13) b = u't/t't
		b	 = (~u) * t;
		b	/= _m.Norm();

		//residuals
		X	-= t * pt;	
		Y	-= t * b * qt;

		for (i = ZERO; i < X.ColNO(); i++)
		{
			P(i, iret) = pt(0, i);
			W(i, iret) = wt(0, i);
		}
		
		for (i = ZERO; i < Y.ColNO(); i++)
			Q(i, iret) = qt(0, i);
		
		B(iret, iret) = b(0, 0); //diagonal, March 27,2008
		//B(iret, 0) = b(0, 0);

		iret++;
	}

	return iret;
}

void QSAR::getCrossMatrix(matrix<REALNUM_TYPE> &dataMtx, matrix<REALNUM_TYPE> &qMtx, REALNUM_TYPE metrPow, UNSIGNED_1B_TYPE mode)
//creates matrix of r-coefficients for the rows of dataMtx
//NB: to do the same for columns just submit transposed dataMtx
//mode: 0 - 3 distance-modes:
//		0 = Euclidean, 1 - Cosine, 2 - Correl, 3 - Tanimoto, 
//mode: > 3 is a normal coefficient mode:
//		4 - r; >4 - r^2
{
	REALNUM_TYPE r;
	UNSIGNED_4B_TYPE DN = dataMtx.RowNO(), CN = dataMtx.ColNO(), f, u, u1;
	apvector<REALNUM_TYPE> dv1(CN), dv2(CN);
	
	qMtx.SetSize(DN, DN);
	for (u = 0; u < DN - 1; u++)
	{
		if (mode < 4) 
			qMtx(u, u) = 0.0; //distance-mode
		else	
			qMtx(u, u) = 1.0; //correlation coefficient

		for (f = 0; f < CN; f++)	dv1[f] = dataMtx(u, f);

		for (u1 = u + 1; u1 < DN; u1++)
		{
			for (f = 0; f < CN; f++)	dv2[f] = dataMtx(u1, f);
			if (mode < 4)
				r = getMetricDistance(dv1, dv2, metrPow, mode);
			else
			{
				r = correl(dv1, dv2);
				if (mode > 4) r *= r;
			}
			qMtx(u, u1) = qMtx(u1, u) = r;
		} //for u1
	} //for u
}

void QSAR::HierCluster(matrix<REALNUM_TYPE> &Map, apvector<apvector_set_type> &Hier, UNSIGNED_4B_TYPE maxNLevels)
//precondition: Map(,) should be symmetric and square matrix
//description: Agglomerative clustering + single-linkage clustering
{
	Hier.resize(ZERO);
	if (!Map.IsSquare()) return;
	if (!Map.IsSymmetric()) return;

	SIGNED_4B_TYPE i, i1, j, j1, MN = Map.ColNO();
	REALNUM_TYPE minV, maxV;

	if (MN < 2) return;

	minV = maxV = Map(0,1); //since the matrix is symmetrical it is faster to find min/max directly
	for (i = 0; i < MN - 1; i++)
		for(j = i + 1; j < MN; j++)
		{
			if ( minV > Map(i, j) ) minV = Map(i, j);
			if ( maxV < Map(i, j) ) maxV = Map(i, j);
		}

	//calculate the size of a scanning step 
	REALNUM_TYPE mapMarge, mapStep = (maxV - minV);
	mapStep	/= 0.99 + maxNLevels; //1 + maxNLevels, but 0.99 -> to go a little over the maxV in the last iteration
	mapMarge = mapStep + minV;
	
	apvector<UNSIGNED_2B_TYPE> clustList1, clustList2; //auxiliary variables
	UNSIGNED_4B_TYPE iterHier = 0;
	Hier.resize(maxNLevels + 1);

	Hier[iterHier].resize(MN);
	for (i = 0; i < MN; i++)	Hier[iterHier][i].PutInSet(i);
	
	//major cycle
	while (iterHier <= maxNLevels)
	{//the clusters A and B are merged if exists at least one pair (a,b) < mapStep
		for (i = 0; i < Hier[iterHier].length() - 1; i++)
		{
			Hier[iterHier][i].GetList(clustList1);
			for(j = i + 1; j < Hier[iterHier].length(); j++)
			{
				Hier[iterHier][j].GetList(clustList2);
				for (i1 = 0; i1 < clustList1.length(); i1++)
					for (j1 = 0; j1 < clustList2.length(); j1++)
					if ( mapMarge  >= Map(clustList1[i1], clustList2[j1]) )//join i-th and j-th clusters
					{
						Hier[iterHier][j] |= Hier[iterHier][i];
						Hier[iterHier][i].Dump();
						goto skipps;
					}
			}// for j
skipps:
			j = j1 = i1 = 0;
		}//for i

		mapMarge += mapStep;
		
		//rearrange the clusters
		j = Hier[iterHier].length();
		i = 0;
		while (i < j)
			if (Hier[iterHier][i].IsEmpty())
				Hier[iterHier][i] = Hier[iterHier][--j];
			else
				i++;
		
		if (Hier[iterHier].length() == j) //no new clusters, skip
			continue;
		else
			if (j == 1)//don't save the last, global, cluster 
				break;
			else
			{
				Hier[iterHier].resize(j);
				iterHier++;
				Hier[iterHier] = Hier[iterHier - 1];
			}
	}//while iterHier

	Hier.resize(iterHier);

	if (iterHier < 2) 
		return;

	//now, sort elements in Hier, so that the layers are aligned
	REALNUM_TYPE *FF = NULL;
	UNSIGNED_4B_TYPE *AA = NULL;
	i = iterHier - 1; 
	while (i > 0) 
	{
		MN = Hier[i - 1].length();
		FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, MN );
		AA = GRAB_MEM_BLOCKS(UNSIGNED_4B_TYPE, MN );

		for (j = 0; j < MN; j++) 
		{//j scans the (i-1)th layer (the one to be sorted)
			AA[j] = j;
			Hier[i - 1][j].GetElement(i1);
			for (j1 = 0; j1 < Hier[i].length(); j1++)
				if ( Hier[i][j1].IsInSet(i1) )//check against the reference i-th layer
				{
					FF[j] = j1;
					break;
				}
		}// for j

		//----- now sort ------
		QSortScore = FF;
		qsort(AA, (size_t)MN, sizeof(UNSIGNED_4B_TYPE), QSortCompareGreater);

		DROP_MEM_BLOCKS(FF);
		apvector_set_type tLayer(MN);
		for (j = 0; j < MN; j++)	tLayer[j] = Hier[i - 1][j];
		for (j = 0; j < MN; j++)	Hier[i - 1][j] = tLayer[ AA[j] ];
		DROP_MEM_BLOCKS(AA);
		//---------------------

		i--;
	} // while i

	QSortScore = NULL;
}

void QSAR::get_conf_mtx(apvector<SIGNED_4B_TYPE> &V_OBS, apvector<SIGNED_4B_TYPE> &V_PRED, matrix<SIGNED_4B_TYPE> &MTX_CF)
//description:	calculates confusion matrix for VALS
//
//precondition:
//				MTX_CF should be square with proper dimensions set, NG X NG
//				V_OBS, V_PRED should be same size, containing integer values from 0 to NG - 1
//
//postcondition:
//				MTX_CF will be: predicted (as rows) vs. observed (as columns)
{
	SIGNED_4B_TYPE i, vN = V_OBS.length(), NG = MTX_CF.RowNO();
	if ( (vN != V_PRED.length()) || (!MTX_CF.IsSquare()) || (NG < 2) )
	{
		MTX_CF.SetSize(0, 0);
		return;
	}

	MTX_CF.Null(NG, NG);
	for (i = 0; i < vN; i++)
	if ( (V_OBS[i] < NG) && (V_PRED[i] < NG) )
		MTX_CF(V_PRED[i], V_OBS[i])++;
	else
	{
		MTX_CF.SetSize(0, 0);
		return;		
	}		
}

UNSIGNED_4B_TYPE QSAR::get_groupN(matrix<SIGNED_4B_TYPE> &MX, SIGNED_4B_TYPE i, bool ObsPred)
//description: 
//				calculates total size of the i-th group in confusion matrix MX
//				if i is out of range, then total size of all groups together is calculated
//
//				ObsPred: true - for observed, false - for predicted
{
	if ( !MX.IsSquare() ) return 0;

	SIGNED_4B_TYPE ix = i, c, NG = MX.RowNO(), RS = 0;
	bool ifTotal = ( (i >= NG) || (i < 0) );
	if (ifTotal) ix = 0;
	for (; ix < NG; ix++)
	{
		for (c = 0; c < NG; c++) if (ObsPred) RS += MX(c, ix); else	RS += MX(ix, c);
		if (!ifTotal) break;
	}
	return (RS);
}

REALNUM_TYPE QSAR::get_mcc(matrix<SIGNED_4B_TYPE> &MX)
//Matthew's Correlation Coefficient: from -1 to 1
//MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
//TP,TN,FP,FN -- TRUE/FALSE POSITIVE/NEGATIVE
//NB:	this is a special case of general Pearson's correlation coefficient for binary data!
//
//NB:	may be misleading if classes are VERY unbalanced by size:
//		i.e. if there are only few datapoints in the smallest class 
//		and 1 of them accidentally predicted right, then TP*TN is still very large!
{
	UNSIGNED_4B_TYPE NG = MX.RowNO();
	if ( (!MX.IsSquare()) || (NG < 2) )	return 0;
	REALNUM_TYPE mcc = MX.Det(), xx;
	if (mcc > 0)
	{
		xx = (MX(0,0)+MX(0,1))*(MX(0,0)+MX(1,0))*(MX(1,1)+MX(0,1))*(MX(1,1)+MX(1,0));
		mcc = sqrt(xx);
	}
	return mcc;
}

REALNUM_TYPE QSAR::get_ccr(matrix<SIGNED_4B_TYPE> &MX, apvector<REALNUM_TYPE> &WS, REALNUM_TYPE PN, UNSIGNED_1B_TYPE MD)
//description:	calculates various performance-parameters for confusion matrix MX
//				MD & 1 -> choice of weighing scheme (0 - CCR, 1 - Accuracy)
//				MD & 2 -> choice of equation (0 - CCR/accuracy, 1 - error based)
//				MD & 4 -> for error based only: to use aver.err instead of max.err
//				MD & 8 -> for error based only: to use class-kind of term: "1" instead of "abs(i-j)"
//preconditions: 
//				some of WS[] > 0; dimensions should match
//
//result:		returns performance parameter ( <= 1), expected (random) performance for different versions:
//				CCR/Accuracy mode normally gives ~1/NG, Max.Error mode ~0.5, Av.Error mode ~ 0
//				ErrClass by itself: ~1/NG, with Av.Error gives ~ 0
{
	SIGNED_4B_TYPE i, j, wN = WS.length(), NG = MX.RowNO(), nn;	
	if ( (!MX.IsSquare()) || (NG < 2) || (wN != NG) )	return 0;

	bool ifPenlt = (PN != 0), ifAvErr = ((MD & 4) == 4), 
		ifErr = ((MD & 2) == 2), ifWgtAlt = ((MD & 1) == 1), ifErrClass = ((MD & 8) == 8);
	REALNUM_TYPE xn = 0, xd = 0, nnn = 0;

	apvector<UNSIGNED_4B_TYPE> gin(NG);
	for (i = 0; i < NG; i++) gin[i] = get_groupN(MX, i);

	if (ifErr) 
	{//calculate error-based
		for (i = 0; i < NG; i++)
		{
			for (nnn = nn = j = 0; j < NG; j++) 
			{
				SIGNED_4B_TYPE df = abs(i - j);
				if (ifErrClass) df = min(df, 1); //will be 0 when i=j; NB: w/o ifAvErr ifErrClass => the same as CCR/Accuracy

				nn	+= df * MX(j, i);
				if (ifAvErr) nnn += df;
			}
			
			//----June - October 2010 fix
			if (ifAvErr) nnn /= NG;	else	nnn = (ifErrClass) ? 1 : max(NG - i - 1, i);
			nnn *= gin[i];
			//----end of mdf

			if (ifWgtAlt)
			{
				xn += WS[i] * nn;				
				xd += WS[i] * nnn;
			}
			else
				if (gin[i])
				{
					xn += (WS[i] * nn) / gin[i];
					xd += (WS[i] * nnn) / gin[i];
				}
		}//for i

		xn = 1.0 - xn / xd;		
	}
	else
	{//calculate normal CCR-based accuracy
		for (i = 0; i < NG; i++)
		if (gin[i])
		{
			if (ifWgtAlt)
			{
				xn += WS[i] * MX(i, i);
				xd += WS[i] * gin[i];
			}
			else
				xn += WS[i] * MX(i, i) / gin[i];
		}
		if (ifWgtAlt) xn /= xd;
	}

	//now xn has the score-parameter
	xd = 0;
	if (ifPenlt)
	{//calculate penalty
		apvector<REALNUM_TYPE> ccr( NG );
		for (i = 0; i < NG; i++) ccr[i] = get_ccri(MX, i);
		for (i = 0; i < NG; i++)
		if (gin[i])
		{
			for (j = i + 1; j < NG; j++)
			if (gin[j] && (fabs(ccr[i] - ccr[j]) > xd) )
				xd = fabs(ccr[i] - ccr[j]);
		}
		xd *= PN;
	}

	return (xn - xd);
}



REALNUM_TYPE QSAR::get_ccri(matrix<SIGNED_4B_TYPE> &MX, SIGNED_4B_TYPE gi, bool MD)
//description:	calculates CCR/Accuracy for class i or the overall one if i is out of range;
//				MD -> choice of scheme (0 - CCR, 1 - Accuracy), 
//				will make difference only for global CCR/Accuracy
//
//result:		returns performance parameter ( < 1)
{	
	SIGNED_4B_TYPE NG = MX.RowNO();
	UNSIGNED_4B_TYPE NIG;
	if ( (!MX.IsSquare()) || (NG < 2) )	return 0;

	if ( (gi >= 0) && (gi < NG) )
	{
		NIG = get_groupN(MX, gi);
		if (NIG) return (REALNUM_TYPE(MX(gi, gi)) / NIG); 
		return 0;
	}

	REALNUM_TYPE xn = 0, xd = 0;
	for (SIGNED_4B_TYPE i = 0; i < NG; i++)
	{
		NIG = get_groupN(MX, i);
		if (NIG == 0) continue;	
		if (MD)
		{//accuracy
			xn += MX(i, i);
			xd += NIG;
		}
		else
			xn += REALNUM_TYPE(MX(i, i)) / NIG; //CCR
	}
	if (!MD) xd = NG;
	xn /= max(xd, 1);
	return xn;
}

REALNUM_TYPE QSAR::get_max_density_cutoff(apvector<REALNUM_TYPE> &dp, REALNUM_TYPE startx)
//NB:	dp() will be treated as non-negative values!
//program finds the dx-value that maximizes f()= #points on [minx..minx+dx] / dx
{
	SIGNED_4B_TYPE i, j, dpN = dp.length();
	if (dpN < 2) return INVALID;
	
	REALNUM_TYPE *FF	= GRAB_MEM_BLOCKS(REALNUM_TYPE, dpN);
	SIGNED_4B_TYPE *AA = GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, dpN);
	for (i = 0; i < dpN; i++)
	{
		AA[i] = i;
		FF[i] = fabs(dp[i]);
	}

	QSortScore = FF;
	qsort(AA, (size_t)dpN, sizeof(UNSIGNED_4B_TYPE), QSortCompareGreater);
	
	//selecting scanning range! ------------------

	//find min non-0 diff b/w datapoints
	REALNUM_TYPE dx, rtv, rtdv;
	rtdv = dx = FF[ AA[dpN-1] ];
	for (i = 1; i < dpN; i++)
	{
		rtv = FF[ AA[i] ] - FF[ AA[i-1] ];
		if  ((rtv > 0) && (rtv < rtdv)) rtdv = rtv;
	}
	rtdv*= 2; //double it
	
	rtv	 = FF[ AA[dpN-1] ] - FF[ AA[0] ];
	rtv	/= dpN << 1;
	rtdv = max(rtv, rtdv); //compare with half aver
	//---------------------------------------------

	j = i = 0;

	rtv = FF[ AA[0] ];
	if ((startx > 0) && (startx < 1)) //start scanning from specified point!
		rtv += (FF[ AA[dpN-1] ] - FF[ AA[0] ]) * startx;

	rtv +=  rtdv;	
	while ( rtv < FF[ AA[dpN-1] ] )
	{
		while ( (i < dpN) && (FF[ AA[i] ] < rtv) ) i++;
		if (j * rtv < i * dx) 
		{
			dx = rtv;
			j = i;
		}
		rtv += rtdv;
	}

	DROP_MEM_BLOCKS(FF);
	DROP_MEM_BLOCKS(AA);
	QSortScore = NULL;

	return dx;
}

REALNUM_TYPE QSAR::SS(apvector<REALNUM_TYPE> &V)
//Sum of Squares
{
	REALNUM_TYPE rtSS = 0;
	for (SIGNED_4B_TYPE i=0; i < V.length(); i++)	rtSS += V[i] * V[i];
	return rtSS;
}

REALNUM_TYPE QSAR::trendline(apvector<REALNUM_TYPE> &Y, apvector<REALNUM_TYPE> &X, REALNUM_TYPE &B)
/*	calculates best fitting trendline for Y = kX + b
	return value is slope, k, the constant term is returned in B

	k = Sum{ (Xj - <X>)(Yj - <Y>) }/Sum{ (Xj - <X>)^2 }
	b = <Y> - a<X>
*/
{
	SIGNED_4B_TYPE N  = X.length();
	REALNUM_TYPE rtK = 0;
	B = 0;
	if (N == Y.length())
	{
		REALNUM_TYPE avX = meanV(X), avY = meanV(Y), rtL = 0;
		for (SIGNED_4B_TYPE i = 0; i < N; i++)
		{
			rtK += (X[i] - avX)*(Y[i] - avY);
			rtL += sqr(X[i] - avX);
		}
		if (rtL > 0)	rtK /= rtL;
		B = avY - rtK*avX;
	}
	return rtK;
}

REALNUM_TYPE QSAR::trendline0(apvector<REALNUM_TYPE> &Y, apvector<REALNUM_TYPE> &X)
/*	calculates best fitting trendline for Y = kX
	return value is slope, k

	k = Sum(XjYj) / Sum(Xj^2)
*/
{
	REALNUM_TYPE rtK = 0;
	if (X.length() == Y.length())
	{
		for (SIGNED_4B_TYPE i = 0; i < X.length(); i++)	rtK += X[i]*Y[i];
		if (rtK > 0)	rtK /= SS(X);
	}
	return rtK;
}

REALNUM_TYPE QSAR::sqrR0(apvector<REALNUM_TYPE> &Y, apvector<REALNUM_TYPE> &X)
/*
	returns R2 for correlation through the origin; Y = kX;
*/
{
	if (X.length() != Y.length()) return 0;
	SIGNED_4B_TYPE N = min(X.length(), Y.length());
	REALNUM_TYPE k = trendline0(Y, X);
	apvector<REALNUM_TYPE> kX(N);
	for (SIGNED_4B_TYPE i = 0; i < N; i++)	kX[i] = k*X[i];
	return q2etc(Y, kX);
}

void QSAR::remove_unpredicted(apvector<REALNUM_TYPE> &e, apvector<REALNUM_TYPE> &p, REALNUM_TYPE MaxV)
//service function to remove unpredicted datapoints
//e - experimental, p - predicted values. 
//MaxV is used to find them.
{
	if (p.length() != e.length()) return;
	SIGNED_4B_TYPE i=0, dpn = p.length();
	while (i < dpn)
	if (p[i] > MaxV)
	{
		p[i] = p[--dpn];
		e[i] = e[dpn];
	}
	else i++;	
	e.resize(dpn);
	p.resize(dpn);
}

REALNUM_TYPE QSAR::erf(REALNUM_TYPE x)
/*	approximation of the error function based on
	Milton Abramowitz and Irene Stegun (1964) 
	Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables.
	
	erf() ~ 1 - (a1*t + a2*t^2 + ... )*exp(-x^2); 	t = 1/(1 + px); 
	1) for 3 terms: max.err 2.5*10^-5 and p=0.47047, a1=0.3480242, a2=-0.0958798, a3=0.7478556
	2) for 5 terms: max.err 1.5*10^-7 and p=0.3275911, a1=0.254829592, a2=-0.284496736, a3=1.421413741, a4=-1.453152027, a5=1.061405429

	another possibility:
	erf() ~ sign(x)*sqrt{ 1-exp(-x^2*(4/pi + ax^2)/(1+ax^2)) }, a = 8(pi - 3)/(3pi*(4 - pi)) ~ 0.140012
	which is very accurate in the neighborhood of 0 and infinity, while the overall max.err is 3.5*10^-4, 
														max.err can be 1.2*10^-4 if alternative a = 0.147 is used
*/
{
	if (x == 0) return 0;
	/*
	const REALNUM_TYPE a	= 0.14001228868666;	//8(pi - 3)/(3pi*(4 - pi))
	const REALNUM_TYPE Fpi	= 1.27323954473516;	//4/pi	
	REALNUM_TYPE x2 = x * x;
	REALNUM_TYPE ret	 = a*x2;
	x2	*= Fpi + ret;
	x2	/= 1 + ret;
	ret = sqrt(1 - exp(-x2)); //for large x2, precision may be lost at 1 - exp() step
	if (x < 0) return (-ret);
	return (ret);
	*/
	
	const REALNUM_TYPE p = 0.47047;
	const REALNUM_TYPE a1= 0.3480242;
	const REALNUM_TYPE a2=-0.0958798;
	const REALNUM_TYPE a3= 0.7478556;
	REALNUM_TYPE t	 = 1;
	if (x > 0) 
		t	/= p*x + 1;
	else
		t	/= 1 - p*x;

	REALNUM_TYPE t2	 = t * t;
	
	REALNUM_TYPE ret = exp(-x*x);
	ret	*= a1*t + a2*t2 + a3*t2*t;
	if (x < 0) return ret;
	return (1.0 - ret);
}

REALNUM_TYPE QSAR::norm_cdf(REALNUM_TYPE xval, REALNUM_TYPE mu, REALNUM_TYPE sigma)
// phi, gaussian cumulative distribution function based on the above erf() approximation
//	NB, direct calculation of CDF as 0.5[1 + erf(x2/sqrt2)] is not accurate for large -x 
//		due to floating number precision issues.
{
	if (sigma <= 0) return INVALID;

	REALNUM_TYPE x = xval - mu;
	if (x == 0) return 0.5;

	x	/= sigma;
	const REALNUM_TYPE p = 0.3326725;	//0.47047 / sqrt(2);
	const REALNUM_TYPE a1= 0.1740121;
	const REALNUM_TYPE a2=-0.0479399;
	const REALNUM_TYPE a3= 0.3739278;
	
	REALNUM_TYPE t	 = 1;
	if (x > 0)	t /= p*x + 1;	else	t /= 1 - p*x;

	REALNUM_TYPE t2	 = t * t, ret = exp(-x*x/2);
	ret	*= a1*t + a2*t2 + a3*t2*t;
	
	if (x < 0) return ret;
	return (1.0 - ret);
}
