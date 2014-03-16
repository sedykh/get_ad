// core.cpp : Defines the entry point for the DLL application.
//
#include "core.h"

//-------   memory leaks catcher for the current source-file  --------
#ifdef ADV_LEAK_CATCHER
#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
#endif
//--------------------------------------------------------------------


//shared variables
apvector<Elements> ATABLE;				//table of atomic elements

apvector <Ring>	   rings, backup_rings;	//structure to store rings, used by molecule-class
apvector<SIGNED_4B_TYPE>  ringSys;		//structure to store ring-systems (it points to rings[])

apvector <Ring>	 delocs, contours;		//stucture to store pi-contours
apvector <atomic_in4cycl> CMAtoms;		//to store cyclic information

//structures for finding and storing hits\subgraphs, etc.
apvector<Cpair> A2Corrs;
apvector<Chain> Corrs;

set OrganicElements, AlkaliMetals, Metals, Halogens, Hbonders,
AromaticityElements, Group5Elements, Group6Elements, 
Loc, Deloc, PiBonds, CovBonds,
HphobResidues, HphilResidues, AmbiResidues, AromResidues;

bool SupressMessages;
STRING_TYPE StartDirName;	//current directory
STRING_TYPE ModulePath;		//a full path to the place where executable and configuration files are
STRING_TYPE LOG_FILENAME;	//log filename

REALNUM_TYPE getMetricDistance(apvector<REALNUM_TYPE> &V1, apvector<REALNUM_TYPE> &V2, REALNUM_TYPE metricPower, UNSIGNED_1B_TYPE metricFunc)
/*
description: calculates various metric distances depending on metricFunc value:
metricFunc = 0 - Euclidean/Minkowski; 1 - Cosine, 2 - Correlation, 3 - Tanimoto

in case of using similarity-based coeffients, they are transformed into distances by using
'1 - () ^ metricPower' equation */
{
	REALNUM_TYPE X, diff, diff1;
	UNSIGNED_4B_TYPE j, D = min( V1.length(), V2.length() );

	if (metricFunc == 0)
	{//Euclidean, general Minkowski
		for (X = j = 0; j < D; j++)	
		{
			diff = fabs(V1[j] - V2[j]);
			if (diff == 0)	continue;
			X += pow(diff, metricPower);
		}
		if (X == 0)	return X;
		return pow(X, 1.0/metricPower);
	}

	//non-Euclidean
	REALNUM_TYPE rtMean = 0, rtMean1 = 0;
	if (metricFunc == 2)
	{//Correlation coff-based. Calculate mean-centered values
		for (j = 0; j < D; j++)	
		{
			rtMean += V1[j]; 
			rtMean1 += V2[j];
		}
		rtMean	/= D;
		rtMean1	/= D;
	}

	for (diff = diff1 = X = j = 0; j < D; j++)
	{
		X += (V1[j] - rtMean)*(V2[j] - rtMean1);
		diff += sqr(V1[j] - rtMean);
		diff1 += sqr(V2[j] - rtMean1);
	}

	if ( (metricFunc < 3) )
	{//cosine similarity coeff, which is r-correlation if vectors are mean-centered
		diff *= diff1;
		if (diff > 0)	X	/= sqrt(diff); 	else X = 0; //make totally dissimilar in exceptional cases
	}

	if (metricFunc == 3)
	{//Tanimoto
		X = fabs(X);
		diff += diff1 - X;
		if (diff != 0)	X	/=	diff; else X = 0; //make totally dissimilar in exceptional cases
	}

	if (X == 0) return 2.0;
	return (1.0 - pow(X, metricPower));	
}

REALNUM_TYPE RRound(REALNUM_TYPE x)
{
	if (fabs(x- floor(x)) > 0.5)
		return ceil(x);
	return floor(x);
}

SIGNED_4B_TYPE Round(REALNUM_TYPE x)
{	
	return ( (SIGNED_4B_TYPE)RRound(x) );
}

#define scrll_nxt_s		for (ii = ii1; F[ii] == ' '; F[ii++] = '_');	ii1 = F.find(' '); xs = F.substr(ii, ii1 - ii)

void Load1El(STRING_TYPE F, UNSIGNED_2B_TYPE szATABLE)
{
	STRING_TYPE xs;
	UNSIGNED_2B_TYPE ii, ii1 = ZERO;
	UNSIGNED_2B_TYPE i, i1, j, k;

	scrll_nxt_s;
	ATABLE[szATABLE].SHORTNAME = xs; 	scrll_nxt_s;
	ATABLE[szATABLE].FULLNAME = xs;		scrll_nxt_s;	
	ATABLE[szATABLE].ATOMIC_NUMBER = atoi(xs.c_str());		scrll_nxt_s;
	ATABLE[szATABLE].ATOMIC_WEIGHT = atof(xs.c_str());		scrll_nxt_s;
	ATABLE[szATABLE].ATOMIC_RADIUS = atof(xs.c_str());		scrll_nxt_s;
	ATABLE[szATABLE].ATOMIC_ENGTV = atof(xs.c_str());		scrll_nxt_s;
	ATABLE[szATABLE].ATOMIC_IP1 = atof(xs.c_str());			scrll_nxt_s;
	ATABLE[szATABLE].ATOMIC_IP2 = atof(xs.c_str());			scrll_nxt_s;
	
	//read valencies	
	j = atoi(xs.c_str());
	ATABLE[szATABLE].STATES.resize(j);

	for (i = 0; i< j; i++)
	{
			scrll_nxt_s;
			ATABLE[szATABLE].STATES[i].VAL = (UNSIGNED_1B_TYPE)atoi(xs.c_str());	
			scrll_nxt_s;			
			i1 = atoi(xs.c_str());	;			//read configurations
			ATABLE[szATABLE].STATES[i].CFG.resize(i1);

			for (k = 0; k< i1; k++)
			{
				scrll_nxt_s;
				ATABLE[szATABLE].STATES[i].CFG[k].CONFG = (UNSIGNED_1B_TYPE)atoi(xs.c_str()); //coordination number
				scrll_nxt_s;
				ATABLE[szATABLE].STATES[i].CFG[k].COV_R = atof(xs.c_str());
				scrll_nxt_s;
				ATABLE[szATABLE].STATES[i].CFG[k].E_NGTV = atof(xs.c_str());				
			};	
	};
}

bool LoadElements()
//Description:  loads information about 
//				elements: names, weights, etc
//				also initializes OrganicElements, AlkaliMetals
//				Metals, Halogens, AromaticityElements, Group6Elements, etc.
//				set-structures
//				and initiates some common-use variables
{
	if (ATABLE.length() == 104)//already loaded
		return true;

	UNSIGNED_2B_TYPE szATABLE = 1; //skip 0 element - for the lone pair
	ATABLE.resize(103 + szATABLE);

	Load1El("H HYDROGEN 1 1.00794 0.371 2.20 13.598 0 1 1 1 1 0.30 2.00", szATABLE++);
	Load1El("He HELIUM 2 4.002602 0.50 0.00 24.587 54.4 0", szATABLE++);
	Load1El("Li LITHIUM 3 6.941 1.45 0.98 5.392 75.6 1 1 1 1 1.45 1.40", szATABLE++);
	Load1El("Be BERYLLIUM 4 9.012182 1.05 1.57 9.322 18.2 1 2 1 2 1.05 1.57", szATABLE++);
	Load1El("B BORON 5 10.811 0.82 2.04 8.298 25.2 1 3 1 3 0.81 2.08", szATABLE++);
	Load1El("C CARBON 6 12.011 0.771 2.55 11.26 24.4 1 4 3 2 0.60 3.15 3 0.67 2.1 4 0.77 2.10", szATABLE++);
	Load1El("N NITROGEN 7 14.00674 0.75 3.04 14.534 29.6 2 3 3 1 0.70 2.59 2 0.7 2.59 3 0.7 2.59 5 3 2 0.7 4.37 3 0.7 4.37 4 0.7 4.37", szATABLE++);
	Load1El("O OXYGEN 8 15.9994 0.74 3.44 13.618 35.1 2 2 2 1 0.62 4.34 2 0.66 3.2 4 1 2 0.62 4.34", szATABLE++);
	Load1El("F FLUORENE 9 18.9984035 0.706 3.98 17.422 35 1 1 1 1 0.64 4.00", szATABLE++);
	Load1El("Ne NEON 10 20.1797 0.65 0.00 21.564 41 0", szATABLE++);
	Load1El("Na NATRIUM 11 22.989767 1.80 0.93 5.139 47.3 1 1 1 1 1.80 0.93", szATABLE++);
	Load1El("Mg MAGNESIUM 12 24.3050 1.50 1.31 7.646 15 1 2 1 2 1.50 1.27", szATABLE++);
	Load1El("Al ALUMINIUM 13 26.981539 1.25 1.61 5.986 18.8 1 3 1 3 1.25 1.61", szATABLE++);
	Load1El("Si SILICON 14 28.0855 1.173 1.90 8.151 16.3 1 4 1 4 1.11 1.99", szATABLE++);
	Load1El("P PHOSPHORUS 15 30.973763 1.10 2.19 10.486 19.7 2 3 3 1 1.10 2.20 2 1.1 2.2 3 1.1 2.2 5 2 4 1.1 2.26 5 1.1 3.13", szATABLE++);
	Load1El("S SULFUR 16 32.066 1.04 2.58 10.36 23.3 3 2 2 1 0.94 4.48 2 1.04 2.74 4 2 3 0.94 2.79 4 1.04 2.58 6 3 3 1.04 2.58 4 1.04 2.58 6 1.04 2.58", szATABLE++);
	
	//Load1El("Cl CHLORINE 17 35.4527 0.994 3.16 12.967 23.8 1 1 1 1 0.99 3.28", szATABLE++);
	Load1El("Cl CHLORINE 17 35.4527 0.994 3.16 12.967 23.8 4 1 1 1 0.99 3.28 3 1 4 0.99 3.28 5 1 4 0.99 3.28 7 1 4 0.99 3.28", szATABLE++);

	Load1El("Ar ARGON 18 39.948 0.95 0.00 15.759 27.6 0", szATABLE++);
	Load1El("K POTASSIUM 19 39.0983 2.20 0.82 4.341 31.6 1 1 1 1 2.20 0.82", szATABLE++);
	Load1El("Ca CALCIUM 20 40.078 1.80 1.00 6.113 11.9 1 2 1 2 1.80 1.00", szATABLE++);
	Load1El("Sc SCANDIUM 21 44.955910 1.60 1.36 6.54 12.8 1 3 1 3 1.60 1.36", szATABLE++);
	Load1El("Ti TITANIUM 22 47.867 1.40 1.54 6.82 13.6 1 4 1 4 1.40 1.54", szATABLE++);
	Load1El("V VANADIUM 23 50.9415 1.35 1.63 6.74 14.6 4 2 1 2 1.35 1.63 3 1 3 1.35 1.63 4 1 3 1.35 1.63 5 1 4 1.35 1.63", szATABLE++);
	Load1El("Cr CHROMIUM 24 51.9961 1.40 1.66 6.766 16.5 3 2 1 2 1.40 1.66 3 1 3 1.4 1.66 6 1 6 1.4 1.66", szATABLE++);
	Load1El("Mn MANGANESE 25 54.93805 1.40 1.55 7.435 15.6 3 2 1 2 1.40 1.55 6 1 6 1.4 1.55 7 1 6 1.4 1.55", szATABLE++);
	
	//Load1El("Fe IRON 26 55.845 1.40 1.83 7.87 16.2 2 2 1 2 1.40 1.83 3 1 3 1.4 1.83", szATABLE++);
	Load1El("Fe IRON 26 55.845 1.40 1.83 7.87 16.2 2 2 3 2 1.40 1.83 4 1.40 1.83 6 1.40 1.83 3 3 3 1.4 1.83 4 1.4 1.83 6 1.4 1.83", szATABLE++);

	Load1El("Co COBALT 27 58.9332 1.35 1.88 7.86 17.1 2 2 1 2 1.35 1.88 3 1 3 1.35 1.88", szATABLE++);
	Load1El("Ni NICKEL 28 58.6934 1.35 1.91 7.635 18.2 1 2 1 2 1.35 1.91", szATABLE++);
	Load1El("Cu COPPER 29 63.546 1.35 1.90 7.726 20.3 1 2 1 2 1.35 1.90", szATABLE++);
	Load1El("Zn ZINC 30 65.39 1.35 1.65 9.394 18 1 2 1 2 1.35 1.65", szATABLE++);
	Load1El("Ga GALLIUM 31 69.723 1.30 1.81 5.999 20.5 1 3 1 3 1.30 1.77", szATABLE++);
	Load1El("Ge GERMANIUM 32 72.61 1.22 2.01 7.899 15.9 1 4 1 4 1.22 2.06", szATABLE++);
	Load1El("As ARSENIC 33 74.92159 1.21 2.18 9.81 18.6 2 3 1 3 1.21 2.38 5 1 4 1.21 2.69", szATABLE++);
	
	//Load1El("Se SELENIUM 34 78.96 1.17 2.55 9.752 21.2 1 2 1 2 1.17 2.54", szATABLE++);
	Load1El("Se SELENIUM 34 78.96 1.17 2.55 9.752 21.2 3 2 2 1 1.17 2.54 2 1.17 2.54 4 2 3 1.17 2.54 4 1.17 2.54 6 3 3 1.17 2.54 4 1.17 2.54 6 1.17 2.54", szATABLE++);
	
	//Load1El("Br BROMINE 35 79.904 1.141 2.96 11.814 21.8 1 1 1 1 1.14 3.13", szATABLE++);
	Load1El("Br BROMINE 35 79.904 1.141 2.96 11.814 21.8 4 1 1 1 1.14 3.13 3 1 4 1.14 3.13 5 1 4 1.14 3.13 7 1 4 1.14 3.13", szATABLE++);
	
	Load1El("Kr KRYPTON 36 83.80 1.10 0.00 13.999 24.4 0", szATABLE++);
	Load1El("Rb RUBIDIUM 37 85.4678 2.35 0.82 4.177 27.3 1 1 1 1 2.35 0.82", szATABLE++);
	Load1El("Sr STRONTIUM 38 87.62 2.00 0.95 5.695 11 1 2 1 2 2.00 0.95", szATABLE++);
	Load1El("Y YTTRIUM 39 88.90585 1.80 1.22 6.38 12.2 1 3 1 3 1.80 1.22", szATABLE++);
	Load1El("Zr ZIRCONIUM 40 91.224 1.55 1.33 6.84 13.1 1 4 1 4 1.55 1.33", szATABLE++);
	Load1El("Nb NIOBIUM 41 92.90638 1.45 1.60 6.88 14.3 1 5 1 4 1.45 1.60", szATABLE++);
	Load1El("Mo MOLYBDENUM 42 95.94 1.45 2.16 7.099 16.5 1 6 1 4 1.45 2.16", szATABLE++);
	Load1El("Tc TECHNETIUM 43 98.00 1.35 1.90 7.28 15.3 3 4 1 2 1.35 1.90 6 1 6 1.35 1.9 7 1 4 1.35 1.90", szATABLE++);
	Load1El("Ru RUTHENIUM 44 101.07 1.30 2.20 7.37 16.8 1 6 1 4 1.30 2.20", szATABLE++);
	Load1El("Rh RHODIUM 45 102.90550 1.35 2.28 7.46 18.1 2 1 1 1 1.35 2.28 3 1 3 1.35 2.28", szATABLE++);
	Load1El("Pd PALLADIUM 46 106.42 1.40 2.20 8.34 19.6 1 2 1 2 1.40 2.20", szATABLE++);
	Load1El("Ag SILVER 47 107.8682 1.60 1.93 7.576 21.5 1 1 1 1 1.60 1.93", szATABLE++);
	Load1El("Cd CADMIUM 48 112.411 1.55 1.69 8.993 16.9 1 2 1 2 1.55 1.69", szATABLE++);
	Load1El("In INDIUM 49 114.818 1.55 1.78 5.786 18.9 1 3 1 3 1.55 1.78", szATABLE++);
	Load1El("Sn TIN 50 118.710 1.41 1.96 7.344 14.6 2 2 1 2 1.41 2.00 4 1 2 1.41 2.00", szATABLE++);
	Load1El("Sb ANTIMONY 51 121.760 1.41 2.05 8.641 16.5 1 3 1 3 1.41 2.19", szATABLE++);
	Load1El("Te TELLURIUM 52 127.60 1.40 2.10 9.009 18.6 1 4 1 4 1.40 2.10", szATABLE++);
	
	//Load1El("I IODINE 53 126.90447 1.333 2.66 10.451 19.1 1 1 1 1 1.33 2.93", szATABLE++);
	Load1El("I IODINE 53 126.90447 1.333 2.66 10.451 19.1 4 1 1 1 1.33 2.93 3 1 4 1.33 2.93 5 1 4 1.33 2.93 7 1 4 1.33 2.93", szATABLE++);
	
	Load1El("Xe XENON 54 131.29 1.30 0.00 12.13 21.2 0", szATABLE++);
	Load1El("Cs CESIUM 55 132.90543 2.60 0.79 3.894 25.1 1 1 1 1 2.60 0.79", szATABLE++);
	Load1El("Ba BARIUM 56 137.327 2.15 0.89 5.212 10 1 2 1 2 2.15 0.89", szATABLE++);
	Load1El("La LANTHANUM 57 138.9055 1.95 1.10 5.58 11.1 1 3 1 3 1.95 1.10", szATABLE++);
	Load1El("Ce CERIUM 58 140.115 1.65 1.12 5.47 10.8 1 3 1 3 1.65 1.12", szATABLE++);
	Load1El("Pr PRASEODYMIUM 59 140.90765 1.55 1.13 5.42 10.6 1 3 1 3 1.55 1.13", szATABLE++);
	Load1El("Nd NEODYMIUM 60 144.24 2.64 1.14 5.49 10.7 1 3 1 3 2.64 1.14", szATABLE++);
	Load1El("Pm PROMETHIUM 61 145 2.62 1.13 5.55 10.9 1 3 1 3 2.62 1.13", szATABLE++);
	Load1El("Sm SAMARIUM 62 150.36 2.59 1.17 5.63 11.1 2 2 1 2 2.59 1.17 3 1 3 2.59 1.17", szATABLE++);
	Load1El("Eu EUROPIUM 63 151.965 2.56 1.20 5.67 11.2 1 3 1 3 2.56 1.20", szATABLE++);
	Load1El("Gd GADOLINIUM 64 157.25 2.54 1.20 6.15 12.1 1 3 1 3 2.54 1.20", szATABLE++);
	Load1El("Tb TERBIUM 65 158.92534 2.51 1.20 5.86 11.5 1 3 1 3 2.51 1.20", szATABLE++);
	Load1El("Dy DYSPROSIUM 66 162.50 2.49 1.22 5.93 11.7 1 3 1 3 2.49 1.22", szATABLE++);
	Load1El("Ho HOLMIUM 67 164.93031 2.47 1.23 6.02 11.8 1 3 1 3 2.47 1.23", szATABLE++);
	Load1El("Er ERBIUM 68 167.26 2.45 1.24 6.101 11.9 1 3 1 3 2.45 1.24", szATABLE++);
	Load1El("Tm THULIUM 69 168.9342 2.42 1.25 6.184 12.1 1 3 1 3 2.42 1.25", szATABLE++);
	Load1El("Yb YTTERBIUM 70 173.04 2.40 1.10 6.254 12.2 1 3 1 3 2.40 1.10", szATABLE++);
	Load1El("Lu LUTETIUM 71 174.967 2.25 1.27 5.43 13.9 1 3 1 3 2.25 1.27", szATABLE++);
	Load1El("Hf HAFNIUM 72 178.49 2.16 1.30 6.65 14.9 1 4 1 4 2.16 1.30", szATABLE++);
	Load1El("Ta TANTALUM 73 180.9479 2.09 1.50 7.89 0 1 5 1 4 2.09 1.50", szATABLE++);
	Load1El("W TUNGSTEN 74 183.84 2.02 2.36 7.98 0 1 6 2 4 2.02 2.36 6 2.02 2.36", szATABLE++);
	Load1El("Re RHENIUM 75 186.207 1.97 1.90 7.88 0 1 7 1 4 1.97 1.90", szATABLE++);
	Load1El("Os OSMIUM 76 190.23 1.92 2.20 8.7 0 2 4 1 2 1.92 2.20 8 1 4 1.92 2.20", szATABLE++);
	Load1El("Ir IRIDIUM 77 192.217 1.87 2.20 9.1 0 3 1 1 1 1.87 2.20 3 1 3 1.87 2.2 4 1 4 1.87 2.20", szATABLE++);
	Load1El("Pt PLATINUM 78 195.08 1.83 2.28 9 18.6 2 2 1 2 1.83 2.28 4 1 4 1.83 2.28", szATABLE++);
	Load1El("Au GOLD 79 196.96654 1.79 2.54 9.225 20.5 2 1 1 1 1.79 2.54 3 1 3 1.79 2.54", szATABLE++);
	Load1El("Hg MERCURY 80 200.59 1.76 2.00 10.437 18.7 1 2 1 2 1.49 1.96", szATABLE++);
	Load1El("Tl THALLIUM 81 204.3833 2.08 2.04 6.108 20.4 2 1 1 1 2.08 2.04 3 1 3 2.08 2.04", szATABLE++);
	Load1El("Pb LEAD 82 207.2 1.81 2.33 7.416 15 2 2 1 2 1.81 2.33 4 1 4 1.81 2.33", szATABLE++);
	Load1El("Bi BISMUTH 83 208.98038 1.63 2.02 7.289 16.7 2 3 1 3 1.63 2.02 5 1 4 1.63 2.02", szATABLE++);
	Load1El("Po POLONIUM 84 209 1.53 2.00 8.42 0 1 2 1 2 1.53 2.00", szATABLE++);
	Load1El("At ASTATINE 85 210 1.43 2.20 9.5 0 2 1 1 1 1.43 2.20 5 1 3 1.43 2.20", szATABLE++);
	Load1El("Rn RADON 86 222 1.34 0.0 10.748 0 0", szATABLE++);
	Load1El("Fr FRANCIUM 87 223 2.5 0.70 4 0 1 1 1 1 2.5 0.70", szATABLE++);
	Load1El("Ra RADIUM 88 226 2.05 0.90 5.279 10.1 1 2 1 2 2.05 0.90", szATABLE++);
	Load1El("Ac ACTINIUM 89 227 1.76 1.10 5.17 12.1 1 3 1 3 1.76 1.10", szATABLE++);
	Load1El("Th THORIUM 90 232.0381 1.65 1.30 6.08 11.5 1 4 1 4 1.65 1.30", szATABLE++);
	Load1El("Pa PROTACTINIUM 91 231.03587 1.65 1.50 5.88 0 1 5 1 4 1.65 1.50", szATABLE++);
	Load1El("U URANIUM 92 238.0289 1.42 1.30 6.05 0 2 4 1 4 1.42 1.30 6 1 4 1.42 1.30", szATABLE++);
	Load1El("Np NEPTUNIUM 93 237 1.65 1.36 6.19 0 4 3 1 3 1.65 1.36 4 1 4 1.65 1.36 5 1 3 1.65 1.36 6 1 4 1.65 1.36", szATABLE++);
	Load1El("Pu PLUTONIUM 94 244 1.65 1.28 6.06 0 3 3 1 3 1.65 1.28 4 1 4 1.65 1.28 6 1 4 1.65 1.28", szATABLE++);
	Load1El("Am AMERICIUM 95 243 1.65 1.30 6 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Cm CURIUM 96 247 1.65 1.30 6.02 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Bk BERKELIUM 97 247 1.65 1.30 6.23 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Cf CALIFORNIUM 98 251 1.65 1.30 6.3 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Es EINSTEINIUM 99 252 1.65 1.30 6.42 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Fm FERMIUM 100 257 1.65 1.30 6.5 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Md MENDELEVIUM 101 258 1.65 1.30 6.58 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("No NOBELIUM 102 259 1.65 1.30 6.65 0 1 3 1 3 1.65 1.30", szATABLE++);
	Load1El("Lr LAWRENCIUM 103 262 1.65 1.30 0 0 1 3 1 3 1.65 1.30", szATABLE++);

	// elements are loaded, proceed with other thigns
	// (in the old version here was the loading from atoms.txt file
	//
	Hbonders.Dump();
	Hbonders.PutInSet(NITROGEN);
	Hbonders.PutInSet(OXYGEN);
	Hbonders.PutInSet(FLUORINE);
	Hbonders.PutInSet(CHLORINE);

	Group5Elements.Dump();
	Group5Elements.PutInSet(NITROGEN);
	Group5Elements.PutInSet(PHOSPHORUS);
	Group5Elements.PutInSet(ARSENIC);
	Group5Elements.PutInSet(ANTIMONY);
	Group5Elements.PutInSet(BISMUTH);
	

	Group6Elements.Dump();
	Group6Elements.PutInSet(OXYGEN);
	Group6Elements.PutInSet(SULFUR);
	Group6Elements.PutInSet(SELENIUM);
	Group6Elements.PutInSet(TELLURIUM);
	Group6Elements.PutInSet(POLONIUM);

	Halogens.Dump();
	Halogens.PutInSet(FLUORINE);
	Halogens.PutInSet(CHLORINE);
	Halogens.PutInSet(BROMINE);
	Halogens.PutInSet(IODINE);
	Halogens.PutInSet(ASTATINE);

	AromaticityElements.Dump();
	AromaticityElements.PutInSet(CARBON);
	AromaticityElements.PutInSet(OXYGEN);
	AromaticityElements.PutInSet(NITROGEN);
	AromaticityElements.PutInSet(SULFUR);
	AromaticityElements.PutInSet(PHOSPHORUS);		//new addition
	AromaticityElements.PutInSet(SELENIUM);		//new addition, May 21 2008
	
	OrganicElements = AromaticityElements;
	OrganicElements.PutInSet(HYDROGEN);		

	AromaticityElements.PutInSet(BORON);

	AlkaliMetals.Dump();
	AlkaliMetals.PutInSet(LITHIUM);
	AlkaliMetals.PutInSet(NATRIUM);
	AlkaliMetals.PutInSet(POTASSIUM);
	Metals.PutInSet(FRANCIUM);	
	Metals.PutInSet(RUBIDIUM);
	Metals.PutInSet(CESIUM);

	Metals = AlkaliMetals;
	//put one high number metal ahead so it will resize set	
	Metals.PutInSet(LAWRENCIUM);	
	Metals.PutInSet(BERYLLIUM);Metals.PutInSet(MAGNESIUM);
	Metals.PutInSet(ALUMINIUM);
	Metals.PutInSet(CALCIUM);Metals.PutInSet(SCANDIUM);
	Metals.PutInSet(TITANIUM);Metals.PutInSet(VANADIUM);
	Metals.PutInSet(CHROMIUM);	Metals.PutInSet(MANGANESE);
	Metals.PutInSet(IRON);Metals.PutInSet(COBALT);
	Metals.PutInSet(ZINC);Metals.PutInSet(GALLIUM);
	Metals.PutInSet(GERMANIUM);
	Metals.PutInSet(STRONTIUM);Metals.PutInSet(YTTRIUM);
	Metals.PutInSet(ZIRCONIUM);Metals.PutInSet(NIOBIUM);
	Metals.PutInSet(MOLYBDENUM);Metals.PutInSet(TECHNETIUM);
	Metals.PutInSet(RUTHENIUM);Metals.PutInSet(RHODIUM);
	Metals.PutInSet(PALLADIUM);Metals.PutInSet(SILVER);
	Metals.PutInSet(CADMIUM);	Metals.PutInSet(INDIUM);
	Metals.PutInSet(TIN);Metals.PutInSet(TELLURIUM);
	Metals.PutInSet(BARIUM);
	Metals.PutInSet(LANTHANUM);	Metals.PutInSet(CERIUM);
	Metals.PutInSet(PRASEODYMIUM);Metals.PutInSet(NEODYMIUM);
	Metals.PutInSet(PROMETHIUM);Metals.PutInSet(EUROPIUM);
	Metals.PutInSet(GADOLINIUM);Metals.PutInSet(SAMARIUM);
	Metals.PutInSet(TERBIUM);Metals.PutInSet(DYSPROSIUM);
	Metals.PutInSet(HOLMIUM);Metals.PutInSet(ERBIUM);
	Metals.PutInSet(THULIUM);Metals.PutInSet(YTTERBIUM);
	Metals.PutInSet(LUTETIUM);Metals.PutInSet(HAFNIUM);
	Metals.PutInSet(TANTALUM);Metals.PutInSet(TUNGSTEN);
	Metals.PutInSet(RHENIUM);Metals.PutInSet(OSMIUM);
	Metals.PutInSet(IRIDIUM);Metals.PutInSet(PLATINUM);
	Metals.PutInSet(LEAD);Metals.PutInSet(GOLD);
	Metals.PutInSet(MERCURY);Metals.PutInSet(THALLIUM);
	Metals.PutInSet(BISMUTH);Metals.PutInSet(POLONIUM);
	Metals.PutInSet(RADIUM);
	Metals.PutInSet(ACTINIUM);Metals.PutInSet(THORIUM);
	Metals.PutInSet(PROTACTINIUM);Metals.PutInSet(URANIUM);
	Metals.PutInSet(NEPTUNIUM);Metals.PutInSet(PLUTONIUM);
	Metals.PutInSet(AMERICIUM);Metals.PutInSet(CURIUM);
	Metals.PutInSet(BERKELIUM);Metals.PutInSet(CALIFORNIUM);
	Metals.PutInSet(EINSTEINIUM);Metals.PutInSet(FERMIUM);
	Metals.PutInSet(MENDELEVIUM);Metals.PutInSet(NOBELIUM);	

		
	Deloc.Dump();
	Deloc.PutInSet(AROMATIC_BOND);
	Deloc.PutInSet(CONJ_BOND);
	
	//note, Loc must not contain TRIPLE_BOND, because of problems with delocalization\localization of it!
	Loc.Dump();
	Loc.PutInSet(DOUBLE_BOND);
	Loc.PutInSet(SINGLE_BOND);

	PiBonds = Deloc;
	PiBonds.PutInSet(DOUBLE_BOND);
	PiBonds.PutInSet(TRIPLE_BOND);

	CovBonds = Loc | Deloc;		//covalent bonds
	CovBonds.PutInSet(TRIPLE_BOND);
	CovBonds.PutInSet(A_D_BOND);

	//various groups of aminoacid residues
	HphobResidues.Dump();
	HphobResidues.PutInSet(LEU1);
	HphobResidues.PutInSet(ILE1);
	HphobResidues.PutInSet(VAL1);
	HphobResidues.PutInSet(PHE1);
	HphobResidues.PutInSet(TRP1);
	HphobResidues.PutInSet(CYS1);
	HphobResidues.PutInSet(MET1);

	HphilResidues.Dump();
	HphilResidues.PutInSet(ARG1);
	HphilResidues.PutInSet(LYS1);
	HphilResidues.PutInSet(ASP1);
	HphilResidues.PutInSet(GLU1);
	HphilResidues.PutInSet(ASN1);
	HphilResidues.PutInSet(GLN1);

	AmbiResidues.Dump();
	AmbiResidues.PutInSet(ALA1);
	AmbiResidues.PutInSet(PRO1);
	AmbiResidues.PutInSet(GLY1);
	AmbiResidues.PutInSet(THR1);
	AmbiResidues.PutInSet(SER1);
	AmbiResidues.PutInSet(TYR1);
	AmbiResidues.PutInSet(HIS1);

	AromResidues.Dump();
	AromResidues.PutInSet(TYR1);
	AromResidues.PutInSet(HIS1);
	AromResidues.PutInSet(PHE1);
	AromResidues.PutInSet(TRP1);

	//other various settings
	SetUpColors();
	SupressMessages = false;
	StartDirName = "";

	return true;
}


void SetUpColors()
{
	UNSIGNED_2B_TYPE i;

	for (i=1; i< ATABLE.length()-1; i++)	
		ATABLE[i].COLOR = DEFAULT_COLOR;

	ATABLE[CARBON].COLOR = CARBON_COLOR;
	ATABLE[NITROGEN].COLOR = NITROGEN_COLOR;
	ATABLE[OXYGEN].COLOR = OXYGEN_COLOR;
	ATABLE[SULFUR].COLOR = SULFUR_COLOR;
	ATABLE[PHOSPHORUS].COLOR = PHOSPHORUS_COLOR;
}

void GeneratePrimeNumbers(apvector<UNSIGNED_4B_TYPE> &PN, UNSIGNED_4B_TYPE N)
{
	UNSIGNED_4B_TYPE k, l, a, b, c;
	UNSIGNED_4B_TYPE Candidat1, Candidat2, MaxD;
	
	l = 2;
	PN.resize(l);
	PN[0] = 2;
	PN[1] = 3;

	a = b = k = 0;
	while (l < N)
	{
		k++;
		Candidat1 = 6*k-1;
		Candidat2 = 6*k+1;
		MaxD = (UNSIGNED_4B_TYPE) sqrt(REALNUM_TYPE(Candidat2));
		
		for (c=2; (c<l)&&(PN[c]<=MaxD)&&(Candidat1 + Candidat2 != ZERO); c++)
		{
			if (Candidat1 != ZERO)
				a = Candidat1 / PN[c];

			if (Candidat2 != ZERO)
				b = Candidat2 / PN[c];

			if (Candidat1 == a*PN[c])
				Candidat1 = ZERO;

			if (Candidat2 == b*PN[c])
				Candidat2 = ZERO;
		};

		if (Candidat1 != ZERO)
		{
			PN.resize(l+1);
			PN[l] = Candidat1;
			l++;
		};

		if (Candidat2 != ZERO)
		{
			PN.resize(l+1);
			PN[l] = Candidat2;
			l++;
		};
	};//while l < N

}

void PutInLogFile(STRING_TYPE S)
{
	STRING_TYPE M(ModulePath);
	if (M.length()) M	+= "\\";
	M	+= LOG_FILENAME;

	FILETYPE_OUT xxfl(M.c_str(), ios_base::app | ios_base::ate); //appending
	xxfl << S << endl;
	xxfl.close();
}

void GetTimeStamp(STRING_TYPE &S)
{
	tm * st;
	time_t secs_t = time(NULL);
	st = localtime(&secs_t);
	S = asctime( st );	
	S.parse_string();
/*
	char sbuffr[256];
	SYSTEMTIME st;
	GetLocalTime(&st);
	sprintf(sbuffr, "%.2d/%.2d/%4d - %.2d:%.2d:%.2d", st.wMonth, st.wDay, st.wYear, st.wHour, st.wMinute, st.wSecond);
	S = sbuffr; S.parse_string();
*/
}


//---------------------------------------------------------
//Windows Interface
#ifdef WIN_PORT2

bool Message(STRING_TYPE S, bool OnlyWarning)
{
	PutInLogFile(S);
	if (SupressMessages)
		return false;

	if (OnlyWarning)
	{
		MessageBox(NULL, S.c_str(), "Warning", MB_OK);
		return true;
	}

	STRING_TYPE M = S + " Quit\\Ignore? (Y/N)";

	if (IDYES == MessageBox(NULL, M.c_str(), "Processing Error", MB_YESNO))
		return true;

	return false;
}

bool GetFile(STRING_TYPE &file,		//default name of file and result
			 STRING_TYPE Title,		//title of operation
			 STRING_TYPE filtername, //filter names separated by |
			 STRING_TYPE extention,	//extention names separated by |
			 HWND hwnd)
{
	char FileName[1024], filters[DEF_STRING_SIZE];
	OPENFILENAME ofn;		
	SIGNED_4B_TYPE i, i1, j = ZERO, j1 = ZERO, result = ZERO;
	STRING_TYPE ws;

	CutStrEnding(file);	//12.09.2008 to prevent conflicts
	sprintf(FileName, "%s", file.c_str());
	FileName[file.length()] = char(0);
	
	strcpy(filters, "");
	do 
	{
		i = filtername.find('|');
		i1 = extention.find('|');
		if (i == INVALID) 	i = filtername.length(); else filtername[i] = '!';
		if (i1 == INVALID) 	i1 = extention.length(); else extention[i1] = '!';
		
		ws = filtername.substr(j,i-j);

		strcat(filters, ws.c_str());
		strcat(filters," (");
		
		ws = extention.substr(j1, i1-j1);
		
		strcat(filters, ws.c_str());
		strcat(filters,")!");
		strcat(filters, ws.c_str());
		strcat(filters, "!");		

		j = i+1;
		j1= i1+1;
	} while ((i != filtername.length()) || (i1 != extention.length()));

	while (strrchr(filters, '!') != NULL)
		*strrchr(filters, '!') = '\0';
	
	memset(&ofn, 0, sizeof(OPENFILENAME));
	ofn.lStructSize = sizeof(OPENFILENAME);
	ofn.hwndOwner = hwnd;
	ofn.nFilterIndex = 1;
	ofn.lpstrFile = (char *)&FileName;
	ofn.nMaxFile = 1024;

	if (extention.find(".*") >= 0)
		ofn.lpstrDefExt = NULL;
	else
	{
		if (strrchr(extention.c_str(), '.') != NULL)
			ofn.lpstrDefExt = strrchr(extention.c_str(), '.')+1;
		else
			ofn.lpstrDefExt = extention.c_str();
	};

	ofn.lpstrCustomFilter = NULL;
	ofn.lpstrFilter = (char *)&filters;	
	ofn.lpstrTitle = Title.c_str();	
	ofn.lpstrFileTitle = NULL;
	ofn.Flags = OFN_HIDEREADONLY | OFN_NOREADONLYRETURN | OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
	ofn.lpstrInitialDir	= StartDirName.c_str();

	result = GetSaveFileName(&ofn);

	file = FileName;

	if ((ofn.nFileOffset > 0) && (ofn.nFileOffset < strlen(FileName)))
	//change Start Directory		
		StartDirName  = file.substr(0, ofn.nFileOffset);		
	
	return (result != 0);
}

bool PutFile(STRING_TYPE &file,		//default name of file and result
			 STRING_TYPE Title,		//title of operation
			 STRING_TYPE filtername, //filter names separated by '|'
			 STRING_TYPE extention,	//extention names separated by '|'
			 HWND hwnd, bool ifMltFiles)
{
	char FileName[1024], filters[DEF_STRING_SIZE];
	STRING_TYPE ws;
	OPENFILENAME ofn;
	SIGNED_4B_TYPE result = ZERO, j = ZERO, j1 = ZERO, i, i1;

	sprintf(FileName, "%s", file.c_str());
	FileName[file.length()] = char(0);
	
	strcpy(filters, "");
	do 
	{
		i = filtername.find('|');
		i1 = extention.find('|');
		if (i == INVALID) 	i = filtername.length(); else filtername[i] = '!';
		if (i1 == INVALID) 	i1 = extention.length(); else extention[i1] = '!';
		
		ws = filtername.substr(j,i-j);

		strcat(filters, ws.c_str());
		strcat(filters," (");
		
		ws = extention.substr(j1, i1-j1);
		
		strcat(filters, ws.c_str());
		strcat(filters,")!");
		strcat(filters, ws.c_str());
		strcat(filters, "!");		

		j = i+1;
		j1= i1+1;
	} while ((i != filtername.length()) || (i1 != extention.length()));

	while (strrchr(filters, '!') != NULL)
		*strrchr(filters, '!') = '\0';
	
	memset(&ofn, ZERO, sizeof(OPENFILENAME));
	ofn.lStructSize = sizeof(OPENFILENAME);
	ofn.hwndOwner = hwnd;
	ofn.nFilterIndex = 1;
	ofn.lpstrFile = (char *)&FileName;
	ofn.nMaxFile = 1024;

	if (extention.find(".*") >= ZERO)
		ofn.lpstrDefExt = NULL;
	else
	{
		if (strrchr(extention.c_str(), '.') != NULL)
			ofn.lpstrDefExt = strrchr(extention.c_str(), '.')+1;
		else
			ofn.lpstrDefExt = extention.c_str();
	};
	

	ofn.lpstrCustomFilter = NULL;
	ofn.lpstrFilter = (char *)&filters;
	
	ofn.lpstrTitle = Title.c_str();	
	ofn.lpstrFileTitle = NULL;
	ofn.lpstrInitialDir	= StartDirName.c_str();
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	if (ifMltFiles)	ofn.Flags |= OFN_ALLOWMULTISELECT | OFN_EXPLORER;
	result = GetOpenFileName(&ofn);

	if (ifMltFiles)
	{//make names tab-separated
		i = 0;
		for (;;)
		{
			if (FileName[i] == '\0') 
			{
				FileName[i] = '\t';
				if (FileName[i+1] == '\0') break;
			}
			i++;
		}
	}

	file = FileName;

	if ((ofn.nFileOffset > 0) && (ofn.nFileOffset < strlen(FileName)))
	//change Start Directory		
		StartDirName  = file.substr(0, ofn.nFileOffset);		
		
	return (result != 0);
}


#endif	//#ifdef WIN_PORT2
//---------------------------------------------------------

bool CheckStrEnding(STRING_TYPE &A, STRING_TYPE ends) 
//can unambiguously check file-types by means of file-extension masks
{
	STRING_TYPE TStr = A.substr( A.length() - ends.length(), ends.length() );
	return ( TStr.find(ends) == 0 );
}

void CutStrEnding(STRING_TYPE &S)
{//cuts the file-extention '.xxx' off the string's end
	char * v = (char *)strrchr(S.c_str(), '.');
	if ( v  > strrchr(S.c_str(), '\\') )
		S = S.substr( 0, SIGNED_4B_TYPE(v - S.c_str()) );
}

REALNUM_TYPE String2Number(STRING_TYPE &L)
//creates a string-representing number. 
//Similar numbers correspond to similar strings (not like CRC32!)
{
	REALNUM_TYPE x = ZERO, ch, div = 1;	
	for (SIGNED_4B_TYPE sc = ZERO; sc < L.length(); sc++)
	{
		ch = UNSIGNED_1B_TYPE(L[sc]);
		ch /= div;
		x += ch;
		div *= 256;
	}
	return x;
}

void SplitString(STRING_TYPE &BASE, STRING_TYPE SEP, apvector<STRING_TYPE> &SPLIT)
{
	SPLIT.resize(0);
	if ( (BASE.length() < SEP.length()) || (SEP.length() == 0) ) return;

	STRING_TYPE B = BASE + SEP;
	SIGNED_4B_TYPE i =0, cx =0, ix =B.find(SEP);
	SPLIT.resize(B.length() / SEP.length());
	do 
	{
		SPLIT[i++] = B.substr(cx, ix - cx);
		cx = ix + SEP.length();
		B[ix]++; //to destroy SEP-match
		ix = B.find(SEP);
	} while (ix >= 0);
	SPLIT.resize(i);
}

SIGNED_4B_TYPE LoadSetAsText(FILETYPE_IN &fip, set &setData, SIGNED_4B_TYPE sh)
{//fip should point in a text file to a string of numbers: set_size element1 element2 ...
	setData.Dump();
	SIGNED_4B_TYPE uu, zz, ii;
	fip >> uu;
	for (ii = 0; ii < uu; ii++) {	fip >> zz;	setData.PutInSet(zz - sh); };
	return uu;
}

SIGNED_4B_TYPE SaveSetAsText(FILETYPE_OUT &fo, set &setData, SIGNED_4B_TYPE sh)
{
	SIGNED_4B_TYPE ii;
	apvector<SIGNED_4B_TYPE> els;
	setData.GetList(els);
	fo << els.length();
	for (ii = 0; ii < els.length(); ii++)	fo << " " << els[ii]+ sh;
	fo << endl;
	return els.length();
}

UNSIGNED_4B_TYPE AtomTypeHash::hashFunction(UNSIGNED_4B_TYPE Fingerprint)
{
	return (Fingerprint);
}
/////





//----------------------------------------------------------//
//--------------- ring handling subroutines ----------------//
//----------------------------------------------------------//

void CleanRing(Ring &r)
{//clean all structures of Ring
	r.done  = ZERO;
	r.lRing.Dump();
	r.setRing.Dump();
	r.pRing.resize(ZERO);
	r.flat  = r.arom  = 
	r.dlcl  = r.pi_s  = false;
}

void InitializeRing(Ring &r)
{	
	r.setRing.Dump();	
	r.done = 0;
	r.lRing.Dump();

	for (SIGNED_4B_TYPE n=0; n<r.pRing.length(); n++)
	{			
			r.setRing.PutInSet(r.pRing[n]);	
			r.lRing.Insert(r.pRing[n]);
			r.lRing.Next();
	}//for n
}

void InitializeRings(apvector<Ring> &rgs)
{	
	SIGNED_4B_TYPE rC = rgs.length();
	for (SIGNED_4B_TYPE m = 0; m<rC; m++)
		InitializeRing(rgs[m]);	
}

void FindRingSystems(apvector<Ring> &rr, apvector<SIGNED_4B_TYPE> &ss, UNSIGNED_1B_TYPE MinIntersection)
//description:		assigns rings from rr[] to several ring systems
//					for that array ss[] is resized to the number of rings in rr[]
//					then ss[i] represents ID of a ring system, to which rr[i] belongs.
//
//MinIntersection:	minimum number of atoms, which two rings should share in order to be merged.
//					default is 2, make it 1, if spiro-connected rings should be in one system
//
//precondition:		rr[] must be initialized by InitializeRings
{
	SIGNED_4B_TYPE m, n, rC = rr.length();
	bool Chan;
	apvector<set> sysets;
	set tS;

	sysets.resize(rC);
	ss.resize(rC);
	for (m = ZERO; m < rC; m++)
	{
		ss[m] = m;
		sysets[m] = rr[m].setRing;
	}
		
	do 
	{
		Chan= false;
		
		for (m = ZERO; m < rC; m++)			
		{
			if (sysets[m].IsEmpty()) 
				continue;

			for (n = m + 1; n < rC; n++)
			{
				if (sysets[n].IsEmpty()) 
					continue;

				if (ss[m] == ss[n]) 
					continue;

				tS = sysets[m] & sysets[n];
				if (tS.Size() < MinIntersection) 
					continue;

				sysets[m] |= sysets[n];
				sysets[n].Dump();
				ss[n] = ss[m];
				
				Chan = true;					
			}//for n
		} //for m
		
	} while (Chan);

	//now, clean up leftovers 
	do 
	{
		Chan = false;
		for (n = 0; n < rC; n++)
			if ((sysets[n].IsEmpty()) && (ss[n] != ss[ss[n]]))
			{
				ss[n] = ss[ss[n]];
				Chan = true;
			}
	}while (Chan);

	//now ss contains assignment of the rings to their ring systems
}

void FinalizeRings(apvector<Ring> &rngs, apvector<SIGNED_4B_TYPE> &ss, SIGNED_4B_TYPE N_MOL)
//description:			minimizes given set of rings, i.e finds
//						a smallest group of smallest rings that still
//						fully characterize the molecule.
//
//N_MOL:				molecule's number of atoms
//
//precondition:			prior to calling this function, rngs[] info 
//						must be filled by InitializeRings() and FindRingSystems()
//						rngs[] should not have internal edges in them:
//						such edge is a bond between ring's atoms that does not belong to the ring
{
	apvector<CEdge> Edges;	//pool of edges in a ring system	
	apvector<set> Ats;		//sets of edges for each atom
	apvector<SIGNED_4B_TYPE> tmp, cuRsys, satoms, smr_stack(DEFAULT_LARGE_CYCLE);
	apvector<set> smr_stack_s(DEFAULT_LARGE_CYCLE);
	
	set alfa, lried, smr;	//work sets
	UNSIGNED_2B_TYPE crsysN, rsysN, rs;
	SIGNED_4B_TYPE m, n, l, newrC, rC = rngs.length(), aC, smrl, ixa;
	bool rfound, cfound;
	
	//first, extract all the stored ring systems ranks
	for (m = ZERO; m < rC; m++)
		alfa.PutInSet(ss[m]);
	
	alfa.GetList(tmp);
	rsysN = tmp.length();
	for (rs = ZERO; rs < rsysN; rs++)
		for (m = ZERO; m < rC; m++)
			if (ss[m] == tmp[rs])
				//explicitly rewrite ranks into the simple increasing order
				ss[m] = rs;	
	//-------------------------------------------------


/*
if (SupressMessages)
{
	STRING_TYPE S1("Number of ring systems detected: "),S2 = rsysN;
	Message(S1+S2,true);
}
//*/

	newrC = rC;
	for (rs = ZERO; rs < rsysN; rs++)
	{		
		alfa.Dump();	//store all the rings of the ring system rs
		for (m = ZERO; m < rC; m++)
			if (ss[m] == rs)
				alfa.PutInSet(m);
		alfa.GetList(cuRsys);
	
		crsysN = cuRsys.length();
		if (crsysN == 1)//only one ring in the ring system, skip it.		
			continue;


		
		MinimizeRings(rngs, cuRsys);

		

		//get all the atoms of the ring system
		alfa.Dump();
		for (m = ZERO; m<crsysN; m++)
			for (l = ZERO; l < rngs[cuRsys[m]].pRing.length(); l++)
				alfa.PutInSet(rngs[cuRsys[m]].pRing[l]);
		
		alfa.GetList(satoms);	//list of atoms in the current ring system
		aC = satoms.length();	//# number of atoms in the current ring system
		//---------------------------------------------------
		


		//---------------------------------------------------------		
		Edges.resize(aC + crsysN);//estimate of # of edges in the current ryng system. +1 is added for the currently analized edge
		Ats.resize(ZERO);
		Ats.resize(N_MOL);		  //info about the edges pertaining to each atom
		l = ZERO;
		//now, generate all rings' edges
		for (m = ZERO; m < crsysN; m++)
		{
			n = rngs[cuRsys[m]].pRing.length();			
			while (n > ZERO)
			{
				n--;
				Edges[l].A1 = rngs[cuRsys[m]].lRing.Curr();
				Edges[l].A2 = rngs[cuRsys[m]].lRing.Next();				

				alfa = Ats[Edges[l].A1] & Ats[Edges[l].A2];
				if (alfa.IsEmpty())
				{
					Ats[Edges[l].A1].PutInSet(l);
					Ats[Edges[l].A2].PutInSet(l);
					l++;
				}
			}
		}// for m
		
		if (l < Edges.length()) //to remove last, temporary, working edge, if necessary
			Edges.resize(l);
		//---------------------------------------------------------
				

		//let us find all small rings in the current ring system		
		lried.Dump();
		for (l = ZERO; l < aC - 1; l++)
		{//no need to explore last atom
			smr.Dump();
			smr.PutInSet(satoms[l]);
			smr_stack[ZERO] = satoms[l];
			smr_stack_s[ZERO] = Ats[satoms[l]];
			smrl = 1;			
				
			while (smrl)
			{
				smr_stack_s[smrl - 1].GetList(tmp);
				rfound = cfound	= false;
				for (n = ZERO; n < tmp.length(); n++)
				{
					if (SIGNED_4B_TYPE(Edges[tmp[n]].A1) == smr_stack[smrl - 1])
						ixa = Edges[tmp[n]].A2;
					else
						ixa = Edges[tmp[n]].A1;

					if (lried.IsInSet(ixa))
						continue;
					
					if (smr.IsInSet(ixa))
					{//ring is found				
						if (smr_stack[smrl - 2] == ixa)
							continue;

						cfound = false;
						
						//check if no internal edges
						if (ixa == satoms[l]) //true ring
							rfound = true;
						else
						{
							rfound = false;	  //needed to prevent 1-edge bridged rings
							break;
						}
						
						continue;					
					};
					cfound = true;
					smr_stack[smrl] = ixa;
					m = tmp[n];
				}//for n

				if (rfound)
				{
					cfound = false;
					//check if it already exists (enough if by vertices, 
					//because of we avoid the cycles with single edge bridges
					//NOTE: we assume the rings in rngs[] do not have internal edges
					//(a bond between ring's atoms that does not belong to the ring)					
					
					//NB(when crsysN was of UNSIGNED_4B_TYPE, it induced
					//type conversion compiler error below, but only in release mode)
					for (m = ZERO; m < crsysN; m++)
						if (rngs[cuRsys[m]].pRing.length() == smrl)
						{
							alfa = smr | rngs[cuRsys[m]].setRing;
							if (alfa.Size() == smrl) //same					
								break;
						}
					
					if (m == crsysN) //new ring
					{
						if (rngs.length() == newrC)
							rngs.resize(newrC + crsysN);

						rngs[newrC].pRing	= smr_stack;
						rngs[newrC].pRing.resize(smrl);
						InitializeRing(rngs[newrC]);

						cuRsys.resize(++crsysN);
						cuRsys[crsysN - 1] = newrC++;
					}
				}

				if (cfound && (smrl < DEFAULT_LARGE_CYCLE - 1))
				{
					smr_stack_s[smrl - 1].RemoveFromSet(m);	//remove this edge

					smr.PutInSet(smr_stack[smrl]);
					smr_stack_s[smrl] = Ats[smr_stack[smrl]];	//copy fresh edges
					smrl++;
					
				}
				else
				{//backtrack
					smrl--;
					smr.RemoveFromSet(smr_stack[smrl]);
				}
			}//while smrl stack, looking for all small rings in the ring system

			lried.PutInSet(satoms[l]);	//remove from further searches
		};//for l		


		//-------------------  sort all rings by size, in ascending order
		do 
		{
			cfound = false;

			for (m = ZERO; m < crsysN - 1; m++)
				if (rngs[cuRsys[m]].pRing.length() > rngs[cuRsys[m + 1]].pRing.length())
				{
					cfound = true;
					Swap(cuRsys[m], cuRsys[m+1], n);
				}

		} while (cfound);
		//----------------------------------------

		//prospective improvements:
		//to remove redundant large rings if any

/*
	if (SuppressMessages)
	{
		STRING_TYPE S1("ring system N"),S2 = rs, S3 = crsysN;
		S1 += S2; 	S1 += ", "; 	S1 += S3; 	S1 += " rings;";
		Message(S1, true);

		for (m = ZERO; m < crsysN; m++)
		{		
			S2  = cuRsys[m];
			S3  = rngs[cuRsys[m]].pRing.length();
			S1  = "ring N"; S1 += S2; 
			S1 += ", length = "; 	S1 += S3; 	S1 += ", atoms:";
			
			for (n = ZERO; n < rngs[cuRsys[m]].pRing.length(); n++)
			{
				S2 = rngs[cuRsys[m]].pRing[n]);
				S1 += " ";
				S1 += S2;
			};
			Message(S1, true);
		}
	}
*/
		
	}//for rs;
	rngs.resize(newrC);
}//FinalizeRings()



void MinimizeRings(apvector<Ring> &rings, apvector<SIGNED_4B_TYPE> &Atomlist)
//description:			minimizes given set of rings, specified by rings[] 
//						and Atomlist[]. The latter usually should point to
//						subset of rings that constitute a single ring-system
//
//precondition:			prior to call of this function, rings[] info 
//						must be filled by InitializeRings()
//
//NOTE:					Normally, should be called from FinalizeRings() only.
{
	SIGNED_4B_TYPE m, n, l, currI, rC = Atomlist.length();
	SIGNED_4B_TYPE OldRingSize, IntersectionSize;
	listItem<SIGNED_4B_TYPE> *pointA, *pointB;
	set Rset;
	UNSIGNED_4B_TYPE lng;
	bool Change;
	
	do 	{
			Change = false;		

			//minimize rings
			for (m=0; m<rC; m++)
				for (n=0; n< rC; n++)
				if ( (n != m) && (rings[Atomlist[m]].pRing.length() <= rings[Atomlist[n]].pRing.length()) )
				{
					//intersect
					Rset.Dump();
					Rset = rings[Atomlist[m]].setRing & rings[Atomlist[n]].setRing;
					IntersectionSize = Rset.Size();

					if (2*IntersectionSize <= rings[Atomlist[m]].setRing.Size())
					//no sufficient overlap
						continue;
					

					////////////////////////////////////////////////////////
					//check the integrity of intersection
					//it should be a simple connected chain without discontinuities
					//////
					//1. ring n (larger one!)
					//move on the atom which belongs to intersection
					while (!Rset.IsInSet(rings[Atomlist[n]].lRing.Next()));					
					//find in the ring's terminal atom of the intersection
					while (Rset.IsInSet(rings[Atomlist[n]].lRing.Next()));
					currI = rings[Atomlist[n]].lRing.Prev();
					
					//synchronize m's ring current pointer accordingly
					while (rings[Atomlist[m]].lRing.Next() != currI);
					//compare intersection atoms step-by-step
					lng = 1;
					while (rings[Atomlist[n]].lRing.Prev() 
							== rings[Atomlist[m]].lRing.Prev())	lng++;
					if (lng == 1)
					{
						rings[Atomlist[n]].lRing.Next();
						rings[Atomlist[m]].lRing.Next();
						while (rings[Atomlist[n]].lRing.Prev() 
						== rings[Atomlist[m]].lRing.Next())	lng++;
					}
					
					if (SIGNED_4B_TYPE(lng) != IntersectionSize) //improper intersection
						continue;
					////////////////////////////////////////////////////////


							
					if (IntersectionSize != rings[Atomlist[m]].pRing.length())
					{
							currI = rings[Atomlist[m]].lRing.Curr();

						//move on atom which belongs only to ring m
						while (rings[Atomlist[n]].setRing.IsInSet(currI))
							currI = rings[Atomlist[m]].lRing.Next();

						pointA = rings[Atomlist[m]].lRing.GetCurItem();
					}
					else					
					{
						currI = rings[Atomlist[n]].lRing.Curr();

						//move on atom which belongs only to ring n
						while (rings[Atomlist[m]].setRing.IsInSet(currI))
							currI = rings[Atomlist[n]].lRing.Next();

						pointA = rings[Atomlist[n]].lRing.GetCurItem();
					};

					//find first joint
					while ((!rings[Atomlist[n]].setRing.IsInSet(pointA->index)) || (!rings[Atomlist[m]].setRing.IsInSet(pointA->index)))
						pointA = pointA->fw_ptr;

					pointB = pointA;

					//find second joint
					while (rings[Atomlist[n]].setRing.IsInSet(pointB->index) && rings[Atomlist[m]].setRing.IsInSet(pointB->index))
						pointB = pointB->fw_ptr; 

					pointB = pointB->bk_ptr;

					
					currI = rings[Atomlist[n]].lRing.Curr();
					while (currI != pointA->index)						
						currI = rings[Atomlist[n]].lRing.Next();

					OldRingSize = rings[Atomlist[n]].pRing.length();
					rings[Atomlist[n]].pRing.resize(1);
					rings[Atomlist[n]].pRing[0] = currI;

					//copy elements that belong to ring n only (between pointA and pointB)
					if (!rings[Atomlist[m]].setRing.IsInSet(rings[Atomlist[n]].lRing.Next()))
						//means pointA->ptr points toward atoms, which belong to ring n only.
					{
						rings[Atomlist[n]].lRing.Prev();
						do
						{
							currI = rings[Atomlist[n]].lRing.Next();
							lng = rings[Atomlist[n]].pRing.length();
							rings[Atomlist[n]].pRing.resize(lng+1);
							rings[Atomlist[n]].pRing[lng] = currI;
						} while(!rings[Atomlist[m]].setRing.IsInSet(currI));
					}
					else						
					{//change direction
						rings[Atomlist[n]].lRing.Prev();
						do
						{
							currI = rings[Atomlist[n]].lRing.Prev();
							lng = rings[Atomlist[n]].pRing.length();
							rings[Atomlist[n]].pRing.resize(lng+1);
							rings[Atomlist[n]].pRing[lng] = currI;

						} while(!rings[Atomlist[m]].setRing.IsInSet(currI));
					};
					 					
					if (rings[Atomlist[m]].pRing.length() != IntersectionSize)
					{
						//copy sequence from ring m (from marginB to marginA)
						currI = rings[Atomlist[m]].lRing.Curr();
						while (currI != pointB->index)
							currI = rings[Atomlist[m]].lRing.Next();
						
						if (!rings[Atomlist[n]].setRing.IsInSet(rings[Atomlist[m]].lRing.Next()))
						{							
							currI = rings[Atomlist[m]].lRing.Curr();
							do
							{
								lng = rings[Atomlist[n]].pRing.length();
								rings[Atomlist[n]].pRing.resize(lng+1);
								rings[Atomlist[n]].pRing[lng] = currI;								

								currI = rings[Atomlist[m]].lRing.Next();
							}while (currI != pointA->index);
						}
						else
						{//change direction
							rings[Atomlist[m]].lRing.Prev();
							currI = rings[Atomlist[m]].lRing.Prev();
							do
							{				
								lng = rings[Atomlist[n]].pRing.length();
								rings[Atomlist[n]].pRing.resize(lng+1);
								rings[Atomlist[n]].pRing[lng] = currI;								

								currI = rings[Atomlist[m]].lRing.Prev();							
							}while (currI != pointA->index);

						};						
						
					};//if (rings[Atomlist[m]].pRing.length() != IntersectionSize)

					if (OldRingSize <= rings[Atomlist[n]].pRing.length())
					{
						//if we gain nothing then put old ring back
						rings[Atomlist[n]].pRing.resize(OldRingSize);
						currI = rings[Atomlist[n]].lRing.Curr();
						for(l=0; l<OldRingSize; l++)
						{
							rings[Atomlist[n]].pRing[l] = currI;
							currI= rings[Atomlist[n]].lRing.Next();
						};
						continue;						
					};
					
					Change = true;							

					//now we have new ring array in rings[Atomlist[n]].pRing, 
					//we should reinitilize auxiliary structures
					InitializeRing(rings[Atomlist[n]]);
					
					//------- end of reorganizing the ring n --------
					
				};//for m,n

			}while (Change);	
}

UNSIGNED_4B_TYPE randomNumber(UNSIGNED_4B_TYPE limit)
{// Assumes limit <= RAND_MAX
    UNSIGNED_4B_TYPE n;
    UNSIGNED_4B_TYPE mask = 0xFFFFFFFF;		
    while (mask > (limit << 1)) mask >>= 1;	//remove unneeded bits

    do {
        n = rand();
        n &= mask;
    } while (n >= limit);

    return n;
}

UNSIGNED_4B_TYPE GetRandomNumber(UNSIGNED_4B_TYPE limit)
//the number generated will be within [0 .. limit - 1]; limit should be < 0x3FFFFFFF
//srand() have to be called earlier to provide unrepetitive randomization
//NB! RangeSize can be bigger than RAND_MAX!
//RAND_MAX on Windows is 0x7FFF and on LINUX 0x7FFFFFFF
{
	UNSIGNED_4B_TYPE n, b;
    UNSIGNED_4B_TYPE mask = 0xFFFFFFFF;
	
	if (limit == RAND_MAX) return rand();
	
	if (limit < RAND_MAX)
	{
		while (mask > (limit << 1)) mask >>= 1;	//remove unneeded bits
		for(;;)
		{
			n  = rand();
			n &= mask;
			if (limit > n) return n;
		};
	}

	//below are the handling mostly for Windows platform
	n = RAND_MAX; 
	b = 0;
	mask = 1;
	while (limit > n) { n <<= 1; n++; b++; mask <<= 1; mask++; };
	
	for(;;)
	{
		n = rand();
		n <<= b;
		n += (rand() & mask);
		if (limit > n) return n;
	}
}

void SortRandomly (apvector<REALNUM_TYPE> &V)
{
	SIGNED_4B_TYPE N = V.length();
	SIGNED_4B_TYPE *A	= GRAB_MEM_BLOCKS(SIGNED_4B_TYPE, N), ic = ZERO;
	REALNUM_TYPE *F		= GRAB_MEM_BLOCKS(REALNUM_TYPE, N);
	QSortScore = F;
	for (; ic < N; ic++)
	{
		A[ic] = ic;	
		F[ic] = (REALNUM_TYPE)rand(); 
	}
	qsort(A, (size_t)N, sizeof(SIGNED_4B_TYPE), QSortCompareGreater);
	DROP_MEM_BLOCKS(F);
	QSortScore = NULL;

	apvector<REALNUM_TYPE> VNEW(N);
	for (ic = ZERO; ic < N; ic++)	VNEW[ic] = V[A[ic]];
	DROP_MEM_BLOCKS(A);
	V = VNEW;
}

//QSort paraphernalia
REALNUM_TYPE *QSortScore = NULL;
int QSortCompareGreater(const void *arg1, const void *arg2)
{   
   if (QSortScore == NULL)
		return ZERO;

   if (QSortScore[* ( UNSIGNED_4B_TYPE* ) arg1] >  QSortScore[* ( UNSIGNED_4B_TYPE* ) arg2])
		return 1;
   
   return -1;
}

int QSortCompareLess(const void *arg1, const void *arg2)
{   
   if (QSortScore == NULL)
		return ZERO;

   if (QSortScore[* ( UNSIGNED_4B_TYPE* ) arg1] <  QSortScore[* ( UNSIGNED_4B_TYPE* ) arg2])
		return 1;
   
   return -1;
}

void BubbleSort(apvector<REALNUM_TYPE> &x, apvector<SIGNED_4B_TYPE> &a)
//description:	sorts array "a" in descending order on the base of "x" ("x" gets sorted too)
//				this function is needed for internal purposes of DefineStereoisomers() and get3d()
//
//Note:			arrays x and a should have same length
{
	bool SwapFlag;
	SIGNED_4B_TYPE j;
	
	if (x.length() != a.length())
		return;

	SIGNED_4B_TYPE tm;
	REALNUM_TYPE   tmr;
	do
	{
		SwapFlag = false;
		for (j=0; j<x.length()-1; j++)
			if (x[j] < x[j+1])
			{
				Swap(x[j], x[j+1], tmr);
				Swap(a[j], a[j+1], tm);
				SwapFlag = true;
			};
	} while (SwapFlag);
}


UNSIGNED_4B_TYPE CRC32Table[256] = 
{ 0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA,
  0x076DC419, 0x706AF48F, 0xE963A535, 0x9E6495A3,
  0x0EDB8832, 0x79DCB8A4, 0xE0D5E91E, 0x97D2D988,
  0x09B64C2B, 0x7EB17CBD, 0xE7B82D07, 0x90BF1D91,
  0x1DB71064, 0x6AB020F2, 0xF3B97148, 0x84BE41DE,
  0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
  0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC,
  0x14015C4F, 0x63066CD9, 0xFA0F3D63, 0x8D080DF5,
  0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172,
  0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B,
  0x35B5A8FA, 0x42B2986C, 0xDBBBC9D6, 0xACBCF940,
  0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59,
  0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116,
  0x21B4F4B5, 0x56B3C423, 0xCFBA9599, 0xB8BDA50F,
  0x2802B89E, 0x5F058808, 0xC60CD9B2, 0xB10BE924,
  0x2F6F7C87, 0x58684C11, 0xC1611DAB, 0xB6662D3D,

  0x76DC4190, 0x01DB7106, 0x98D220BC, 0xEFD5102A,
  0x71B18589, 0x06B6B51F, 0x9FBFE4A5, 0xE8B8D433,
  0x7807C9A2, 0x0F00F934, 0x9609A88E, 0xE10E9818,
  0x7F6A0DBB, 0x086D3D2D, 0x91646C97, 0xE6635C01,
  0x6B6B51F4, 0x1C6C6162, 0x856530D8, 0xF262004E,
  0x6C0695ED, 0x1B01A57B, 0x8208F4C1, 0xF50FC457,
  0x65B0D9C6, 0x12B7E950, 0x8BBEB8EA, 0xFCB9887C,
  0x62DD1DDF, 0x15DA2D49, 0x8CD37CF3, 0xFBD44C65,
  0x4DB26158, 0x3AB551CE, 0xA3BC0074, 0xD4BB30E2,
  0x4ADFA541, 0x3DD895D7, 0xA4D1C46D, 0xD3D6F4FB,
  0x4369E96A, 0x346ED9FC, 0xAD678846, 0xDA60B8D0,
  0x44042D73, 0x33031DE5, 0xAA0A4C5F, 0xDD0D7CC9,
  0x5005713C, 0x270241AA, 0xBE0B1010, 0xC90C2086,
  0x5768B525, 0x206F85B3, 0xB966D409, 0xCE61E49F,
  0x5EDEF90E, 0x29D9C998, 0xB0D09822, 0xC7D7A8B4,
  0x59B33D17, 0x2EB40D81, 0xB7BD5C3B, 0xC0BA6CAD,

  0xEDB88320, 0x9ABFB3B6, 0x03B6E20C, 0x74B1D29A,
  0xEAD54739, 0x9DD277AF, 0x04DB2615, 0x73DC1683,
  0xE3630B12, 0x94643B84, 0x0D6D6A3E, 0x7A6A5AA8,
  0xE40ECF0B, 0x9309FF9D, 0x0A00AE27, 0x7D079EB1,
  0xF00F9344, 0x8708A3D2, 0x1E01F268, 0x6906C2FE,
  0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7,
  0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC,
  0xF9B9DF6F, 0x8EBEEFF9, 0x17B7BE43, 0x60B08ED5,
  0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252,
  0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B,
  0xD80D2BDA, 0xAF0A1B4C, 0x36034AF6, 0x41047A60,
  0xDF60EFC3, 0xA867DF55, 0x316E8EEF, 0x4669BE79,
  0xCB61B38C, 0xBC66831A, 0x256FD2A0, 0x5268E236,
  0xCC0C7795, 0xBB0B4703, 0x220216B9, 0x5505262F,
  0xC5BA3BBE, 0xB2BD0B28, 0x2BB45A92, 0x5CB36A04,
  0xC2D7FFA7, 0xB5D0CF31, 0x2CD99E8B, 0x5BDEAE1D,

  0x9B64C2B0, 0xEC63F226, 0x756AA39C, 0x026D930A,
  0x9C0906A9, 0xEB0E363F, 0x72076785, 0x05005713,
  0x95BF4A82, 0xE2B87A14, 0x7BB12BAE, 0x0CB61B38,
  0x92D28E9B, 0xE5D5BE0D, 0x7CDCEFB7, 0x0BDBDF21,
  0x86D3D2D4, 0xF1D4E242, 0x68DDB3F8, 0x1FDA836E,
  0x81BE16CD, 0xF6B9265B, 0x6FB077E1, 0x18B74777,
  0x88085AE6, 0xFF0F6A70, 0x66063BCA, 0x11010B5C,
  0x8F659EFF, 0xF862AE69, 0x616BFFD3, 0x166CCF45,
  0xA00AE278, 0xD70DD2EE, 0x4E048354, 0x3903B3C2,
  0xA7672661, 0xD06016F7, 0x4969474D, 0x3E6E77DB,
  0xAED16A4A, 0xD9D65ADC, 0x40DF0B66, 0x37D83BF0,
  0xA9BCAE53, 0xDEBB9EC5, 0x47B2CF7F, 0x30B5FFE9,
  0xBDBDF21C, 0xCABAC28A, 0x53B39330, 0x24B4A3A6,
  0xBAD03605, 0xCDD70693, 0x54DE5729, 0x23D967BF,
  0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94,
  0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D};

UNSIGNED_4B_TYPE GetCRC32 (const GENERIC_POINTER p, UNSIGNED_4B_TYPE len)
//For each byte:
// 1.	XOR  the input byte with the low-order byte of
//		the CRC register to get a table index t
//
// 2.	Shift the CRC register eight bits to the right
//
// 3.	XOR the CRC register with the contents of CRC32Table[t]
{
  
  UNSIGNED_4B_TYPE CRCr = 0xFFFFFFFF;
  UNSIGNED_1B_TYPE * q = (UNSIGNED_1B_TYPE *)p, t;
  
  for (UNSIGNED_4B_TYPE i = ZERO;i < len; i++)
  {
	t  = *q++;
	t ^= (CRCr & 0xFF);
    CRCr >>= 8;
	CRCr ^= CRC32Table[t];    
  }

  return CRCr;
} //CalcCRC32


SIGNED_4B_TYPE FindArrPoz(apvector<SIGNED_4B_TYPE> &ARR, SIGNED_4B_TYPE V)
//description:		returns position in KEYS where hashkey V should be
//perconditon:		array must be sorted!
//postcondition:	returns found position or -1 otherwise
{
	UNSIGNED_4B_TYPE a = 0, e = ARR.length(), m;
	while (e > a + 1)
	{
		m = (a + e) >> 1;
		if ( ARR[m] < V) a = m;	else e = m;
	};

	if (ARR[a]> V) return INVALID;
	if (ARR[a] == V) return a;	
	return e;
}

bool GetCombination(set &setBase, set &setCmb, UNSIGNED_2B_TYPE k)
//alternative interface for the main version below
{
	apvector<SIGNED_4B_TYPE> b, c;
	setBase.GetList(b);
	setCmb.GetList(c);
	
	bool res = GetCombination(b, c, k);
	if (res) setCmb = c; else setCmb.Dump();
	return res;
}

bool GetCombination(apvector<SIGNED_4B_TYPE> &Base, apvector<SIGNED_4B_TYPE> &Cmb, UNSIGNED_2B_TYPE k)
/*description:	returns in Cmb a combination of k elements from Base, 
				if Cmb is not empty and has k-size then it is used
				as a current combination to get the next one
 postcondition:
				returns true if new combination was created in Cmb
 precondition:				
				first call should have 0 < k < set-size.
				consecutive calls should have k = 0 and Cmb of k-size;
*/
{
	SIGNED_4B_TYPE N = Base.length();
	if (k)
	{//first call, populate combinations
		if (k == N) {Cmb = Base; return true; };
		Cmb.resize(0);
		if (k > N) return false;
		Cmb.resize(k);
		for (UNSIGNED_2B_TYPE i = 0; i < k; i++) Cmb[i] = Base[i];
		return true;
	}

	SIGNED_4B_TYPE C = Cmb.length();
	if ((C > N) || (C == 0))
	{
		Cmb.resize(0);
		return false;
	}

	SIGNED_4B_TYPE B, C1 = 0;
	while (C1 < C)
	{
		C1++;
		B = FindArrPoz(Base, Cmb[C - C1]);
		if (B < 0) break; //should not happen!
		if ( (B + C1) < N ) 
		{
			while (C1)	Cmb[C - C1--] = Base[++B];
			return true;
		}
	}

	Cmb.resize(0);
	return false; //all possibilities are exhausted
}

char CheckResidue(STRING_TYPE residue)
//converts 3-letter aminoacid codes into 1-letter ones
//i.e. ARNDCQEGHILKMFPSTWYV
{
	STRING_TYPE r = residue;
	r.parse_string();
	r.touppercase();

	//unknown
	if (r.find("UNK") == 0) return char(0);

	//aminoacids
	if (r.find(ALA_) == 0) return ALA1;
	if (r.find(ARG_) == 0) return ARG1;
	if (r.find(ASP_) == 0) return ASP1;
	if (r.find(ASN_) == 0) return ASN1;
	if (r.find(CYS_) == 0) return CYS1;
	if (r.find(GLN_) == 0) return GLN1;
	if (r.find(GLU_) == 0) return GLU1;
	if (r.find(GLY_) == 0) return GLY1;
	if (r.find(HIS_) == 0) return HIS1;
	if (r.find(ILE_) == 0) return ILE1;
	if (r.find(LEU_) == 0) return LEU1;
	if (r.find(LYS_) == 0) return LYS1;
	if (r.find(MET_) == 0) return MET1;
	if (r.find(PHE_) == 0) return PHE1;
	if (r.find(PRO_) == 0) return PRO1;
	if (r.find(SER_) == 0) return SER1;
	if (r.find(THR_) == 0) return THR1;
	if (r.find(TRP_) == 0) return TRP1;
	if (r.find(TYR_) == 0) return TYR1;
	if (r.find(VAL_) == 0) return VAL1;

	///nucleotides, 
	//modified bases have '+' in front
	if (r.find("A") == 0) return N_A;
	if (r.find("C") == 0) return N_C;
	if (r.find("G") == 0) return N_G;
	if (r.find("T") == 0) return N_T;
	if (r.find("U") == 0) return N_U;

	return char(0);
}
