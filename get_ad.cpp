/* get_ad.cpp Sept 19 2011 - present
Console application (Linux, Windows). 

Based on dataset-class
TODO:

DONE: 

(history list of recent changes)
1.0			Sep 20 2011		Basic version completed
			
Tested with:
MRP5x_moe_n_mdl1 -K=5 -F=A -2PART -Z=0.5
MRP5x_moe_n_mdl1 -4PRED=MRP5x_moe_n_ext1.xa -K=3 -F=H -2PART -Z=0.0
MRP5x_moe_n_mdl1 -4PRED=MRP5x_moe_n_ext1.xa -K=3 -F=H -2PART
*/

#include "dataset.h"

#define Version "1.0"
#define COMMENT		"#"
#define	AD_FILE		".ad"

//Global variables ------------------------------------
//distance metric
REALNUM_TYPE METRIC_V = 2.0;
UNSIGNED_1B_TYPE METRIC_K = 0;

//AD settings, etc.
REALNUM_TYPE AD_dist = INVALID, AD_Z = 1.0;
UNSIGNED_1B_TYPE AD_K = 1, AD_MODE = 1;
bool D1dist = false, ExplicitReport = false;

STRING_TYPE stInput, stOutput, stScreen;
//End of global variables ------------------------------

void SaveDataset(dataset &dtstD, STRING_TYPE stName, SIGNED_2B_TYPE mtxTYPE)
{
	STRING_TYPE stO = stName;
	if ( !(CheckStrEnding(stO, ".x")||CheckStrEnding(stO, ".xa")||CheckStrEnding(stO, ".svm")) )
	switch (mtxTYPE)
	{
		case 1:	stO += ".x";		break;
		case 2:	stO += ".svm";		break;
		case 0:	default:
				stO += ".xa";		break;
	}

	FILETYPE_OUT foD(stO.c_str());
	if (mtxTYPE == 1)
	{
		CutStrEnding(stO);
		stO += ".a";
		FILETYPE_OUT foAfile(stO.c_str());
		dtstD.save(foD, &foAfile, 1);
		foAfile.close();
	}
	else
		dtstD.save(foD, NULL, (UNSIGNED_1B_TYPE)mtxTYPE);
	foD.close();
}

SIGNED_2B_TYPE LoadDataset(dataset &dtstD, STRING_TYPE &stName)
{
	//input matrix format, new one by default
	SIGNED_2B_TYPE nTYPE = 0;
	STRING_TYPE stMTX = stName;
	stMTX.tolowercase();

	if (!CheckStrEnding(stMTX, ".xa") && !CheckStrEnding(stMTX, ".x") && !CheckStrEnding(stMTX, ".svm"))
	{//if no extension then by default assume a new format
		CutStrEnding(stName);
		stName += ".xa";
		FILETYPE_IN fiTest(stName.c_str());
		if (fiTest.eof() || fiTest.fail())
		{//revert to old format
			CutStrEnding(stName);
			stName	 += ".x";
			stMTX	 += ".x";
		}
		else
			stMTX	 += ".xa";
		fiTest.close();
	}
	
	FILETYPE_IN fiXAfile(stName.c_str());
	if (fiXAfile.eof() || fiXAfile.fail())
	{
		cout << "Can not open input file: '" << stName << "'" << endl;
		return -1;
	}

	bool ifRead = true;
	if (CheckStrEnding(stMTX, ".x"))
	{//old format
		nTYPE = 1;
		STRING_TYPE stA = stName;
		CutStrEnding(stA);
		stA += ".a";
		FILETYPE_IN fiAfile(stA.c_str());
		if (fiAfile.eof() || fiAfile.fail())
			ifRead = dtstD.load(fiXAfile, NULL, 10);
		else
			ifRead = dtstD.load(fiXAfile, &fiAfile, 1);
		fiAfile.close();
	}
	else
	{
		if (CheckStrEnding(stMTX, ".svm")) nTYPE = 2;
		ifRead = dtstD.load(fiXAfile, NULL, (UNSIGNED_1B_TYPE)nTYPE);
	}
	fiXAfile.close();
	if (ifRead) return (nTYPE);
	return -1;
}

void PrintHelp()
{//print help
	cout << endl << "get_ad v" << Version << " - calculates global Applicability Domain (AD) of a given data set" << endl;
	cout << "(see A.Tropsha, A.Golbraikh, Curr.Pharm.Des.(2007), 13, pp3494-3504)." << endl;
	cout << "Usage:  get_ad mdlfile [flags] or get_ad mdlfile -4PRED=tstfile [flags]" << endl;
	cout << "Input: .xa, .x or .svm model files; Output: " << AD_FILE << " file" << endl << endl;

	cout << "Possible flags are:    '-OUT=' output file name" << endl;	
	cout << "'-4PRED=' test file to calc AD against modeling set, " << endl;
	cout << "if test not set, self-AD of the modeling set is reported." << endl;
	cout << "'-2PART' - outputs in and out of AD parts as separate files" << endl << endl;

	cout << "'-K=' - number of nearest neighbors to use <def. 1>" << endl;
	cout << "'-Z=' - z-score cut-off to use <def. 1.0>" << endl;		
	cout << "'-M=dx' = Metric to use; x - is a power coff. <def. 2.0>" << endl;
	cout << "d: E -Euclidean <def.>, T -Tanimoto, R -Correl, C -Cosine" << endl << endl;
	
	cout << "'-D2AD' - to use squared distances <def.>; '-D1AD' - direct distances" << endl;
	cout << "'-F=..' - applicability domain filtering modes: " << endl;				
	cout << "'M' - av.dist to k neighbors should be within AD (traditional)" << endl;
	cout << "'A' - all k neighbors should be within AD" << endl;
	cout << "'H' - half neighbors within AD, 'L' - at least 1 within AD <def.>" << endl << endl;
}

void ProcessArgumentString(STRING_TYPE &S)
{
	REALNUM_TYPE rtX;
	
	S.parse_string();
	STRING_TYPE stX = S;
	S.touppercase();

	SIGNED_4B_TYPE intX, intU = S.find("-OUT="); 
	if (intU == 0)
	{//output file
		stOutput = stX.substr(intU + 5, stX.length());
		return;
	}

	intU = S.find("-4PRED=");
	if (intU == 0)
	{//output file
		stScreen = stX.substr(intU + 7, stX.length());
		stScreen.parse_string();
		return;
	}

	if (S.find("-2PART") == 0)
	{//to generate explicit descriptor matrices for in-AD and out-of-AD potions of the input file
		ExplicitReport = true;
		return;
	}

	if (S.find("-D1AD") == 0)
	{//to use direct distances for AD evaluation
		D1dist = true;
		return;
	}

	intU = S.find("-Z="); //Z-cut-off
	if (intU == 0)
	{
		stX = S.substr(intU + 3, S.length());		
		rtX = atof( stX.c_str() );
		AD_Z  = rtX;
		return;
	}

	intU = S.find("-M="); //distance metric to use
	if (intU == 0)
	{
		if (S[intU + 3] == 'C') METRIC_K = 1;
		if (S[intU + 3] == 'R') METRIC_K = 2;
		if (S[intU + 3] == 'T') METRIC_K = 3;
		if (S.length() > intU + 4)
		{
			rtX = atof( S.substr(intU + 4, S.length()).c_str() );
			if (rtX > 0)	METRIC_V = rtX;
		}
		return;
	}

	intU = S.find("-K="); //sets k or R; e.g. -KR=1@9 -KR=0.5
	if (intU == 0)
	{
		stX = S.substr(intU + 3, S.length());		
		intX = atoi( stX.c_str() );
		if (intX > 0) AD_K = intX;		
		return;
	}

	intU = S.find("-F="); //applicability domain settings, e.g. -AD=0.5 or -AD=0.5d1 or -AD=0.5_avk
	if (intU == 0)
	{
		stX = S.substr(intU + 3, S.length());		
		if (S[intU + 3] == 'M') AD_MODE = 0;
		if (S[intU + 3] == 'A') AD_MODE = 3;
		if (S[intU + 3] == 'H') AD_MODE = 2;
		if (S[intU + 3] == 'L') AD_MODE = 1; //default
		return;
	}
}

bool VerifyDescriptors(dataset &A, dataset &B)
//checks if descriptors are the same in A and B datasets:
{
	SIGNED_4B_TYPE dn = A.get_Ndimensions();
	if (B.get_Ndimensions() != dn) return false;

	for (SIGNED_4B_TYPE i = 0; i < dn; i++)
		if (i == B.get_dscr_pos(A.get_dscr(i))) 
			continue;
		else
			return false;
	
	return true;
}

int main(int argc, char* argv[])
{	
	if (argc < 2) 
	{
		PrintHelp();	
		return 0;
	}
	
	stInput = argv[1];
	
	//temporary service variables
	STRING_TYPE stTime, stArg, stJ = "command line: get_ad ";
	GetTimeStamp(stTime);

	dataset datasetX, datasetC;	
	SIGNED_2B_TYPE inMTX = LoadDataset(datasetX, stInput); //save input matrix format
		
	if (inMTX == -1)
	{
		cout << "Can not load data from input file." << endl;
		return -1;
	}
	
	SIGNED_4B_TYPE i = 1;	
	stJ += argv[i];
	//---------------------------------------------------------
	while (argc > ++i)
	{ 
		stArg = argv[i];
		stJ += BLANK;
		stJ += stArg;
		ProcessArgumentString(stArg);
	}//while (argc > ++nArg) loop	
	//---------------------------------------------------------

	dataset * pD = NULL;
	bool SelfAD = false;
	if (stScreen.length() == 0)
	{
		SelfAD = true;
		stScreen = stInput;
		pD = &datasetX;
	}
	else
	{
		if (LoadDataset(datasetC, stScreen) == -1)
		{
			cout << "Cannot load data from screening file: " << stScreen << endl;
			return -1;
		}
		
		if ( !VerifyDescriptors(datasetX, datasetC) )
		{
			cout << "Descriptors are mismatched in the input and screening files." << endl;
			return -1;
		}
		pD = &datasetC;
	}

	//log time and settings:
	if (stOutput.length() == 0)
	{
		if (SelfAD) 
			stOutput = "SelfAD"; 
		else 
		{
			stOutput = stScreen;
			CutStrEnding(stOutput);			
		}
		stOutput += "_by_" + stInput;
		CutStrEnding(stOutput);
	}

	if (!CheckStrEnding(stOutput, AD_FILE))
		stOutput += AD_FILE;

	FILETYPE_OUT foAD(stOutput.c_str());
	foAD << COMMENT << stJ << endl;
	foAD << COMMENT << "Started at " << stTime << endl;	

	//calculate Global AD cut-off
	datasetX.calc_dist(0, METRIC_V, METRIC_K);

	QSAR qsarT;
	apvector<REALNUM_TYPE> knn_stats, dsx, kneib;
	if (D1dist)
		datasetX.get_NearNeibDistances(knn_stats, AD_K, 0, 2);
	else
		datasetX.get_NearNeibDistances(knn_stats, AD_K, 0, 3);
	
	REALNUM_TYPE mnd, f = knn_stats[0] + AD_Z*knn_stats[2];
	if (D1dist) AD_dist = f; else AD_dist = sqrt(f);
	foAD << COMMENT << "AD cut-off=" << AD_dist << endl;
	foAD << "ID\tSID\tDist\tZ-score\tWITHIN_AD" << endl;
	char buff[200] = "";	//to format numerical output

	//calculate distances for the screened set
	SIGNED_4B_TYPE j, k, l, n = pD->get_Ndatapoints(), bn = datasetX.get_Ndatapoints();
	kneib.resize(AD_K);
	pD->test.Dump();
	pD->train.Dump();

	for (i = 0; i < n; i++)
	{
		if (SelfAD)
		{
			apvector<SIGNED_4B_TYPE> neibs;
			datasetX.get_NearNeib(i, neibs, AD_K, 0);
			for (k = 0; k < AD_K; k++)	kneib[k] = datasetX.get_Distance(i, neibs[k]);
		}
		else
		{
			for (k = 0; k < AD_K; k++)	kneib[k] = INVALID;	//invalidate k-neighbor distances

			pD->get_DimValues(i, dsx);
			datasetX.add_dp(dsx);			

			for (j = 0; j < bn; j++)
			{//store k nearest distances
				f = datasetX.get_indDistance(bn, j, NULL, METRIC_V, METRIC_K);
				for (k = 0; k < AD_K; k++)
				if ( (kneib[k] > f) || (kneib[k] == INVALID) )
				{
					if (kneib[k] > f)
					{
						for (l = k+1; l < AD_K; l++)
							kneib[l] = kneib[l-1];
					}
					kneib[k] = f; 
					break;	
				}// for k
			}// for j

			datasetX.remove_dp(bn);
		}//else...if SelfAD

		switch (AD_MODE)
		{//calculate AD distance for the tested point
			case 0: //average dist of k neighbors should be within AD
				mnd = qsarT.meanV(kneib);
			break;
			
			case 2: //half of the neighbors should be within AD
				k = (AD_K >> 1);
				if ( AD_K == (k << 1) )
					mnd = kneib[k - 1];
				else
					mnd = kneib[k];
			break;
			
			case 3: //all of the neighbors should be within AD
				mnd = kneib[AD_K - 1];
			break;

			case 1: //the closest neighbor should be within AD
			default:
				mnd = kneib[0];
			break;
		}//switch (AD_MODE)

		//calculate actual z-score (can be compared with the AD_Z cutoff to define "in" or "out" of AD)
		if (D1dist) f = mnd; else f = mnd*mnd;		
		f -= knn_stats[2];
		f /= knn_stats[0];
		
		stJ = pD->get_sid(i);
		sprintf(buff, "%d\t%s\t%6.3f\t%6.3f\t%1d", i, stJ.c_str(), mnd, f, UNSIGNED_1B_TYPE(mnd < AD_dist));
		stJ = buff;
		stJ.parse_string();
		foAD << stJ << endl;
		
		if (ExplicitReport)
		{
			if (mnd < AD_dist) 
				pD->train.PutInSet(i); 
			else 
				pD->test.PutInSet(i);
		}
	}//for i
	
	
	if (ExplicitReport)
	{//save "within-AD" and "out-of-AD" portions into separate files
		CutStrEnding(stScreen);
		if (!pD->test.IsEmpty())
		{
			dataset X = pD->get_test_set();
			SaveDataset(X, stScreen + "_out_ad", inMTX);
		}

		if (!pD->train.IsEmpty())
		{
			dataset X = pD->get_training_set();
			SaveDataset(X, stScreen + "_in_ad", inMTX);
		}
	}

	foAD.close();
	cout << "Done." << endl;
	return 0;
}

