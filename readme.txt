Get_ad estimates global applicability domain, i.e., expected coverage in the descriptor space,
1) for a given dataset A by a given dataset B or 2) self-coverage for a given dataset A

Get_ad operates on descriptor matrices, supporting the same file formats as in datasplit:
.svm, .xa, .x (see input_demo* files)
NB: Descriptor matrices, supplied as input, should have the same dimensions and scaling.

Get_ad calculates distribution of distances within the reference dataset A, 
then based on that, calculates a distance threshold = mean + st.dev * Z (an input parameter, 1.0 by default),
this distance threshold is then compared with the distance 
from each point in the test dataset B to its nearest neighbor in reference dataset A. 
If the threshold is higher, the respective test datapoint is "within AD", and "out" otherwise.
NB: This comparison can also be done with z-scores, by comparing (testpoint distance - mean) / st.dev with Z

The output is a tabulated file with actual and z-scored distances reported for the test dataset.

Usage example:
	"get_ad input_demo1.xa -4PRED=input_demo2.xa -OUT=report"


To print detailed help, run get_ad without parameters.
---------------------------------------------
get_ad v1.20 - calculates global Applicability Domain (AD) of a given data set
(see A.Tropsha, A.Golbraikh, Curr.Pharm.Des.(2007), 13, pp3494-3504).
Usage:  get_ad mdlfile [flags] or get_ad mdlfile -4PRED=tstfile [flags]
Input: .xa, .x or .svm model files; Output: .gad file

Possible flags are:    '-OUT=' output file name
'-4PRED=' test file to calc AD against modeling set, 
if test not set, self-AD of the modeling set is reported.
'-2PART' - outputs in and out of AD parts as separate files

'-K=' - number of nearest neighbors to use <def. 1>
'-Z=' - z-score cut-off to use <def. 1.0>
'-M=dx' = Metric to use; x - is a power coff. <def. 2.0>
d: E -Euclidean <def.>, T -Tanimoto, R -Correl, C -Cosine

'-D2AD' - to use squared distances <def.>; '-D1AD' - direct distances
'-F=..' - applicability domain filtering modes: 
'M' - av.dist to k neighbors should be within AD (traditional)
'A' - all k neighbors should be within AD
'H' - half neighbors within AD, 'L' - at least 1 within AD <def.>

