//bonds, elements, aminoacids, and nucleotides definitions.

//bonds
#define		NO_BOND			0
#define		SINGLE_BOND		1
#define		DOUBLE_BOND		2
#define		TRIPLE_BOND		3
#define		AROMATIC_BOND	4
#define		CONJ_BOND		5
#define		TAUTAMERIC_BOND 8
#define		HYDROGEN_BOND	10
#define		A_D_BOND		20			//donor-acceptor
#define		IONIC_BOND		30
#define		UNKNOWN_BOND    40
#define		DELETE_BOND     100			//only if existing connection is there



//some atomic elements
#define		HYDROGEN		1
#define		HELIUM			2
#define		LITHIUM			3
#define		BERYLLIUM		4
#define		BORON			5
#define		CARBON			6
#define		NITROGEN		7
#define		OXYGEN			8
#define		FLUORINE		9
#define		NEON			10
#define		NATRIUM			11
#define		MAGNESIUM		12
#define		ALUMINIUM		13
#define		SILICON			14
#define		PHOSPHORUS		15
#define		SULFUR			16
#define		CHLORINE		17
#define		ARGON			18
#define		POTASSIUM		19
#define		CALCIUM			20
#define		SCANDIUM		21
#define		TITANIUM		22
#define		VANADIUM		23
#define		CHROMIUM		24
#define		MANGANESE		25
#define		IRON			26
#define		COBALT			27
#define		NICKEL			28
#define		COPPER			29
#define		ZINC			30
#define		GALLIUM			31
#define		GERMANIUM		32
#define		ARSENIC			33
#define		SELENIUM		34
#define		BROMINE			35
#define		KRYPTON			36
#define		RUBIDIUM		37
#define		STRONTIUM		38
#define		YTTRIUM			39
#define		ZIRCONIUM		40
#define		NIOBIUM			41
#define		MOLYBDENUM		42
#define		TECHNETIUM		43
#define		RUTHENIUM		44
#define		RHODIUM			45
#define		PALLADIUM		46
#define		SILVER			47
#define		CADMIUM			48
#define		INDIUM			49
#define		TIN				50
#define		ANTIMONY		51
#define		TELLURIUM		52
#define		IODINE			53
#define		XENON			54
#define		CESIUM			55
#define		BARIUM			56
#define		LANTHANUM		57
#define		CERIUM			58
#define		PRASEODYMIUM	59
#define		NEODYMIUM		60
#define		PROMETHIUM		61
#define		SAMARIUM		62
#define		EUROPIUM		63
#define		GADOLINIUM		64
#define		TERBIUM			65
#define		DYSPROSIUM		66
#define		HOLMIUM			67
#define		ERBIUM			68
#define		THULIUM			69
#define		YTTERBIUM		70
#define		LUTETIUM		71
#define		HAFNIUM			72
#define		TANTALUM		73
#define		TUNGSTEN		74
#define		RHENIUM			75
#define		OSMIUM			76
#define		IRIDIUM			77
#define		PLATINUM		78
#define		GOLD			79
#define		MERCURY			80
#define		THALLIUM		81
#define		LEAD			82
#define		BISMUTH			83
#define		POLONIUM		84
#define		ASTATINE		85
#define		RADON			86
#define		FRANCIUM		87
#define		RADIUM			88
#define		ACTINIUM		89
#define		THORIUM			90
#define		PROTACTINIUM	91
#define		URANIUM			92
#define		NEPTUNIUM		93
#define		PLUTONIUM		94
#define		AMERICIUM		95
#define		CURIUM			96
#define		BERKELIUM		97
#define		CALIFORNIUM		98
#define		EINSTEINIUM		99
#define		FERMIUM			100
#define		MENDELEVIUM		101
#define		NOBELIUM		102
#define		LAWRENCIUM		103


#define	S_HYDROGEN			"H"
#define	S_HELIUM			"He"
#define	S_LITHIUM			"Li"
#define	S_BERYLLIUM			"Be"
#define	S_BORON				"B"
#define	S_CARBON			"C"
#define	S_NITROGEN			"N"
#define	S_OXYGEN			"O"
#define	S_FLUORINE			"F"
#define	S_NEON				"Ne"
#define	S_NATRIUM			"Na"
#define	S_MAGNESIUM			"Mg"
#define	S_ALUMINIUM			"Al"
#define	S_SILICON			"Si"
#define	S_PHOSPHORUS		"P"
#define	S_SULFUR			"S"
#define	S_CHLORINE			"Cl"
#define	S_ARGON				"Ar"
#define	S_POTASSIUM			"K"
#define	S_CALCIUM			"Ca"
#define	S_SCANDIUM			"Sc"
#define	S_TITANIUM			"Ti"
#define	S_VANADIUM			"V"
#define	S_CHROMIUM			"Cr"
#define	S_MANGANESE			"Mn"
#define	S_IRON				"Fe"
#define	S_COBALT			"Co"
#define	S_NICKEL			"Ni"
#define	S_COPPER			"Cu"
#define	S_ZINC				"Zn"
#define	S_GALLIUM			"Ga"
#define	S_GERMANIUM			"Ge"
#define	S_ARSENIC			"As"
#define	S_SELENIUM			"Se"
#define	S_BROMINE			"Br"
#define	S_KRYPTON			"Kr"
#define	S_RUBIDIUM			"Rb"
#define	S_STRONTIUM			"Sr"
#define	S_YTTRIUM			"Y"
#define	S_ZIRCONIUM			"Zr"
#define	S_NIOBIUM			"Nb"
#define	S_MOLYBDENUM		"Mo"
#define	S_TECHNETIUM		"Tc"
#define	S_RUTHENIUM			"Ru"
#define	S_RHODIUM			"Rh"
#define	S_PALLADIUM			"Pd"
#define	S_SILVER			"Ag"
#define	S_CADMIUM			"Cd"
#define	S_INDIUM			"In"
#define	S_TIN				"Sn"
#define	S_ANTIMONY			"Sb"
#define	S_TELLURIUM			"Te"
#define	S_IODINE			"I"
#define	S_XENON				"Xe"
#define	S_CESIUM			"Cs"
#define	S_BARIUM			"Ba"
#define	S_LANTHANUM			"La"
#define	S_CERIUM			"Ce"
#define	S_PRASEODYMIUM		"Pr"
#define	S_NEODYMIUM			"Nd"
#define	S_PROMETHIUM		"Pm"
#define	S_SAMARIUM			"Sm"
#define	S_EUROPIUM			"Eu"
#define	S_GADOLINIUM		"Gd"
#define	S_TERBIUM			"Tb"
#define	S_DYSPROSIUM		"Dy"
#define	S_HOLMIUM			"Ho"
#define	S_ERBIUM			"Er"
#define	S_THULIUM			"Tm"
#define	S_YTTERBIUM			"Yb"
#define	S_LUTETIUM			"Lu"
#define	S_HAFNIUM			"Hf"
#define	S_TANTALUM			"Ta"
#define	S_TUNGSTEN			"W"
#define	S_RHENIUM			"Re"
#define	S_OSMIUM			"Os"
#define	S_IRIDIUM			"Ir"
#define	S_PLATINUM			"Pt"
#define	S_GOLD				"Au"
#define	S_MERCURY			"Hg"
#define	S_THALLIUM			"Tl"
#define	S_LEAD				"Pb"
#define	S_BISMUTH			"Bi"
#define	S_POLONIUM			"Po"
#define	S_ASTATINE			"At"
#define	S_RADON				"Rn"
#define	S_FRANCIUM			"Fr"
#define	S_RADIUM			"Ra"
#define	S_ACTINIUM			"Ac"
#define	S_THORIUM			"Th"
#define	S_PROTACTINIUM		"Pa"
#define	S_URANIUM			"U"
#define	S_NEPTUNIUM			"Np"
#define	S_PLUTONIUM			"Pu"
#define	S_AMERICIUM			"Am"
#define	S_CURIUM			"Cm"
#define	S_BERKELIUM			"Bk"
#define	S_CALIFORNIUM		"Cf"
#define	S_EINSTEINIUM		"Es"
#define	S_FERMIUM			"Fm"
#define	S_MENDELEVIUM		"Md"
#define	S_NOBELIUM			"No"
#define	S_LAWRENCIUM		"Lr"



//aminoacids
#define		ALA1		'A'
#define		ARG1		'R'
#define		ASN1		'N'
#define		ASP1		'D'
#define		CYS1		'C'
#define		GLN1		'Q'
#define		GLU1		'E'
#define		GLY1		'G'
#define		HIS1		'H'
#define		ILE1		'I'
#define		LEU1		'L'
#define		LYS1		'K'
#define		MET1		'M'
#define		PHE1		'F'
#define		PRO1		'P'
#define		SER1		'S'
#define		THR1		'T'
#define		TRP1		'W'
#define		TYR1		'Y'
#define		VAL1		'V'

#define		ALA_		"ALA"
#define		ARG_		"ARG"
#define		ASN_		"ASN"
#define		ASP_		"ASP"
#define		CYS_		"CYS"
#define		GLN_		"GLN"
#define		GLU_		"GLU"
#define		GLY_		"GLY"
#define		HIS_		"HIS"
#define		ILE_		"ILE"
#define		LEU_		"LEU"
#define		LYS_		"LYS"
#define		MET_		"MET"
#define		PHE_		"PHE"
#define		PRO_		"PRO"
#define		SER_		"SER"
#define		THR_		"THR"
#define		TRP_		"TRP"
#define		TYR_		"TYR"
#define		VAL_		"VAL"

//ARNDCQEGHILKMFPSTWYV

///nucleotides
#define		N_A				char(21);	//'\a'
#define		N_C				char(22);	//'\c'
#define		N_G				char(23);	//'\g'
#define		N_T				char(24);	//'\t'
#define		N_U				char(25);	//'\u'

#define		DEFAULT_LARGE_CYCLE	9				//cycles of this length or greater are not considered to be affecting chemical properties of an atom

//constants for stereoisomers
#define		CIS					2
#define		TRANS				4

#define		Z_ISOMER			20
#define		E_ISOMER			40

#define		CLOCKWISE			2
#define		COUNTERCLOCKWISE	1

#define		R_ISOMER			6
#define		S_ISOMER			7
