//Basic collection of type and macro definitions.

#if !defined(CORE_TYPES)
#define CORE_TYPES

//C headers
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

//---------------------------------------------------------
//use WIN_PORT2 to enable Windows-specific parts of code
//#define WIN_PORT2

#ifdef WIN_PORT2
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#include <windows.h>
//#include <windowsx.h>
//#include <mmsystem.h>
#include <commctrl.h>
#include <commdlg.h>
#include <iomanip>
#include <malloc.h>
#endif
//---------------------------------------------------------



#ifdef _DEBUG
//#define MEMLEAK_MODE		//uncomment for memory leak catching
//#define ADV_LEAK_CATCHER	//uncomment for the advanced memory debugging stuff
#endif 

#ifdef MEMLEAK_MODE
#include "leax.h"		//memory leak catch + declaration of some macros
#endif


//#define _NO_NAMESPACE
#define _NO_EXCEPTION

#include "matrix.h"

#include "apvector.h"
#include "apstring.h"

#define		ZERO				0
#define		SET					1
#define		INVALID				-1
#define		LABEL				10000
#define     BIG_NUMBER			1000
#define		SMALL_NUMBER		0.000000001	//0.00001		//small number	- very important, affects precision of geometry resolution by get3d and molecular mechanics
#define		SMALL				0.1							//minimum interatomic distance

#define		DEF_STRING_SIZE		256
#define		DIGITtoCHAR			48

#define		ALLBITS1BYTE		0xFF

#define	DEF_FREAD_MOD			ios_base::in				//ios::nocreate | ios::in
#define sqr(X)		((X)*(X))
#define pi						3.141592653589793
#define D2RAD					0.01745329					//(pi/180)//for conversion degrees to radians

#define	MODE_WORK				0
#define MODE_ABORT				1
#define MODE_PAUSE				2

//swapping elements MACRO
#define Swap(a, b, t)		{(t) = (a); (a) = (b); (b) = (t);}

//output formatting spacers
#define BLANK	" "
#define TAB		"\t"

#ifdef WIN_PORT2
//a macro to use inside subroutines
#define CLEAR_WINMESSAGES	MSG msg;	while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)){\
	TranslateMessage(&msg);		DispatchMessage(&msg); }
#else
#define CLEAR_WINMESSAGES		//do nothing
#endif

//types
typedef unsigned char	UNSIGNED_1B_TYPE;		//1 byte
typedef signed char		SIGNED_1B_TYPE;			//1 byte
typedef signed short	SIGNED_2B_TYPE;			//2 bytes
typedef unsigned short 	UNSIGNED_2B_TYPE;		//2 bytes
typedef unsigned int	UNSIGNED_4B_TYPE;		//4 bytes
typedef signed int		SIGNED_4B_TYPE;			//4 bytes
//typedef char *			STRING_TYPE;			//can be Cstring but then MFC support will be needed
typedef apstring		STRING_TYPE;			//can be Cstring but then MFC support will be needed
typedef double			REALNUM_TYPE;			//8 bytes
//typedef float			SREALNUM_TYPE;			//4 bytes
typedef void *			GENERIC_POINTER;
typedef ofstream		FILETYPE_OUT;
typedef ifstream		FILETYPE_IN;


#define GRAB_MEM_BLOCK(x)		new x					//(x*)(malloc(sizeof(x)))
#define DROP_MEM_BLOCK(p)		delete p				//free(p)

#define GRAB_MEM_BLOCKS(x, n)	new x[n]				//(x*)(malloc(n*sizeof(x)))
#define DROP_MEM_BLOCKS(p)		delete [] p				//free(p)


#endif		//#define CORE_TYPES
