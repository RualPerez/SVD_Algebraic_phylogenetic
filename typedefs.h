#if ! defined (TYPEDEFS_H) 
#define TYPEDEFS_H

/* the headers and typedefs used throughout the code */

#include <map>
#include <vector>
#include <iostream>

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

using namespace std;
typedef map<string,string> SShash;
typedef map<string,int> SIhash;
typedef vector<int> split;
typedef pair<char, int> MyPairType;
	

extern int VERBOSITY;
extern unsigned int NUM_SPECIES;

#endif
