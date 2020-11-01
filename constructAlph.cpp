// It takes all the different characters of the sequence in order to determine the alphabet
// of the alignment. Note that the alignment can induce an alphabet with any size.


#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;	//using namespace std;

bool find_ChrExclude(string chr_exclude, char to_find) {
	int len = chr_exclude.size();
	for(int i=0; i<len; ++i) {
		if(chr_exclude[i]==to_find) return true;
	}
	return false;	
}

map<char,int> constructAlph ( Alignment &align, string chr_exclude ) {
	int value = 0;
	int al_length = align.seq_len;
	map<char,int> alph;
	map<string,string> myseqs = align.seqs;
	map<string,string>::const_iterator pos;
	
	// Go through the entire alignment in search of
	// different characters
	for ( int i =0; i < al_length; i ++) {
		for ( pos = myseqs.begin(); pos != myseqs.end(); ++ pos ) {
			bool exclude = find_ChrExclude( chr_exclude, pos->second[i] );
			if (not exclude and alph.find( pos->second[i] ) == alph.end() ) {
				alph[ pos->second[i] ] = value;
				++value;
			}
		}
	}
	return alph;
}