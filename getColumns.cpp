// From an alignment, it obtains the absolute frequencies of each column; it performs
 //the joint distribution vector (which may not contain all possible patterns).

#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;	//using namespace std;


bool findChrExclude ( string chr_exclude, char to_find ) {
	/* Check that the considered character is not equal to one
	   of those that the user wanted to exclude ( usually gaps ) */

	int len = chr_exclude.size();
	for (int i =0; i < len; ++ i ) {
		if( chr_exclude[i]== to_find ) return true;
	}
	return false;
}

map<string, int> getColumns ( Alignment &align, string chr_exclude ) {
	int al_length = align.seq_len;
	map<string,string> myseqs = align.seqs;
	map<string,string>::const_iterator pos;
	map<string,int> columnscount;
	
	// Loop for all columns of the alignment
	for ( int i =0; i < al_length; i++ ) {
		bool exclude = false;
		string col = "";
		
		// Loop for all the taxa sequences to capture each column
		for ( pos = myseqs.begin(); pos != myseqs.end() and not exclude; ++ pos ) {
				col.push_back(pos->second[i]);
				exclude = findChrExclude(chr_exclude, pos->second[i]);
			}
			if( not exclude ) columnscount[col]++;
		}
	return (columnscount);
}