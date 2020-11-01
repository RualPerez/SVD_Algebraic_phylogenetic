// It reads an alignment to obtain the so-called Alignment struct.

#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;	//using namespace std;
using namespace arma; // for 'armadillo';

/*struct Alignment {
	unsigned int num_taxa;
	unsigned int seq_len;
	vector<string> taxa; // the taxa names
	map<string, string> seqs; // dictionary : taxa name - sequence
};*/

Alignment readFASTA (string fname) {

	Alignment align;
	string line;
	fstream file;
	string speciesname;
	string seq = "";
	int numchar;
	vector<int> al_length;
	
	// Opening file
	file.open(fname.c_str(), fstream::in);
	if ( file == NULL ) {
		cout << " cannot open file " << fname.c_str() << "\n";
		exit(1);
	}
	
	// Get the first line from the FASTA file
	getline(file, line, '\n');
	if ( line[0]!= '>') {
		cout << " error : not a FASTA file \n";
		exit(0);
	}
	numchar = line.length()-1;

	speciesname = line.substr(1, numchar);
	align.taxa.push_back( speciesname );
	while (! file.eof()) {
		getline(file, line, '\n');
		// Storage all the lines until it finds the next name specie
		if ( line[0] != '>') seq += line;
		else {
			align.seqs[speciesname] = seq;
			al_length.push_back( seq.length() );
			numchar = line.length()-1;
			speciesname = line.substr(1, numchar);
			align.taxa.push_back( speciesname);
			seq.clear();
			if ( align.seqs[speciesname].length() > 0)
				cerr << " Warning - found 2+ sequences with same name " << speciesname << endl;
		}
	}
	
	align.seqs[speciesname] = seq;
	al_length.push_back(seq.length());
	file.close();
	
	align.num_taxa = align.taxa.size();
	if( al_length.size() != align.num_taxa ) cout << " error in getting the alignment ";
	for ( int i=0; i < al_length.size()-1; i++) {
		if ( al_length[i] != al_length[i+1]) {
			cout << " sequence " << align.taxa[i+1] << " does not have same length \n";
			exit(0);
		}
	}
	align.seq_len = al_length.back();
	return align;
}