// It writes the taxa names as 1, 2, 3, 4, following the order of appearance in the FASTA file

#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;	//using namespace std;

vector <int> realOrder( vector<string> taxa ) {

	int nb_taxa = taxa.size();
	vector<int> order_taxa( nb_taxa );
	
	map<string,int> sequences;
	map<string,int>:: iterator it;
	
	for (int i=0; i < nb_taxa; ++i ) sequences[ taxa[i] ] = 0;
	
	int alph_order = 0;
	for (it = sequences.begin() ; it != sequences.end() ; ++it) {
		sequences[ it->first ] = alph_order;
		++alph_order;
	}
	
	for (int i=0; i < nb_taxa; ++i ) {
		it = sequences.find( taxa[i] );
		order_taxa[i] = it->second;
	}
	
	return order_taxa;
}