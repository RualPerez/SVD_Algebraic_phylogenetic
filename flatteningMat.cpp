// It constructs the flattening matrix of a given split.


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


vector<int> power_K;

vector<int> powersK(int numb_sp, int nb_alph) {
	int aux =1;
	for (int i=0; i <= numb_sp; i++) {
		power_K.push_back(aux);
		aux = nb_alph*aux ;
	}

	return(power_K);
}

// Calculate the coordinates of the flattening matrix for a given column
unsigned int getCoord( string it, vector<int> set_species, 
	map<char,int> alph , vector<int> power_K ) {
	
	int nb_set = set_species.size();
	long int coord = 0;

	vector<long int> coords;	
	for (int i =0; i < nb_set; i ++) {
		char nucleotide = it[ set_species[i] ];
		coord += alph[nucleotide]*power_K[ nb_set-1-i ];
	}

	return coord;
}

mat flatteningMat( map<string,int> tensor, map<char,int> alph,
	vector<int> split, vector<int> othersp, unsigned int seq_len ) {

	unsigned int coord_row, coord_col;

	int n_rows = pow( alph.size(), split.size() );
	int n_cols = pow( alph.size(), othersp.size() );
	vector<int> power_K = powersK( max(split.size(), othersp.size()), alph.size() );
	
	mat flatmatrix = zeros<mat>(n_rows, n_cols);
	
	for ( map<string,int>::iterator it = tensor.begin(); it != tensor.end(); ++it ) {
		coord_row = getCoord( it->first, split, alph, power_K );
		coord_col = getCoord( it->first, othersp, alph, power_K );
		flatmatrix( coord_row, coord_col ) = it->second;
	}

	return flatmatrix;
}