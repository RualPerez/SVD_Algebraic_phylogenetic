// It computes the Frobenius-distance between the given flattening matrix and the space
// of matrices with rank ≤ nb_partititons*nb_alph

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

double distFrobenius ( mat flat, int nb_partitions, int nb_alph ) {

	double tail = 0;
	vec singval = svd(flat); // SVD function from Armadillo library

	// Frobenius distance
	for (int k = nb_alph*nb_partitions; k < singval.size(); k++) 
		tail += singval[k]*singval[k];
	
	return sqrt(tail);
}