/* ###################################
#				 TFG 2018 			 #
# 			  Code : SVD.cpp 		 #
#################################### */

/*
Copyright (C) 2018 Raül Pérez, Casanellas and Fernandez-Sanchez.
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <unistd.h>

using namespace std;
using namespace arma; // for ’armadillo ’

int main(int argc, char *argv[]) {

	// Fasta filename
	string filename = argv[1];
	// Characters that do not belong to the alphabet, i.e., gaps
	string chr_exclude = argv[2];
	// each partition has evolved according to a different evolutionary model (valid input: 1, 2 or 3)
	int nb_partitions = stoi(argv[3]);
	
	double dist_F1, dist_F2, dist_F3;
	
	Alignment align = readFASTA(filename);
	// Vector of joint distribution
	map<string,int> tensor = getColumns(align, chr_exclude);

	// Construct the alphabet from the data
	map<char,int>alph = constructAlph(align, chr_exclude);
	int nb_alph = alph.size();
	
	// Sort the name of taxa by the order of appearance of the fasta file
	vector<int> order_taxa = realOrder(align.taxa);
	vector<int> split_1, split_2;
	
	// Flattening matrix for 12|34 and its Frobeniusdistance to E_{rk}
	split_1 = realSplit({0,1}, order_taxa);
	split_2 = realSplit({2,3}, order_taxa);
	mat flat_01 = flatteningMat(tensor, alph, split_1, split_2, align.seq_len);
	dist_F1 = distFrobenius(flat_01, nb_partitions, nb_alph);
	
	// Flattening matrix 13|24 and its Frobenius distance to E_{rk}
	split_1 = realSplit({0,2}, order_taxa);
	split_2 = realSplit({1,3}, order_taxa);
	mat flat_02 = flatteningMat(tensor, alph, split_1, split_2, align.seq_len);
	dist_F2 = distFrobenius(flat_02,nb_partitions, nb_alph);
	
	// Flattening matrix 14|23 and its Frobenius distance to E_{rk}
	split_1 = realSplit({0,3}, order_taxa);
	split_2 = realSplit({1,2}, order_taxa);
	mat flat_03 = flatteningMat(tensor, alph, split_1, split_2, align.seq_len);
	dist_F3 = distFrobenius(flat_03, nb_partitions, nb_alph);
	
	double min_score = min(dist_F1, dist_F2);
	min_score = min(min_score, dist_F3);
	
	cout << endl;
	cout << " Taxa 1: " << align.taxa[0] << endl;
	cout << " Taxa 2: " << align.taxa[1] << endl;
	cout << " Taxa 3: " << align.taxa[2] << endl;
	cout << " Taxa 4: " << align.taxa[3] << endl;
	cout << endl;
	cout << " Topologies : \t" << " 12|34 \t \t" << " 13|24 \t \t" << " 14|23 " << endl;
	cout << " Frobenius distance : \t" << dist_F1 << "\t \t" << dist_F2 << "\t \t" << dist_F3 << endl;
	cout << endl;

	cout << " Inferred topology : \t";
	if( min_score == dist_F1 ) cout << " 12|34 " << endl;
	else if( min_score == dist_F2 ) cout << " 13|24 " << endl;
	else cout << " 23|14 " << endl;
	cout << endl;
	
	return(0);
}
