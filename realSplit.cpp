// It writes the taxa names of the split as 1, 2, 3, 4, following the order determined in function realOrder.

#include "Clean_declaration.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;	//using namespace std;

vector<int> realSplit( vector<int> split, vector<int> order_taxa ) {
	int nb_split = split.size();
	vector <int> order_split( nb_split );
	
	for (int i=0; i < nb_split; ++i ) order_split[i] = order_taxa[ split[i] ];
	
	return order_split;
}