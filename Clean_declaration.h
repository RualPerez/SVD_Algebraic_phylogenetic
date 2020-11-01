#if ! defined (ALIGNUTILS_H) 
#define ALIGNUTILS_H

#include "typedefs.h"
#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <armadillo>

using namespace arma;

struct Alignment {
	unsigned int num_taxa;
	unsigned int seq_len;
	vector<string> taxa; // the taxa names
   map<string,string> seqs;  // the sequences
  //	unsigned int num_states;
//	map<char,int> table;
};

struct CompareSecond {
	bool operator()(const MyPairType& left, const MyPairType& right) const {
		return left.second < right.second;
	}
};

//................................................................................................................................
// Clean_functions.cpp
//................................................................................................................................

// AUXILIARS
vector<int> powers_K(int numb_sp, int nb_alph);
vector<string> makepatterns(vector<string> &patterns);
vector<string> createpatterns(int n);

// MATRICES & TENSORS
mat projection(mat A, int k);
map <string,float> matrix2tensor(mat A); 
float frobnorm(mat spl_mat);
map <string,float> prob_teor(string filename);
int write_relfreq(string filename, map<string,float> relfreq);
map <string,float> read_freqs(string filename, int i);

// ALIGNMENTS  
Alignment readFASTA (string fname);
bool find_chr_exclude(string chr_exclude, char to_find); 	
map<string,int> getColumns(Alignment &align, string chr_exclude);
map<char, int> constructAlph (Alignment &align, string chr_exclude);
map<char,int> Freq_chr (Alignment &align);
vector<int> realOrder(vector<string> taxa);
vector<int> realSplit(vector<int> split, vector<int> order_taxa);

// NOT_USED
Alignment gen_align(int par_a,int par_b,int par_num, int length);
void printAlign (Alignment const & align);

// ERIK+2 
mat flatteningMat(map<string,int> tensor, map<char,int> alph, vector<int> split, vector <int> othersp, unsigned int seq_len);
mat Flattening_time(map<string,int> tensor, map<char,int> alph, vector<int> split, vector <int> othersp);
double Erik2(mat flat, int n_partitions, int nb_alph, vector<int> split, vector<int> count_split);
double Erik(mat flat, int n_partitions, int nb_alph);
vector<long int> Decomp(string it, vector <int> split, vector <int> othersp); 
double distFrobenius( mat flat, int nb_partitions, int nb_alph);

// CONVERSIONS
string int2string(int number);
string longint2string(long int number);
string float2string(float number);
string treu_coma(string seq);
int nuc2int(string x);
int cnuc2int(char x);
long int translator(string x);
pair<char, int> min_pair(map<char,int> mymap);

// VOIDS
void Print_Freq_sort(map<char,int> mymap, unsigned int all_freq);
void construct_tree(vector<string> taxa,  vector <int> set_species,  vector <int> other_set);

#endif


