Remember to install armadillo and openBlas to be able to compile

#ssh hercules05
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/lib/armadillo-8.6rc2/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/soft/lib/openBlas2.20

1) Compile the different file.o scripts based on its file.cpp
The files Clean_declaration.h (declare the functions) and typedefs.h are needed:
g++ -std=c++11 -c constructAlph.cpp distFrobenius.cpp flatteningMat.cpp getColumns.cpp readFASTA.cpp realOrder.cpp realSplit.cpp -O1 -larmadillo

2) Compiling the main script (SVD.cpp) based on the other functions
g++ -std=c++11 SVD.cpp constructAlph.o distFrobenius.o flatteningMat.o getColumns.o readFASTA.o realOrder.o realSplit.o  -o SVD -O1 -larmadillo

3) To run the SVD just indicate the input fasta file, the characters that you want to exclude (like gaps: -) and the number of partitions ()
./SVD fasta_file.fa characters_to_exclude nb_partitions

	Example: ./SVD example.fa - 3
