/*
 * ZorroInterface.h
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 *
 *  // Isolates all the mess in Zorro globals the here
 */

#ifndef ZORROINTERFACE_H_
#define ZORROINTERFACE_H_

extern "C" {
#define EXTERN extern
#include "hmm.h"
#undef EXTERN
}

void InitialiseZorro(std::vector <CSequence> *cpp_seq);		// Get the memory that Zorro wants
double ** RunHMM(std::vector <CSequence> *cpp_seqs, std::string outFile, bool forceOverwrite = false);		// C++ interaction function with Zorro


// Kmer stuff from PaHMM-Tree by Marcin Bogusz (https://github.com/marbogusz/paHMM-Tree)
void MakeKmers(std::vector <CSequence> *cpp_seqs);
void extractKmers(std::string& seq, std::unordered_map<std::string,short>* umap);
std::vector <double> MakeKmerDistanceRow(int CurSeq, std::vector <CSequence> *data);
unsigned int commonKmerCount(unsigned int i, unsigned int j);
std::vector <int> GetClosest(int Sequence, std::vector <CSequence> *data);


#endif /* ZORROINTERFACE_H_ */
