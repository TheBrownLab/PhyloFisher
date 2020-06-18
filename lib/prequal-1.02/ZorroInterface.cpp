/*
 * ZorroInterface.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

extern "C" {
#define EXTERN extern
#include "hmm.h"
#undef EXTERN
}

#include "Sequence.h"
#include "prequal.h"
#include "ZorroInterface.h"

using namespace::std;

extern COptions *options;			// Global options; would be better as a singleton, but this works

void InitialiseZorro(vector <CSequence> *cpp_seq) {
	Nseq = cpp_seq->size(); // Allocates the number of sequences in utils.h
	sequence = new char*[cpp_seq->size()]; // Allocates sequence array from utits.h
	lens = new int[cpp_seq->size()]; // Gets the memory for lengths from utils.h
	for (int i = 0; i < cpp_seq->size(); i++) {
		sequence[i] = new char[cpp_seq->at(i).length()];
		lens[i] = cpp_seq->at(i).length();
		for (int j = 0; j < cpp_seq->at(i).length(); j++) {
			if (cpp_seq->at(i).Seq(j).find("-") != std::string::npos) {
				cout << "\nERROR: sequence contains gaps and will break when pushed through a pairHMM...\n";
				exit(-1);
			}
			sequence[i][j] = pep2num(cpp_seq->at(i).Seq(j)[0]); // Do the weird character to int to character conversion that Zorro demands...
		}
	}
}

#define DEBUG_PP 0
// Run the HMM stuff if required
double ** RunHMM(vector <CSequence> *cpp_seq, string outFile, bool forceOverwrite) {
	double **retPP;
	assert(cpp_seq != NULL);
	// Initialise the memory for the PPs
	retPP = new double*[cpp_seq->size()]; // Initialise Posterior probabilities
	for (int i = 0; i < cpp_seq->size(); i++) {
		retPP[i] = new double[cpp_seq->at(i).length()];
		for (int j = 0; j < cpp_seq->at(i).length(); j++) { retPP[i][j] = 0.0; }
	}
	// Initialise the Zorro stuff
	InitialiseZorro(cpp_seq);

	// If file exists read the PPs
	if(file_exist(outFile) && !forceOverwrite) {
		cout << "\nPosterior probability file <" << outFile << "> exists. \n\tReading file";
		vector <string> Toks;
		string tmp;
		ifstream in(outFile);
		// Check sequence number
		tmp = read_line(in);
		Toks = Tokenise(tmp);
		if(Toks.size() != 1) { cout << "\nNumber of sequences on first line contains multiple tokens in PP file\n"; exit(-1); }
		if(atoi(Toks[0].c_str()) != cpp_seq->size()) { cout << "\nNumber of sequences in first line of PP file does not match data\n"; exit(-1); }
		cout << " ... checking names" << flush;
		// Check sequence names
		for(int i = 0 ;i < cpp_seq->size(); i++) {
			tmp = read_line(in);
			Toks = Tokenise(tmp);
			if(Toks.size() != 1) { cout << "\nName of sequence " << i << " has white spaces. This is not allowed.\n"; exit(-1); }
			if(Toks[0] != cpp_seq->at(i).Name()) { cout << "\nName of sequence " << i << " does not match that in data. This is not allowed\n"; exit(-1); }
		}
		cout << " ... checking PPs" << flush;
		// Input PPs
		double value;
		for(int i = 0 ;i < cpp_seq->size(); i++) {
			tmp = read_line(in);
			Toks = Tokenise(tmp);
			if(Toks.size() != cpp_seq->at(i).length()) { cout << "\nLength of sequence " << i << " does not match PPs. This is not allowed and can be caused by repeat masking.\n\tTry removing PP file?\n"; exit(-1); }
			for(int j = 0; j < cpp_seq->at(i).length(); j ++) {
				value = (double) atof(Toks[j].c_str());
				if(value < 0.0 || value > 1.0) { cout << "\nThe PPs for sequence " << i << " contains a non probability (" << value<< ")"; exit(-1); }
				retPP[i][j] = value;
			}
		}
		in.close();
		cout << " ... success" << flush;
	} else {
		// Calculate max length (repeat, but that's okay...)
		vector <int> lengthDist;
		for(int i = 0; i < cpp_seq->size(); i++) {
			lengthDist.push_back(cpp_seq->at(i).length());
		}
		cout << "\nPrepping pairHMM" << flush;
		initHMM(CSequence::MaxLength());
		cout << " ... done" << flush;
		///////////////////////////////// Get the PPs according to options ////////////////////////////////////"
		if(cpp_seq->size() - options->PPnumber() < 5 || options->AllPP()) { // Get all the posteriors
			cout << "\nCollecting all posterior probabilities. This may take some time for large data sets:\n";
			SimonGetPosteriors(CSequence::MaxLength(), retPP, options->DoApprox(), options->DefaultBound());
		} else if(options->ClosePP()) { // The closest
			cout << "\nCollecting subset of posterior probabilities based on closest " << options->PPnumber() << " sequences determined by Kmers\n\tThis may take time for larger data sets:\n" << flush;
			int rLsize = -1;
			MakeKmers(cpp_seq);
			cout << "\nCreating collection sets of PPs based on Kmer distances\n" << flush;
			// Calculate the median sequence length to ensure that some good sequence coverage is possible
			vector <int> lengths;
			for(auto &s : *cpp_seq) { lengths.push_back(s.length()); }
			sort(lengths.begin(),lengths.end());
			int median_length = lengths[lengths.size() / 2]; 							// The median length
			int reqOverMedian = (0.3 * my_min(options->PPnumber(),cpp_seq->size()));	// The number of sequences that must be over the median length
			// Now get the run list
			int **runList = new int*[cpp_seq->size()];
			for(int i = 0; i < cpp_seq->size(); i++) {
				int numOverMedian = reqOverMedian, pos = 0;
				vector <int> list = GetClosest(i, cpp_seq);
				int listSize = my_min(options->PPnumber(),list.size());
				ProgressSpinner(i + 1, cpp_seq->size());
				runList[i] = new int[listSize+1];
				if(rLsize < 0) { rLsize = listSize; }
				else if(rLsize != listSize) { cout << "\nProblem with different sized lists of running groups..."; }
				// Complex logic that picks the closest, while ensure at least reqOverMedian sequences are of at least median size
				for(int j = 0; j < list.size() && listSize > 0; j++) {
					if(numOverMedian <= 0 || listSize > numOverMedian || cpp_seq->at(list[j]).length() >= median_length || options->PPnumber() > list.size()) {
						if(cpp_seq->at(list[j]).length() >= median_length) { numOverMedian --; }
						runList[i][pos++] = list[j];
						listSize--;
					}
				}
				if(listSize != 0) {
					cout << "\nERROR: Failed to build a suitable set of PPs to sample for sequence["<<i<<"] " << cpp_seq->at(i).Name() << "\n\n";
					cout << "\nList["<<listSize <<"]: "; for(auto l : list) { cout << ":" << l << " "; }
					exit(-1);
				}
			}
			assert(rLsize > 0);
			cout << " ... done\nGetting posterior probabilities:\n";
			SimonGetSubsetPosteriors(CSequence::MaxLength(), retPP, runList, rLsize, options->DoApprox(), options->DefaultBound());
			for(int i = 0; i < cpp_seq->size(); i++) { delete [] runList[i]; } delete [] runList;
		} else if(options->LongPP()) { // The longest
			cout << "\nCollecting subset of posterior probabilities based on longest " << options->PPnumber() << " sequences:\n" << flush;
			int *coreList = NULL;
			vector <int> newL = ordered(lengthDist);
			coreList = new int[options->PPnumber()];
			for(int i = 0; i < options->PPnumber(); i++) {
				coreList[i] = newL[cpp_seq->size() - 1 - i];
			}
			int **runList = new int*[cpp_seq->size()];
			for(int i =0 ; i < cpp_seq->size(); i++) {
				runList[i] = coreList;
			}
			SimonGetSubsetPosteriors(CSequence::MaxLength(), retPP, runList, options->PPnumber(), options->DoApprox(), options->DefaultBound());
			for(int i =0 ; i < cpp_seq->size(); i++) { runList[i] = NULL; }
			delete [] runList;
			delete [] coreList;
		} else {
			cout << "\nUnknown option for calculating posteriors...\n"; exit(-1);
		}
		// Output PPs
		if (options->DoPPs()) {
			cout << "\nOutputting posterior probabilities to " << outFile
					<< flush;
			ofstream out(outFile.c_str());
			out << cpp_seq->size();
			for (int i = 0; i < cpp_seq->size(); i++) {
				out << "\n" << cpp_seq->at(i).Name();
			}
			for (int i = 0; i < cpp_seq->size(); i++) {
				out << "\n";
				for (int j = 0; j < cpp_seq->at(i).length(); j++) {
#if DEBUG_PP == 1
					out <<"[" << i << ","<< j<<"]";
#endif
//					if(!InRange(retPP[i][j],-0.001,1.001)) { cout << "\nERROR: obtained a PP of " << retPP[i][j] << " for residue " << j << " in sequence [" << i << "]  " << cpp_seq->at(i).Name(); exit(-1); }
					if(retPP[i][j] > 1) { retPP[i][j] = 1.0; }
					else if(retPP[i][j] < 0) { retPP[i][j] = 0.0; }
					out << retPP[i][j] << "\t";
				}
			}
			out << "\n";
			out.close();
			cout << "... done" << flush;
		}
	}
	// Note: There's a lot of memory unallocated here, but that's Zorro's problem for now...

	return retPP;
}


//////////////////////////////////////// Kmer stuff from PaHMM tree //////////////////////////////////////////
vector<unordered_map<string,short>*>* kmers;
unsigned int kmerSize = 4;

void MakeKmers(vector <CSequence> *data) {
	assert(data != NULL);
	// Get the memory
	kmers = new vector<unordered_map<string,short>*>(data->size());

	for(int i = 0; i < data->size(); i++) {
		(*kmers)[i] = new unordered_map<string,short>();
		string seq = data->at(i).Seq();
		extractKmers(seq, (*kmers)[i]);
	}
}
void extractKmers(string& seq, unordered_map<string, short>* umap)
{
	string kmer;

	if(seq.size() < kmerSize)
		return;

	for(unsigned int i = 0; i< (seq.size() - kmerSize); i++)
	{
		kmer = seq.substr(i, kmerSize);
		++((*umap)[kmer]);
	}
}

vector <int> GetClosest(int Sequence, vector <CSequence> *data) {
	vector <double> dist = MakeKmerDistanceRow(Sequence, data);
	vector <int> indices = ordered(dist);
	for(int i = 0; i < indices.size(); i++) {
		if(indices[i] == Sequence) {
			indices.erase(indices.begin() + i);
			break;
		}
	}
	return indices;
}

vector <double> MakeKmerDistanceRow(int CurSeq, vector <CSequence> *data) {
	vector <double> RetDist;
	for(int i = 0; i < data->size(); i++) {
		if(i == CurSeq) { RetDist.push_back(1.0); continue; }	// For the current sequence return error
		double identity = 1.0 - commonKmerCount(i,CurSeq)/(double)my_min(data->at(i).length(),data->at(CurSeq).length());
		identity = (pow(100, identity-1.04)+0.01)/0.6;	// Marcin's AA function
		RetDist.push_back(identity);
	}
	return RetDist;
}

unsigned int commonKmerCount(unsigned int i, unsigned int j)
{
	unordered_map<string, short>* m1 = (*kmers)[i];
	unordered_map<string, short>* m2 = (*kmers)[j];
	unsigned int commonCount = 0;

	for(auto it = m1->begin(); it != m1->end(); it++)
	{
		commonCount += std::min((short)((*m2)[it->first]), (short)(it->second));
	}
	//DEBUG("Common k-mer count between seq. " << i << " and " << j << " is " << commonCount);
	return commonCount;


}

