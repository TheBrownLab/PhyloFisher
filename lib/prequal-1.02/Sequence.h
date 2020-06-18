 /* Sequence.h
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

#ifndef SEQUENCE_H_
#define SEQUENCE_H_
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <numeric>
#include <regex>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

extern std::stringstream warningStream;

// Alphabet checker
inline bool IsABET(char &c, const std::string &ABET) { if(find(ABET.begin(),ABET.end(),toupper(c)) == ABET.end()) { return false; } return true; };
// Gap alphabet
const std::string GAP_ABET = "-X?"; 		// The list of gaps
inline bool IsGap(char &c) { return IsABET(c, GAP_ABET); };
// DNA alphabet
const std::string DNA_ABET = "ACGTN";	// Keep n since it's common
inline bool IsDNA(char &c) { return IsABET(c, DNA_ABET); };
 // Amino acid Alphabet
 const std::string AA_ABET  = "ARNDCQEGHILKMFPSTWYV";
inline bool IsAA(char &c) { return IsABET(c, AA_ABET); };
// Codon alphabet
const std::string COD_ABET = "AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT";
inline int GetCodon(std::string codon) {
	int pos;
	if(codon.size() != 3) { std::cout << "\nWrong codon size. Developer problem?\n"; exit(-1); }
	for(pos = 0; pos < COD_ABET.size(); pos+=3) {
		if(COD_ABET[pos+0] == toupper(codon[0]) && COD_ABET[pos+1] == toupper(codon[1]) && COD_ABET[pos+2] == toupper(codon[2])) { break; }
	}
	if(pos == COD_ABET.size()) { std::cout << "\nUnknown codon. Is it DNA? Developer problem?\n"; exit(-1); }
	return pos/3;
};
// Genetic codes
 // Definitions of the genetic code
 const int NumGenCode = 11;
 // The genetic codes are:
 //	0:  Universal code
 //	1:  Vertebrate mt
 //	2:  Yeast mt
 //	3:  Mould mt
 //	4:  Invertebrate mt
 //	5:  Ciliate nuclear
 //	6:  Echinoderm mt
 //	7:  Euplotid mt
 //	8:  Alternative yeast nuclear
 //	9:  Ascidian mt
 //  10: Blepharisma nuclear
 //	11: Fake code where everything codes
 const std::string GenCodeName[] = {
 		"Universal",					// [0]
 		"Vertebrate mt",				// [1]
 		"Yeast mt",						// [2]
 		"Mould mt",						// [3]
 		"Invertebrate mt",				// [4]
 		"Ciliate nuclear",				// [5]
 		"Echinoderm mt",				// [6]
 		"Euplotid mt",					// [7]
 		"Alternative yeast nuclear",	// [8]
 		"Ascidian mt",					// [9]
 		"Blepharisma",					// [10]
 		"Fake with everything coding"	// [11]
 };

 const int GenCodes[][64] = { // Universal
 	{	11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
 	},{ // Vert mt
 		11, 2,11, 2,16,16,16,16,-1,15,-1,15,12, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{ // Yeast mt
 		11, 2,11, 2,16,16,16,16,-1,15,-1,15,12, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{ // Mould mt
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15,12, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,16,16,16,16,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
 	},{ // Intvertebtate mt
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{ // Ciliate nuclear
 		11, 2,11, 2,16,16,16,16,15,15,15,15,12, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		 5,18, 5,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
 	},{
 		 2, 2,11, 2,16,16,16,16,15,15,15,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15, 4, 4,17, 4,10,13,10,13
 	},{
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,15,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
 	},{
 		11, 2,11, 2,16,16,16,16, 7,15, 7,15,12, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18,-1,18,15,15,15,15,17, 4,17, 4,10,13,10,13
 	},{
 		11, 2,11, 2,16,16,16,16, 1,15, 1,15, 9, 9,12, 9,
 		 5, 8, 5, 8,14,14,14,14, 1, 1, 1, 1,10,10,10,10,
 		 6, 3, 6, 3, 0, 0, 0, 0, 7, 7, 7, 7,19,19,19,19,
 		-1,18, 5,18,15,15,15,15,-1, 4,17, 4,10,13,10,13
 	},{
 		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 	     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 		 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
 }	};


// Basic class for sequences
class CSequence {
public:
	CSequence() { };								// Blank constructor
	CSequence(std::string name, std::string seq);	// Standard constructor
//	std::string name;
//	std::string seq;
	std::vector <bool> Inside;					// Whether the character is on the inside or outside
	std::vector <bool> Remove;					// Whether to remove in the filter (true = remove)
	double PropInside;							// The proportion of the sequence labeled inside
	double PropRemoved;							// The proportion of hte sequence labeled to be removed
	bool AllRemoved() { return _allRemoved; }
	void AddSequence(std::string seq);
	void AddName(std::string name);
	static void SetFilter(char filterOut) { _filterOut = filterOut; };
	int length() { return _seq.size(); }
	static int MaxLength() { return _maxLength; }
	std::string RealSeq(int pos = -1);				// Outputs the unfiltered seq
	std::string Seq(int pos = -1, bool filter = true, bool showOutside = false);		// Output the sequence (or pos i of sequence)
	std::string Name() { return _name; }
	bool Filter(int pos);						// Whether pos should be filtered/removed in any way
	// Checks whther a sequence has a repeat; Filter = whether to set for removal; repeatLength = number of characters specified for repeat
	bool CleanRepeat(int repeatLength = 20);
	// Resets the maxLength if needed
	static void ResetMaxLength(std::vector <CSequence> *seqs) {
		_maxLength = 0;
		for(auto &seq : *seqs) { if(seq.length() > _maxLength) { _maxLength = seq.length(); } }
	};
	// Handle translations
	bool MakeTranslation(bool forceUniversal = false);
	bool HasDNA() { if (_genCode == -1) { return false; } return true; }
	std::string RealDNA(int pos = -1);	// Matches the RealSeq command, but for DNA (codons)
	std::string DNA(int pos = -1, bool filter = true, bool showOutside = false); // Matches the Seq command but for DNA (codons)
	std::string GenCode() { if(_genCode == -1) { return "Not translated"; } return GenCodeName[_genCode]; }
	// Output stuff
	std::string out() { return _name + " " + _seq; }
	void CalculateSummary() {
		int in = 0, rem = 0;
		for(int i = 0; i < length(); i++) {
			if(Inside[i]) { in++; }
			if(Remove[i]) { rem++; }
		}
		PropInside = (double) in / (double) length();
		PropRemoved = (double) rem / (double) length();
		if(rem == length()) { _allRemoved = true; }
	}
private:
	static int _maxLength;		// Maximum length of the sequences examined
	std::string _name;			// The name
	std::string _seq;			// The amino acid sequence
	int _genCode = -1;			// The genetic code the translation has
	std::string _dna_seq;		// The dna seq (if it exists)
	static char _filterOut;		// The string output on filtering
	bool _allRemoved = false;	// Whether the sequence is fully removed

	void InitialiseFlags() {
		assert(Inside.empty() && Remove.empty());
		Inside.assign(_seq.size(),true);
		Remove.assign(_seq.size(),false);
	}

	bool TryTranslation(int genCode, bool force = false);	// Translate to genetic code genCode, if force then internal stops set to X

};

// File reader
std::vector <CSequence> *FASTAReader(std::string SeqFile, bool forceUniversal = false);

// Other minor tools
template <class TRange> bool InRange(TRange Val, TRange LowerBound, TRange UpperBound) { return ( !(Val < LowerBound) && ( Val < UpperBound) ); }
#define my_min(a,b) ((a)<(b)?(a):(b))
#define my_max(a,b) ((a)>(b)?(a):(b))
std::string RemoveWhiteSpace(std::string s);
std::vector <std::string> Tokenise(std::string line);	// Tokenise a string
std::vector <std::string> Tokenise(std::string line, std::string Delim);		// Tokenise a string according to delimiter Delim
inline void ProgressSpinner(int suffix = -1,int suffix_total = -1,std::string prefix = "") {
        static int count = 0;
        static char progress_spinner [] = "/-\\|";
        std::cout << "\r" << prefix << progress_spinner[count++];
        if(suffix >= 0) { std::cout << " " << suffix; }
        if(suffix_total >= 0) { std::cout << " / " << suffix_total; }
        std::cout << std::flush;
        if(count == 4) { count = 0; }
};
inline bool file_exist (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

inline std::string read_line(std::istream &in) {
	std::string tmp;
	getline(in,tmp);
	if(!in.good()) { std::cout << "\nError reading file..."; exit(-1); }
	return tmp;
}

template <typename T>
std::vector<int> ordered(std::vector<T> const& values) {
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));

    std::sort(
        begin(indices), end(indices),
        [&](int a, int b) { return values[a] < values[b]; }
    );
    return indices;
}

#endif /* SEQUENCE_H_ */
