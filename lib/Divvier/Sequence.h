/*
 * Sequence.h
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

const std::string AA_ABET  = "ARNDCQEGHILKMFPSTWYV";
#define BIG_NUMBER 1000000

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
	bool Filter(int pos);					// Whether pos should be filtered/removed in any way

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
	std::string _name;			// The sequence
	std::string _seq;			// The name
	static char _filterOut;		// The string output on filtering
	bool _allRemoved = false;	// Whether the sequence is fully removed

	void InitialiseFlags() {
		assert(Inside.empty() && Remove.empty());
		Inside.assign(_seq.size(),true);
		Remove.assign(_seq.size(),false);
	}

};

// File readers
enum EFileType { FASTA, MSF, Phylip, Interleaved };
inline std::string FileTypeName(EFileType type) {
	switch(type) {
	case FASTA:
		return "FASTA";
	case MSF:
		return "MSF";
	case Phylip:
		return "Phylip";
	case Interleaved:
		return "Interleaved";
	default:
		std::cout << "\nUknown FileTypeName..."; exit(-1);
	}
}
EFileType TestFile(std::string seqFile);
std::vector <CSequence> *ReadSequences(std::string seqFile);
std::vector <CSequence> *FASTAReader(std::string seqFile);
std::vector <CSequence> *MSFReader(std::string seqFile);
std::vector <CSequence> *PhylipReader(std::string seqFile);
std::vector <CSequence> *InterleavedReader(std::string seqFile);

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

inline bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

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

inline bool IsGap(char c) {
	std::string gaps = ".*-X?";
	if(std::find(gaps.begin(),gaps.end(),c) != gaps.end()) { return true; }
	return false;
}

inline bool IsSeq(char c) {
	if(std::find(AA_ABET.begin(),AA_ABET.end(),c) != AA_ABET.end()) { return true; }
	return false;
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
