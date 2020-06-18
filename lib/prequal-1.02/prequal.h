/*
 * prequal.h
 *
 *  Created on: 12 Sep 2017
 *      Author: simon
 */

#ifndef SEQFILTER_H_
#define SEQFILTER_H_

#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <numeric>
#include <regex>
#include <cstring>
#include <math.h>


/////////////// General naming stuff
const std::string PROGRAM_NAME = "PREQUAL";
const std::string VERSION_NUMBER = "1.02";
#define DEFAULT_THRESHOLD 0.994

extern std::stringstream warningStream;

// Class holding the options for the program. Implemented in Options.cpp
class COptions {
public:
	COptions(int argc, char * argv[]);			// Creates the options from reading the command line. Only done once

	// Access functions
	// 1. File handling
	std::string Infile() { return _inputFile; }
	std::string OutSuffix() { return _outputSuffix; }
	bool DoSummary() { return _doSummary; }
	std::string SummarySuffix() { return _summarySuffix; }
	bool DoDetail() { return _doDetail; }
	std::string DetailSuffix() { return _detailSuffix; }
	bool DoPPs() { return _doPPs; }
	std::string PPSuffix() { return _ppSuffix; }
	bool Overwrite_PP() { return _overwritePP; }
	// 2. Core regions
	int RunBeforeInside() { return _runBeforeInside; }
	bool Remove2Core() { return _remove2Core; }
	bool RemoveAll() { return _removeAll; }
	char CoreFilter() {
		if(_removeAll) {return '\0'; }		// RemoveAll is equivalent of a null for the CoreFilter character
		return _coreFilter;
	}
	// 3. Filtering options
	bool IgnoreSequence(std::string seq_name) {
		if(_noFilterList.empty() && _noFilterWord.empty()) { return false; }
		if(find(_noFilterList.begin(), _noFilterList.end(),seq_name) != _noFilterList.end()) { return true; }
		for(auto &w : _noFilterWord) {
			if(seq_name.find(w) != std::string::npos) { return true; }
		}
		return false;
	}
	bool DoKeepProportion() { if(_filterThreshold < 0) { return true; } else { return false; } }
	double KeepProportion() { return _keepProportion; }
	double KeepThreshold() { return _filterThreshold; }
	int FilterRange() { return _joinFilterRange; }
	// Extended HMM filtering
	bool DoHMMpostFilter() { return _doHMMpostFilter; }
	bool ExpErrorProb() { return _expErrorProb; }
	bool ExpErrorLength() { return _expErrorLength; }
	// 4. PP calculation options
	int PPcalcNumber() { return _ppCalcs; }
	bool AllPP() { if(_ppCalcs == -1) { return true; } return false; }
	bool LongPP() { if(_ppCalcs == 1) { return true; } return false; }
	bool ClosePP() { if(_ppCalcs == 0) { return true; } return false; }
	int PPnumber() { return _ppCalcNumber; }
	// 5. Approximations for pairHMM calculations
	bool DoApprox() { return _approxPP; }
	int DefaultBound() { return 25; }
	// 6. Repeat handling
	bool RemoveRepeat() { return _removeRepeat; }
	int RepeatLength() { return _repeatLength; }
	// 7. DNA handling
	bool AllowDNA() { return _allowDNA; }
	bool AlwaysUniversal() { return _alwaysUniversal; }
private:
	// Behaviour defining core and how the different parts are treated/filtered
	int _runBeforeInside = 3;			// Number of positively defined homologies before a sequence is consider part of the core
	bool _remove2Core = true;			// Whether parts outside the core are removed
	bool _removeAll = false;			// Whether all fails are removed
	char _coreFilter = 'X';			// The character used for filtering in the core
	bool _removeRepeat = true;			// Whether repeats are classified as removed prior to analysis
	int _repeatLength = 20;				// The default length to specify a repeat
	// DNA handling
	bool _allowDNA = true;				// Whether DNA sequences are allowed or not
	bool _alwaysUniversal = false;      // Whether to always use the universal genetic code for translations
	// File handling
	std::string _inputFile;
	std::string _outputSuffix = ".filtered";
	bool _doSummary = false;				// General summary
	std::string _summarySuffix = ".summary";
	bool _doDetail = false;				// Detailed output per site
	std::string _detailSuffix = ".detail";
	bool _doPPs = true;					// Grid of PPs for fast downstream computation
	std::string _ppSuffix = ".PP";
	// Filtering options
	std::vector <std::string> _noFilterList;		// List of taxa names that will not be filtered at all
	std::vector <std::string> _noFilterWord;		// List of taxa regular expressions that will not be filtered at all
	double _keepProportion = 0.95;		// The proportion of residues to be kept
	double _filterThreshold = DEFAULT_THRESHOLD;		// Hard set the filter threshold
	int _joinFilterRange = 10;			// The maximum gap between filtered characters before the whole segment is filtered
	// post pairHMM length filtering based on another HMM
	bool _doHMMpostFilter = false;
	double _expErrorProb = 0.001;			// Expected probability of errors
	int _expErrorLength = 10;			// Expected error length
	// Posterior probability calculation options
	int _ppCalcs = 0;					// Options are 0: closest _ppCalcNumber; 1: longest _ppCalcNumber; -1: all
	int _ppCalcNumber = 10;
	// Options that are currently disabled
	bool _overwritePP = false;
	bool _approxPP = true;				// Whether to use full or approximate PP computations according to my hack
	// Internal functions
	std::vector <std::string> GetStringList(std::string FilterListFile);
};

double mean(std::vector <double> vec);
double stdev(std::vector <double> vec);
double TargetCutoff(double prop2Keep, std::ostream &os = std::cout);	// Computes the target cutoff given a specific proportion to keep
void DoFiltering(double Threshold);										// Applies the filtering method

#endif /* SEQFILTER_H_ */
