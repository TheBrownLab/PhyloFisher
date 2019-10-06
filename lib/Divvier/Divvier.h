/*
 * Divvier.h
 *
 *  Created on: 24 Oct 2017
 *      Author: simon
 */

#ifndef DIVVIER_H_
#define DIVVIER_H_

#include <string>
#include <vector>
#include <sstream>
#include "Tree.h"	// Mucky dependency, but this has the tools in it

const string VERSION_NUMBER = "1.0 (release)";		// Version number

CTree MakeTree(std::vector <std::string> Names, std::vector <std::string> Seqs);	// Make a distance based guide tree
double GetPercentDiff(std::string seq1, std::string seq2);
double AAJCdist(double p);
bool IsGap(char c);

std::string DoBioNJ(std::vector <double> PWdists, std::vector <std::string> Names, bool DoNumbers = false);
bool CheckNext(int i, int argc, char *argv[]);	// For options

void GetPosteriors(std::string File);

// Class for the posterior probabilities
class CPostP {
public:
	CPostP(int X, int Y, double *postP, int len) {
		_X = X; _Y = Y;
		for(int i = 0; i < len; i++) {
			_PP.push_back(my_max(0.0, my_min( 1.0 , postP[i] )));
		}
	}
	CPostP(std::string input) {		// Reads in same format as out
		std::vector<std::string> Toks = Tokenise(input);
		if(Toks.size() < 6) { cout << "\nError in CPostP::CPostP(input) -- line needs at least 1 posterior probability: "<< input << endl; exit(-1); }
		if(Toks[0][0] != '[' || Toks[4][0] != ']') { cout << "\nSuspect opening to posterior probability line in CPostP::CPostP(input): " << input << endl; exit(-1); }
		_X = atoi(Toks[1].c_str()); _Y = atoi(Toks[3].c_str());
		assert(_X >= 0 && _Y >= 0);
		for(int i = 5; i < Toks.size() ; i++)  {
			_PP.push_back(my_max(0.0, my_min( 1.0 , atof(Toks[i].c_str()))));
		}
	}
	std::string out() {
		std::stringstream ss;
		ss <<"[ " << _X << " , " << _Y <<" ]";
		for(int i = 0 ; i < _PP.size(); i++) { ss<< "\t" << _PP[i]; }
		return ss.str();
	}
	int x() { return _X; }
	int y() { return _Y; }
	double PP(int i) { return _PP[i]; }

private:

	int _X, _Y;					// The sequences
	std::vector<double> _PP;		// The posterior probabilities
};


#endif /* DIVVIER_H_ */
