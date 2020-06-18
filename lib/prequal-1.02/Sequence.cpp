/*
 * Sequence.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

#include "Sequence.h"

int CSequence::_maxLength = 0;
char CSequence::_filterOut = 'X';

using namespace::std;

//////////////// CSequence
CSequence::CSequence(std::string name, std::string seq) {
	AddName(name);
	AddSequence(seq);
	InitialiseFlags();
}
void CSequence::AddName(std::string name) {
	if(!_name.empty()) {
		std::cout << "\nCSequence ERROR: Trying to added name to non-empty named sequence\n";
		exit(-1);
	}
	_name = name;
}
void CSequence::AddSequence(std::string seq) {
	if (!_seq.empty()) {
		std::cout << "\nCSequence ERROR: Trying to added sequence to non-empty named sequence\n";
		exit(-1);
	}
	_seq = seq;
	if(_seq.size() > _maxLength) { _maxLength = _seq.size(); }
}
std::string CSequence::Seq(int pos, bool filter, bool showOutside) {
	std::stringstream ss;
	if(pos != -1) {
		if(!showOutside && !Inside[pos]) { ss << "0"; }
		else if(filter && Remove[pos]) { if(_filterOut != '\0') { ss <<_filterOut; } }
		else { ss << _seq[pos]; }

	} else {
		for( int i = 0 ; i < _seq.size(); i++) {
			if(!showOutside && !Inside[i]) { continue; }
			if(filter && Remove[i]) { if(_filterOut != '\0') { ss <<_filterOut; } }
			else { ss << _seq[i]; }
		}
	}
	return ss.str();
}
std::string CSequence::RealSeq(int pos) {
	if(pos != -1) {
		return _seq.substr(pos,1);
	}
	return _seq;
}
bool CSequence::Filter(int pos) {
	if(Remove[pos] || !Inside[pos]) { return true; }
	return false;
}
string CSequence::RealDNA(int pos) { // Real codon at position pos
	if(pos == -1) { return _dna_seq; }
	if(pos > length()) { cout << "\nPosition must be within amino acid sequence for RealDNA output"; }
	return _dna_seq.substr(pos*3,3);
}
string CSequence::DNA(int pos, bool filter, bool showOutside) { // Filtered codon at position pos
	std::stringstream ss;
	if(pos != -1) {
		if(pos > length()) { cout << "\nPosition must be within amino acid sequence for RealDNA output"; }
		if(!showOutside && !Inside[pos]) { ss << "0" << "0" << "0"; }
		else if(filter && Remove[pos]) { if(_filterOut != '\0') { ss <<_filterOut << _filterOut << _filterOut; } }
		else { ss << _dna_seq.substr(pos*3,3); }
	} else {
		for( int i = 0 ; i < length(); i++) {
			if(!showOutside && !Inside[i]) { continue; }
			if(filter && Remove[i]) { if(_filterOut != '\0') { ss <<_filterOut << _filterOut << _filterOut; } }
			else { ss << _dna_seq.substr(i*3,3); }
		}
	}
	return ss.str();
}


// Function that examines a sequence and looks for self repeats
// ---
// This might be better served by a suffix tree implementation, but I'm just doing the dumb thing for now
bool CSequence::CleanRepeat(int repeatLength) {
	string seq = _seq;
	string dna = _dna_seq;
	bool foundRepeat = false;
	int end;
	for(int pos = 0; pos < (int) seq.size() - repeatLength; pos ++) {
		string repeat = seq.substr(pos,repeatLength);
		// Repeat can never contain X
		if(find(repeat.begin(),repeat.end(),'X') != repeat.end()) { continue; }
		size_t next_pos= pos + repeatLength;
		while(next_pos != string::npos) {
			next_pos = seq.find(repeat,next_pos);
			if(next_pos == string::npos) { break; }
			foundRepeat=true;
			for(end = 0; end < (int) seq.size() - repeatLength - next_pos; end++) { // Extend the repeat while there is perfect match
				if(seq[pos + end + repeatLength] != seq[next_pos + end + repeatLength]) { break; }
			}
			seq.erase(next_pos,repeatLength + end);
			if(!dna.empty()) {
				dna.erase(next_pos*3,(repeatLength + end)*3);
			}
		}
	}
	if(foundRepeat) { _seq = seq; _dna_seq = dna; }
	return foundRepeat;
}

// Translation functions
bool CSequence::MakeTranslation(bool forceUniversal) {
	// Check the sequence is in triplets and suitable for translation
	if(length() % 3 != 0) { cout << "\nTrying to translate " << Name() << ", but not in a multiple of three\n\n"; exit(-1); }

	if(forceUniversal) {
		if(TryTranslation(0,true)) { return true; }
	} else {
		for(int i = 0; i < NumGenCode; i++) {
			if(TryTranslation(i)) { return true; }
		}
	}
	return false;
}
bool CSequence::TryTranslation(int genCode, bool force) {
	if(_genCode != -1) { cout << "\nTrying translation when sequence is already translated\n"; exit(-1); }
	string codon;
	string aa_seq;
	bool okay;
	int cod_num;
	for(int i = 0 ; i < length(); i+=3) {
		codon = _seq.substr(i,3);
		okay = true;
		for(auto &s : codon) {
			if(s == '?' || toupper(s )== 'X' || toupper(s) == 'N') { // Checks for ambiguity characters, but may miss some weirder ones
				okay = false; aa_seq.push_back('X'); break;
			} else if(IsGap(s)) { // If it's a gap then push out a gap
				okay = false; aa_seq.push_back('-'); break;
			}
		} // Check whether codon can be translated
		if(!okay) { continue; }
		cod_num = GenCodes[genCode][GetCodon(codon)];
		if(cod_num == -1) {
			if(force) { aa_seq.push_back('X'); continue; }
			else { // Last codon or force allowed to be a stop codon
				if(i + 3 >= length()) {
					_seq.erase(i,3);
					warningStream << "\nWARNING: " << Name() << " has has stop codon at end of sequence removed";
					break;
				}
				else { return false; }
			}
		}
		aa_seq.push_back(AA_ABET[cod_num]);
	}
	_dna_seq = _seq;
	_seq = aa_seq;
	_genCode = genCode;
	return true;
}


/////////////// Minor functions
// File reader
std::vector <CSequence> *FASTAReader(std::string SeqFile, bool forceUniversal) {
	std::vector <CSequence> *RetSeq = new std::vector<CSequence>();

	std::ifstream input(SeqFile.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }

    std::string line, name, content;

    while(getline( input, line ) ){
    	line = RemoveWhiteSpace(line);
    	if(line.empty()) { continue; }		// Skip blank lines
        if(line[0] == '>' ){ // Identifier marker
            if( !line.empty() ){
            	if(!name.empty()) { // Push the sequence back
            		RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
            		name.clear(); content.clear();
            	}
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            content += RemoveWhiteSpace(line);
        }
    }
    // Add the final sequence
    RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
    // Check whether all DNA
    bool allDNA = true;
    for(auto & seq : *RetSeq) {
    	for(auto & s : seq.Seq()) {
    		if(IsGap(s)) { continue; }
    		if(!IsDNA(s)) { allDNA = false; break; }
    	}
    	if(!allDNA) { break; }
    }
    if(allDNA) {
    	cout << "\nFound only DNA sequences. Doing translations.";
    	for(auto & seq : *RetSeq) {
    		if(!seq.MakeTranslation(forceUniversal)) {
    			cout << "\nFound DNA sequences, but cannot find a successful translation... abandoning!";
    			cout << "\nSequence failed: " << seq.Name() << "\n\n"; exit(-1);
    		}
    	}
    }
    return RetSeq;
}

std::string RemoveWhiteSpace(std::string s) {
	s.erase( std::remove_if( s.begin(), s.end(), ::isspace ), s.end() );
	return s;
}

std::vector <std::string> Tokenise(std::string line) {
	std::string buf;
	std::stringstream in(line);
	std::vector <std::string> Toks;
	Toks.~vector();
	while(in >> buf) { Toks.push_back(buf); }
	return Toks;
}
std::vector <std::string> Tokenise(std::string line, std::string Delim)	{
	size_t i = 0, j,j1;
	std::vector <std::string> Toks;
	while(i != (int)line.size())	{
		j = line.find(Delim,i+1);
		if(j == std::string::npos) { j = j1 = (int)line.size(); } else { j1 = j+1; }
		Toks.push_back(line.substr(i,j-i));
		i = j1;
	}
	return Toks;
}

