/*
 * Sequence.cpp
 *
 *  Created on: Oct 13, 2017
 *      Author: simon
 */

#include "Sequence.h"
#include <cstdlib>

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
		else if(filter && Remove[pos]) { ss <<_filterOut; }
		else { ss << _seq[pos]; }

	} else {
		for( int i = 0 ; i < _seq.size(); i++) {
			if(!showOutside && !Inside[i]) { continue; }
			if(filter && Remove[i]) { ss <<_filterOut; }
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


/////////////// Minor functions
// File reader
std::vector <CSequence> *ReadSequences(std::string seqFile) {
	vector <CSequence> *ret = NULL;
	switch(TestFile(seqFile)) {
	case FASTA:
		ret = FASTAReader(seqFile); break;
	case MSF:
		ret = MSFReader(seqFile); break;
	case Interleaved:
		ret = InterleavedReader(seqFile); break;
	case Phylip:
		ret = PhylipReader(seqFile); break;
	default:
		std::cerr << "\nError: cannot work out format of file " << seqFile;
	}
	if(ret->size() == 0) {
		cerr << "\nError in reading file " << seqFile << "? Couldn't find sequences when looking under format " << FileTypeName(TestFile(seqFile)) << endl;
		exit(-1);
	}
	return ret;
}
// File tester
EFileType TestFile(std::string SeqFile) {
	std::ifstream input(SeqFile.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }
    std::string line;
    std::vector<std::string> Toks;
    while(getline( input, line ) ) {
    	// Looks like FASTA
    	if(line[0] == '>') {
    		while(getline( input, line ) ) {
    			if(!line.empty()) {
    				line = RemoveWhiteSpace(line);
    				if(IsSeq(line[0]) || IsGap(line[0])) { return FASTA; }
    				break;
    	} 	}	}
    	// MSF
    	if(line.find("MULTIPLE_ALIGNMENT") != string::npos || line.find("PileUp") != string::npos) {
    		while(getline( input, line ) ) {
    			if(line[0] == '/' && line[1] == '/') { return MSF; }
    		}
    		cerr << "\nFound something that looks like MSF format, but doesn't appear to have sequences\n"; exit(-1);
    	}
    	// Phylip/Interleaved
    	Toks = Tokenise(line);
    	if(Toks.size() == 2) {
    		int NoSeq = atoi(Toks[0].c_str()), length = atoi(Toks[1].c_str());
    		if(NoSeq > 0 && NoSeq < BIG_NUMBER && length > 0 && length < BIG_NUMBER) { // Confirms it looks like one of these formats
    			while(getline( input, line ) ) {
    				Toks = Tokenise(line);
    				if(Toks.size() > 0) {
    					if(Toks.size() == 1) { return Interleaved; } // If next line contains only a single token then it's interleaved
    					else {
    						bool flag = true;
    						for(int i = 1; i < Toks.size() ; i++) {
    							for(auto &c : Toks[i]) { if(!IsSeq(c) && !IsGap(c)) { flag = false; } }
    						}
    						if(flag) { return Phylip; }
    					}
    					cerr << "\nWeird format on line: " << line;
    					cerr << "\nExpected either a sequence name or a sequence name followed by sequence, but got unexpected characters";
    }	}	}	}	}
    // If not sure, then assume FASTA
    cout << "\nWARNING: Unknown sequence format. Cowardly assuming FASTA...\n";
    return FASTA;
}

// File readers
std::vector <CSequence> *FASTAReader(std::string SeqFile) {
	std::vector <CSequence> *RetSeq = new std::vector<CSequence>();

	std::ifstream input(SeqFile.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }
    std::string line, name, content;
    while(getline( input, line ) ){
    	line = RemoveWhiteSpace(line);
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !line.empty() ){
            	if(!name.empty()) { // Push the sequence back
            		RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
            		name.clear(); content.clear();
            	}
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    // Add the final sequence
    RetSeq->push_back(CSequence(RemoveWhiteSpace(name),RemoveWhiteSpace(content)));
    return RetSeq;
}
std::vector <CSequence> *MSFReader(std::string SeqFile) {
	vector <CSequence> *RetSeq = new std::vector<CSequence>();
	vector <std::string> Names;
	vector<std::string> Toks;
	ifstream input(SeqFile.c_str(), std::ifstream::in);
	string line;
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }
    // Read the file names from the header
    while(getline( input, line ) ){
    	if(line.size() < 2) { continue; }
    	if(line[0] == '/' && line[1] == '/') { break; }
    	Toks = Tokenise(line);
    	if(Toks[0] == "Name:") { Names.push_back(Toks[1]); }
    }
    if(Names.empty()) { std::cerr << "\nCouldn't find names in MSF file header. They should occur after Name:\n"; exit(-1); }
    for(auto &n1 : Names) {
    	int count = 0;
    	for(auto &n2 : Names) {
    		if(n1 == n2) { count ++; }
    	}
    	if(count != 1) { std::cerr << "\nMultiple copies of name " << n1 << "?\n"; exit(-1); }
    }
    std::vector <std::stringstream> Seqstream(Names.size());
    while(getline( input, line ) ){
    	Toks = Tokenise(line);
    	if(Toks.size() < 2) { continue; }
    	// Find the name
    	for(int i = 0 ; i < Names.size(); i++) {
    		if(Toks[0] == Names[i]) {
    			for(int j = 1; j < Toks.size(); j++) { Seqstream[i] << Toks[j]; }
    			break;
    		}
    	}
    }
    for(int i = 0 ; i < Names.size(); i++) {
    	if(Seqstream[i].str().size() != Seqstream[0].str().size()) {
    		std::cerr << "\nERROR: Some sequences seem of different lengths?\n"; exit(-1);
    	}
    	RetSeq->push_back(CSequence(RemoveWhiteSpace(Names[i]),RemoveWhiteSpace(Seqstream[i].str())));
    }
    return RetSeq;
}
std::vector <CSequence> *InterleavedReader(std::string SeqFile) {
	int noSeq = -1, length;
	vector <CSequence> *RetSeq = new std::vector<CSequence>();
	vector <std::string> Names;
	vector<std::string> Toks;
	ifstream input(SeqFile.c_str(), std::ifstream::in);
	string line;
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }
    // Get the line specifying number of sequences and their length
    while(getline( input, line ) ){
    	if(line.size() < 1) { continue; }
    	if(line[0] == '#') { continue; }
    	Toks = Tokenise(line);
    	if(Toks.size() != 2) { cerr << "\nError: first line in a phylip format file should be #seq length"; exit(-1); }
    	noSeq = atoi(Toks[0].c_str()); assert(noSeq > 0);
    	length = atoi(Toks[1].c_str()); assert(length > 0);
    	break;
    }
    if(noSeq < 0) { cerr << "\nCouldn't find number of sequences in file: " << SeqFile << "?\n"; exit(-1); }
    Names.assign(noSeq,"");
    std::vector <std::stringstream> Seqstream(Names.size());
    for(int i = 0; i< Names.size(); i++) {
        while(getline( input, line ) ){
        	if(line.size() < 1) { continue; }
        	if(line[0] == '#') { continue; }
        	Names[i] = line;
        	break;
        }
        while(getline( input, line ) ){
			if(line.size() < 1) { continue; }
			if(line[0] == '#') { continue; }
			Seqstream[i] << RemoveWhiteSpace(line);
			Seqstream[i].seekg(0, ios::end);
			if(Seqstream[i].tellg() >= length) { break; }
		}
    }
    for(int i = 0 ; i < Names.size(); i++) {
     	if(Seqstream[i].str().size() != Seqstream[0].str().size()) {
     		std::cerr << "\nERROR: Some sequences seem of different lengths?\n"; exit(-1);
     	}
     	RetSeq->push_back(CSequence(RemoveWhiteSpace(Names[i]),RemoveWhiteSpace(Seqstream[i].str())));
     }
    return RetSeq;
}
std::vector <CSequence> *PhylipReader(std::string SeqFile) {
	int noSeq = -1, length;
	vector <CSequence> *RetSeq = new std::vector<CSequence>();
	vector <std::string> Names;
	vector<std::string> Toks;
	ifstream input(SeqFile.c_str(), std::ifstream::in);
	string line;
    if(!input.good()){
        std::cerr << "Error opening '"<< SeqFile <<"'. Provide a valid file." << std::endl;
        exit(-1);
    }
    // Get the line specifying number of sequences and their length
    while(getline( input, line ) ){
    	if(line.size() < 1) { continue; }
    	if(line[0] == '#') { continue; }
    	Toks = Tokenise(line);
    	if(Toks.size() != 2) { cerr << "\nError: first line in a phylip format file should be #seq length"; exit(-1); }
    	noSeq = atoi(Toks[0].c_str()); assert(noSeq > 0);
    	length = atoi(Toks[1].c_str()); assert(length > 0);
    	break;
    }
    if(noSeq < 0) { cerr << "\nCouldn't find number of sequences in file: " << SeqFile << "?\n"; exit(-1); }
    Names.assign(noSeq,"");
    std::vector <std::stringstream> Seqstream(Names.size());
    // Go through the rest of the file
    while(getline( input, line ) ){
    	if(line.size() < 1) { continue; }
    	if(line[0] == '#') { continue; }
    	// Hard fixed loop of noSeq with no gaps
    	for(int i = 0; i < noSeq; i++) {
    		Toks = Tokenise(line);
    		if(Toks.empty()) { cerr << "\nError: unexpected line between sequences where " << Names[i] << " is expected"; }
    		if(Names[i].empty()) { Names[i] = Toks[0]; }
    		if(Toks[0] == Names[i]) { Toks.erase(Toks.begin()); }
    		for(auto & s : Toks) { Seqstream[i] << s; }
    		if(i < noSeq - 1) { if(!getline( input, line )) { break; } }
    	}
    }
    for(int i = 0 ; i < Names.size(); i++) {
     	if(Seqstream[i].str().size() != Seqstream[0].str().size()) {
     		std::cerr << "\nERROR: Some sequences seem of different lengths?\n"; exit(-1);
     	}
     	RetSeq->push_back(CSequence(RemoveWhiteSpace(Names[i]),RemoveWhiteSpace(Seqstream[i].str())));
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
	if(Toks.size() > 1) {
		if(Toks[0].size() == 0) { Toks.erase(Toks.begin()); }
	}
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

