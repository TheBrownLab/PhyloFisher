/* Options.cpp
 * **
 * Sets the options from the command line
 *
 *  Created on: Oct 4, 2017
 *      Author: simon
 */

#include "prequal.h"
#include "Sequence.h"

using namespace::std;

string ShowUsage();

string ShowUsage() {
	stringstream t_name;
	for(char c : PROGRAM_NAME) { t_name << (char) tolower(c); }
	return "\n\nUsage: \n\t./" + t_name.str() + " [options] input_file\n\t-h [all] for [full] options\n\nTypical usage (should do a good job with most sequences):\n\t ./" + t_name.str() + " input_file\n\n";

}

COptions::COptions(int argc, char *argv[]) {
	bool setThresh = false, setProp = false;
	cout << "----------------------------------------------------\n\t"
			<< PROGRAM_NAME << " v." << VERSION_NUMBER
			<< "  by Simon Whelan"
			<< "\n----------------------------------------------------";
	// Check basic input is matched
	if (argc < 2) {
		cout << ShowUsage();
		exit(-1);
	}
	// Get the input file
	if (strcmp(argv[1], "-h") != 0) {
		if (argv[argc - 1][0] == '-') {
			cout << "\nInput file must be specified last. you tried input of <" << argv[argc - 1] << ">\n";
			exit(-1);
		}
		_inputFile = argv[argc-1];
		argc--;
	}
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-h") == 0) {
			if(i + 1 < argc) {
				if(strcmp(argv[i+1],"all") == 0) {
					cout << "\n\nOptions affecting the core region and filtering:";
					cout << "\n\t-corerun X       \t: X number of high posterior residues at beginning and end before \n\t\t\t\t\ta core region is defined [DEFAULT "
						 << _runBeforeInside << "]";
					cout << "\n\t-nocore          \t: No core region will be defined [DEFAULT uses core region]";
					cout << "\n\t-removeall       \t: Remove all residues rather than those outside the core region \n\t\t\t\t\t[DEFAULT residues in core filtered with " << CoreFilter() <<"]";
					cout << "\n\t-corefilter X    \t: The character outputted in the core for the  filtered alignment [DEFAULT is "<< CoreFilter() << "]";
					cout << "\n\t-noremoverepeat  \t: Do not remove repeated regions of length >" << _repeatLength << ". Repeats break PPs!";
					cout << "\n\nOptions dealing with protein coding DNA sequences:";
					cout << "\n\t-nodna           \t: Enforces banning of dna sequences. The program will fail if only DNA is input.";
					cout << "\n\t-forceuniversal  \t: Enforces the usage of the universal code. Internal stop codons will be replaced with X for other codes.";
					cout << "\n\nOptions affecting output formats:";
					cout << "\n\t-outsuffix X     \t: Output file will be the original name with X as a suffix [DEFAULT " << OutSuffix() << "]";
					cout << "\n\t-dosummary       \t: Output summary statistics to file suffixed with " << SummarySuffix();
					cout << "\n\t-dodetail        \t: Output detailed statistics to file suffixed with " << DetailSuffix();
					cout << "\n\t-noPP            \t: Stop outputting the posterior probability matrix";
					cout << "\n\nOptions affecting posterior probabilities and filtering:";
					cout << "\n\t-pptype X [Y]       \t: Specify the algorithm used to calculate posterior probabilities\n\t\t\t\t\tX = all : for all against all sequence comparisons";
					cout << "\n\t\t\t\t\tX = closest : for Y closest relatives [DEFAULT; Y = " << PPnumber() << "]\n\t\t\t\t\tX = longest : for comparing the Y longest sequences [Y = " << PPnumber() << "]";
					cout << "\n\t-filterprop X       \t: Filter the sequences so that in total X proportion (range 0.0 - 1.0) \n\t\t\t\t\tof the sequences are maintained. This can be used to tune filtering.";
					cout << "\n\t-filterthresh X     \t: Filter the sequences to the posterior probabilities threshold X [DEFAULT = " << KeepThreshold() << "]\n\t\t\t\t\t(range 0.0 - 1.0). DEFAULT filtering option with threshold";
					cout << "\n\t-filterjoin X       \t: Extend filtering over regions of unfiltered sequence less than X [DEFAULT X = " << FilterRange() << "]";
					cout << "\n\t-nofilterlist X     \t: Specify a file X that contains a list of taxa names that will \n\t\t\t\t\tnot be filtered. In X one name per line.";
					cout << "\n\t-nofilterword X     \t: Specify a file X that contains a list of words and sequence names that contain \n\t\t\t\t\tthose words will not be filtered. In X one word per line.";
					cout << ShowUsage();
					exit(-1);
				}
			}
			cout << "\n\nSimple options (-h all for full options)";
			cout << "\n\t-filterthresh X     \t: Filter the sequences to the posterior probabilities threshold X [DEFAULT = " << KeepThreshold() << "]\n\t\t\t\t\t(range 0.0 - 1.0). DEFAULT filtering option with threshold";
			cout << "\n\t-corerun X       \t: X number of high posterior residues at beginning and end before \n\t\t\t\t\ta core region is defined [DEFAULT "
				 << _runBeforeInside << "]";
			cout << "\n\t-pptype X [Y]       \t: Specify the algorithm used to calculate posterior probabilities\n\t\t\t\t\tX = all : for all against all sequence comparisons";
			cout << "\n\t\t\t\t\tX = closest : for Y closest relatives [DEFAULT; Y = " << PPnumber() << "]\n\t\t\t\t\tX = longest : for comparing the Y longest sequences [Y = " << PPnumber() << "]";
			cout << "\n\t-filterjoin X       \t: Extend filtering over regions of unfiltered sequence less than X [DEFAULT X = " << FilterRange() << "]";
			cout << "\n\t-nofilterlist X     \t: Specify a file X that contains a list of taxa names that will \n\t\t\t\t\tnot be filtered. In X one name per line.";
			cout << ShowUsage();
			exit(-1);
		}
		/////////// Core stuff
		else if (strcmp(argv[i], "-corerun") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -corerun. Specify number afterwards";
				exit(-1);
			}
			_runBeforeInside = atoi(argv[++i]);
			if (!InRange(_runBeforeInside, 0, 100)) {
				cout << "\n-corerun " << _runBeforeInside
						<< " outside good limits (0,100)";
			}
		} else if (strcmp(argv[i], "-nocore") == 0) {
			_remove2Core = false;
		} else if (strcmp(argv[i], "-removeall") == 0) {
			_removeAll = true;
		} else if (strcmp(argv[i], "-corefilter") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -corefilter. Specify a character afterwards";
				exit(-1);
			}
			_coreFilter = argv[++i][0];
		} else if(strcmp(argv[i], "-noremoverepeat") == 0) {
			_removeRepeat = false;
		} else if(strcmp(argv[i], "-nodna") == 0) {
			_allowDNA = false;
		} else if(strcmp(argv[i], "-forceuniversal") == 0) {
			_alwaysUniversal = true;
		}
		//////////// File stuff
		else if (strcmp(argv[i], "-outsuffix") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -outsuffix. Specify a suffix afterwards";
				exit(-1);
			}
			_outputSuffix = argv[++i];
			if (_outputSuffix.size() < 1) {
				cout << "\nMust specify a meaningful outsuffix";
			}
		} else if (strcmp(argv[i], "-dosummary") == 0) {
			_doSummary = true;
		} else if (strcmp(argv[i], "-dodetail") == 0) {
			_doDetail = true;
		} else if (strcmp(argv[i], "-noPP") == 0) {
			_doPPs = false;
		}
		///////////// Filtering and PP stuff
		else if (strcmp(argv[i], "-pptype") == 0) {
			if (i + 1 == argc ) {
				cout << "\nError when definining -pptype. Specify type afterwards (all; closest; longest)";
				exit(-1);
			}
			string pp_option = argv[++i];
			if(pp_option == "all") {
				_ppCalcs = -1;
			} else if(pp_option == "closest") {
				_ppCalcs = 0;
			} else if(pp_option == "longest") {
				_ppCalcs = 1;
			} else {
				cout << "\nUnknown option when defining -pptype : " << pp_option; exit(-1);
			}
			if(i + 1 < argc && _ppCalcs != -1) {
				if(argv[i+1][0] != '-') {
					_ppCalcNumber = atoi(argv[++i]);
					if(!InRange(_ppCalcNumber,1,1000)) {
						cout << "\n-pptype X Y : Y is not in expected range: " << _ppCalcNumber << endl; exit(-1);
					}
				}
			}
		} else if (strcmp(argv[i], "-nofilterlist") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -nofilterlist. Specify a file afterwards";
				exit(-1);
			}
			_noFilterList = GetStringList(argv[++i]);
			cout << "\nFilter list consists of " << _noFilterList.size()
					<< " entries";
		} else if (strcmp(argv[i], "-nofilterword") == 0) {
			if (i + 1 == argc) {
				cout << "\nError when definining -nofilterword. Specify a file afterwards\n";
				exit(-1);
			}
			_noFilterWord = GetStringList(argv[++i]);
			cout << "\nFilter word consists of " << _noFilterWord.size()
					<< " entries";
		} else if (strcmp(argv[i], "-filterthresh") == 0) {
			if (i + 1 < argc) {
				_filterThreshold = atof(argv[++i]);
			} else {
				_filterThreshold = -1;
			}
			if (!InRange(_filterThreshold, 0.0, 1.0)) {
				cout << "\nError when defining -filterthresh. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
			setThresh = true;
		} else if (strcmp(argv[i], "-filterprop") == 0) {
			if (_filterThreshold > 0) {
				_filterThreshold *= -1;	// Turn off the threshold
			}
			if (i + 1 < argc) {
				_keepProportion = atof(argv[++i]);
			} else {
				_keepProportion = -1;
			}
			if (!InRange(_keepProportion, 0.0, 1.0)) {
				cout << "\nError when defining -filterprop. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
			setProp = true;
		} else if (strcmp(argv[i], "-filterjoin") == 0) {
			if (i + 1 < argc) {
				_joinFilterRange = atof(argv[++i]);
			} else {
				_joinFilterRange = -1;
			}
			if(!InRange(_joinFilterRange,0,1000)) {
				cout << "\nNeed to define a value after -filterjoin of 0 or greater \n";exit(-1);
			}
			if (!InRange(_keepProportion, 0.0, 1.0)) {
				cout << "\nError when defining -filterprop. Specify a valid number (0,1) afterwards\n";
				exit(-1);
			}
		} else if(strcmp(argv[i], "-newHMM") == 0) {
			_doHMMpostFilter = true;
		}
		///////////// Default option
		else {
			cout << ShowUsage();
			exit(-1);
		}

	}
	if(setThresh && setProp) {
		cout << "\nERROR: user defined both filter threshold and a filter proportion in commandline?\n\n"; exit(-1);
	}
}

vector<string> COptions::GetStringList(string FilterListFile) {
	ifstream input(FilterListFile.c_str());
	if (!input.good()) {
		std::cerr << "Error opening '" << FilterListFile
				<< "'. Provide a valid file." << std::endl;
		exit(-1);
	}
	string line, name;
	vector < string > FilterList;
	while (getline(input, line)) {
		if (line.size() < 1) {
			continue;
		}
		line = RemoveWhiteSpace(line);
		if(line[0] == '>') { line = line.substr(1); } // Skip the first character if it's a '>' from fasta format
		FilterList.push_back(line);
	}
	return FilterList;
}
