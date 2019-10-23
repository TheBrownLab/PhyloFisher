/*
 * prequal.cpp
 *
 *  Created on: 12 Sep 2017
 *      Author: Simon Whelan
 *	Collaboration and testing: Fabien Burki and Iker Irisarri
 */


#include <algorithm>
#include <iomanip>
#include "prequal.h"
#include "Sequence.h"
#include "ZorroInterface.h"

using namespace::std;

// Global variables. Ugly, but easiest quick solution
vector <CSequence> *data = NULL;
double **PP = NULL;
COptions *options = NULL;			// Global options; would be better as a singleton, but this works
stringstream warningStream;

// Some local functions
void DoHMMfiltering();
void GetHMMLabelling(double er, double leave, CSequence &seq, int seqNum);

int main(int argc, char * argv[]) {

	// Collect options
	options = new COptions(argc, argv);
	CSequence::SetFilter(options->CoreFilter());
	double threshold;
	vector<double> values;							// Generic vector for moving around values

	// Read data and sort initialisation
	data = FASTAReader(options->Infile(),options->AlwaysUniversal()); // Reads the sequences
	cout << "\nThere are " << data->size() << " sequences of max length " << CSequence::MaxLength();
	// Some error checking with regard to DNA sequences
	for(auto & seq : *data) {
		if(seq.HasDNA() && !options->AllowDNA()) {
			cout << "\nHave a DNA sequence while the -nodna option is active. Exiting...\n\n"; exit(-1);
		}
		if(seq.HasDNA() != data->at(0).HasDNA()) {
			cout << "\nError: potential mix of DNA and amino acid sequences. This is not handled?\n\n"; exit(-1);
		}
	}
//	for(int i = 0; i < data->size(); i++) { cout << "\ni=" << i << "\t" << data->at(i).out(); }

	// Do the repeat filtering if required
	bool repeated = false;
	if(options->RemoveRepeat()) {
		for(int i = 0; i < data->size(); i++) {
			if(options->IgnoreSequence(data->at(i).Name())) { continue; }
			int oriLength = data->at(i).length();
			if(data->at(i).CleanRepeat(options->RepeatLength())) {
				repeated = true;
				warningStream << "\nWARNING: repeat found in sequence["<<i<<"] " << data->at(i).Name();
				warningStream << "\n\tLength changed from " << oriLength << " -> " << data->at(i).length();
			}
		}
	}
	if(repeated) {
		CSequence::ResetMaxLength(data);
		warningStream << "\nWARNING: Repeats found. If you do not want them removed us -noremoverepeat option";
	}

	// Run the HMM if needed
	PP = RunHMM(data,options->Infile() + options->OutSuffix() + options->PPSuffix(), options->Overwrite_PP());

	// Define the threshold
	if (options->DoKeepProportion()) {
		cout << "\n\nExamining posterior probabilities to determine appropriate thresholds to retain " << options->KeepProportion() * 100 << "% of sequence" << flush;
		threshold = TargetCutoff(options->KeepProportion());
	} else {
		threshold = options->KeepThreshold();
		if(threshold != DEFAULT_THRESHOLD) { cout << "\n\nThreshold set to " << threshold << flush; }
	}
	assert(InRange(threshold, 0.0, 1.0));

	// Do the filtering
	DoFiltering(threshold);

	////////////////////////////////////////////////////////////////
	cout << "\n\nOutputting results: ";
	if (options->DoDetail()) {
		cout << "\n\tDoing detailed output to " << options->Infile() << options->DetailSuffix() << flush;
		ofstream detail_out(options->Infile() + options->DetailSuffix());
		detail_out << std::fixed;
		detail_out << std::setprecision(4);
		detail_out << "# [seq_pos]seq_character\tmaxPP\tToRemove\tInside\n";
		for (int i = 0; i < data->size(); i++) {
			values.clear();
			double max = 0.0;
			for(int j = 0; j < data->at(i).length(); j++) {
				if(PP[i][j] > max) { max = PP[i][j]; }
				values.push_back(PP[i][j]);
			}
			detail_out << ">" << data->at(i).Name();
//			detail_out << "\nmean= " << mean(values) << " : stdev= " << stdev(values) << " : poscut= " << mean(values) - (4*stdev(values)) << " : max= " << max;
			for(int j = 0; j < data->at(i).length(); j++) {
				detail_out << "\n["<< j<< "]" << data->at(i).RealSeq(j) << "\t" << PP[i][j] << "\t" << data->at(i).Remove[j] << "\t" << data->at(i).Inside[j];
			}
			detail_out << endl;
		}
		detail_out.close();
		cout << " ... done" << flush;
	}
	if(options->DoSummary()) {
		cout << "\n\tDoing summary output to " << options->Infile() << options->SummarySuffix() << flush;
	}
	cout << "\n\tOutputting filtered amino acid sequences to " << options->Infile() << options->OutSuffix();
	int total_char = 0;
	int output_char = 0;
	int output_seq = 0;
	bool hasDNA = false;
	ofstream sequence_out(options->Infile() + options->OutSuffix());
	for(int i = 0; i < data->size(); i++) {
		if(data->at(i).HasDNA()) { hasDNA = true; }
		total_char += data->at(i).length();
		if(data->at(i).AllRemoved()) {
			warningStream << "\nWARNING: Fully removed sequence [ "<<i<<"] " << data->at(i).Name();
			continue;
		}
		if(data->at(i).PropRemoved > 0.25) {
			warningStream << "\nWARNING: " << fixed << setprecision(2) <<  data->at(i).PropRemoved * 100 << "% of sequence removed/missing for ["<<i<<"] " << data->at(i).Name();
		}
		output_seq++;
		sequence_out << ">" << data->at(i).Name() << endl;
		// The filtered style sequence determined by the CSequence class with some added counting stuff
		string output = data->at(i).Seq();
		for(int j = 0; j < output.size(); j++) { if(output[j] != options->CoreFilter()) { output_char ++; } }
		sequence_out << output << endl;
	}
	// If DNA sequences output these and also some translation information
	if(hasDNA) {
		string outDNAfile = options->Infile() + ".dna" + options->OutSuffix();
		cout << "\n\tOutputting filtered DNA sequences to " << outDNAfile;
		ofstream outDNA(outDNAfile);
		for(auto &seq : *data) {
			if(seq.AllRemoved()) { continue; }
			outDNA << ">" << seq.Name();
			outDNA << endl << seq.DNA() << endl;
		}
		outDNA.close();
		outDNAfile = options->Infile() + ".translation";
		cout << "\n\tOutputting DNA translation information to " << outDNAfile;
		ofstream outTrans(outDNAfile);
		for(auto &seq : *data) {
			outTrans << ">" << seq.Name() << " translated using " << seq.GenCode();
			outTrans << "\n" << seq.RealDNA();
			outTrans << "\n" << seq.RealSeq() << "\n";
		}
		outTrans << "\n\nNOTE: If you are unhappy with some of these translations it is probably because a DNA sequence contains a stop codon"
				<< "\n  in your preferred genetic code. Try editing the file before you filter. You'll have to delete the .PP file to run again";
		outTrans.close();

	}

	cout << "\n\nComputation complete";
	// Make nice summary of information
	cout << "\n\n===================== Summary ======================";
	cout << "\n              " << std::setw(12) << "Original" << std::setw(12) << "Filtered" << std::setw(12) << "%Retained";
	cout << "\n#Sequences    " << std::setw(12) << data->size() << std::setw(12) << output_seq << std::setw(12) << std::setprecision(3) << 100 * (double)output_seq/(double)data->size() << "%";
	cout << "\n#Residues     " << std::setw(12) << total_char << std::setw(12) << output_char << std::setw(12) << std::setprecision(3) << 100 * (double)output_char/(double)total_char << "%";
	cout << "\n\n";
	sequence_out.close();

	if(options->DoSummary()) {
		ofstream summary_out(options->Infile() + options->SummarySuffix());
		summary_out << std::fixed;
		summary_out << std::setprecision(4);
		// Calculate statistics
		double rem_mean = 0, rem_max = 0, in_mean = 0, in_min = 1.0;
		int rem_index = -1, in_index = -1;
		for(int i = 0; i < data->size(); i++) {
			// Removed
			rem_mean += data->at(i).PropRemoved;
			if(data->at(i).PropRemoved > rem_max) { rem_max = data->at(i).PropRemoved; rem_index = i; }
			// Inside
			in_mean += data->at(i).PropInside;
			if(data->at(i).PropInside < in_min) { in_min = data->at(i).PropInside; in_index = i; }
		}
		rem_mean /= (double) data->size();
		in_mean /= (double) data->size();
		// Output
		summary_out << "\nThere are " << data->size() << " sequences";
		summary_out << "\nRemoval:\n\tOn average " << rem_mean * 100 << "% of sequence retained";
		if(rem_index >= 0) { summary_out << "\n\tSequence with most removed (" << rem_max * 100 << "%) is [" << rem_index << "] = "<< data->at(rem_index).Name(); }
		summary_out << "\nCore regions:\n\tOne average " << in_mean * 100 << "% of sequence is in the core region";
		if(in_index >= 0) { summary_out << "\n\tSequence with least in core (" << in_min * 100 << "%) is [" << in_index << "] = "<< data->at(in_index).Name(); }
		summary_out << "\n##";
		for(int i = 0; i < data->size(); i++) {
			summary_out << "\n["<<i<<"] " << data->at(i).Name() << " has " << data->at(i).PropRemoved * 100 << "% removed and " << data->at(i).PropInside* 100 << "% in the core";
		}
		// Make a nice table of the RHS of the PP distribution
		TargetCutoff(0.1,summary_out);
		// Make nice summary of information
		summary_out << "\n\n===================== Summary ======================";
		summary_out << "\n              " << std::setw(12) << "Original" << std::setw(12) << "Filtered" << std::setw(12) << "%Retained";
		summary_out << "\n#Sequences    " << std::setw(12) << data->size() << std::setw(12) << output_seq << std::setw(12) << std::setprecision(3) << 100 * (double)output_seq/(double)data->size() << "%";
		summary_out << "\n#Residues     " << std::setw(12) << total_char << std::setw(12) << output_char << std::setw(12) << std::setprecision(3) << 100 * (double)output_char/(double)total_char << "%";
		summary_out << "\n\n";
		summary_out.close();
	}
	// Output warnings if required
	if(warningStream.rdbuf()->in_avail() != 0) {
		string warningFile = options->Infile() + ".warning";
		cout << "Analysis may have some problems. Warnings output to " << warningFile << "\n\n";
		ofstream warnOut(warningFile);
		warnOut << warningStream.str();
		warnOut.close();
	}
	cout << "Filtering complete!\n";
	// Clean up memory
	for(int i = 0; i < data->size(); i++) { delete [] PP[i];  } delete [] PP;
	delete data;
}


// Returns the cutoff based on the empirical set of PPs in PP[][]
double TargetCutoff(double prop2Keep, ostream &os) {
	int total_length = 0;
	vector <double> tmp_PP;
	for(int i = 0; i < data->size(); i++)  {
		total_length += data->at(i).length();
		for(int j = 0; j < data->at(i).length(); j++)  {
			tmp_PP.push_back(PP[i][j]);
		}
	}
	std::sort(tmp_PP.begin(),tmp_PP.end());
	int count_stop;
	os << "\n\nHelpful cut-offs ([PropRetained] Cutoffs):";
	os << std::fixed;
	os << std::setprecision(4);

	int spacer = 5;
	for(double x = 1.0; x >= 0.75; x-= 0.01,spacer ++) {
		if(spacer >= 5) { os << "\n"; spacer = 0; }
		count_stop = (int)((1.0 - x) * (double) total_length);
		os << "\t[" << x << "] " << tmp_PP[count_stop];
	}
	count_stop = (int)((1.0 - prop2Keep) * (double) total_length);
	return tmp_PP[count_stop];
}

void DoFiltering(double threshold) {
	cout << "\n\nPerforming filtering:";
	cout << "\n\tApplying standard threshold " << threshold;
	int thresholdCount = 0;
	int ignoreCount = 0;
	// If required rework the posterior probabilities according to a length distribution HMM
	if(options->DoHMMpostFilter()) { DoHMMfiltering(); }

	// Apply the threshold in a simple way
	for (int i = 0; i < data->size(); i++) {
		// 1. Skip sequences on the filter lists
		if(options->IgnoreSequence(data->at(i).Name())) { ignoreCount++; continue; }
		// 2. Find the residues to be filtered
		//		cout << " filtered residues ..." << flush;
		for (int j = 0; j < data->at(i).length(); j++) {
			if (PP[i][j] < threshold) {
				thresholdCount++;
				data->at(i).Remove[j] = true;
			}
		}
	}
	cout << " resulting in " << thresholdCount << " residues removed" << flush;
	// Do the joining of filtered/outside regions if options require so
	if(options->FilterRange() > 0) {
		cout << "\n\tExtending filtered regions with width of " << options->FilterRange() << " ";
		int filter_count = 0;
		for(int i = 0 ; i < data->size(); i++) {
			if(options->IgnoreSequence(data->at(i).Name())) { continue; }
			int lastFilter = 0;
			for(int j = 0; j < data->at(i).length(); j++) {
				if(data->at(i).Filter(j)) {
					if(j - lastFilter < options->FilterRange() && j - lastFilter > 1) {
//						cout << "\nRemoving seq["<<i<<"]["<<j<<"] " << data->at(i).Name() << "== >"<< data->at(i).Seq(j) << "< : due to value " << j - lastFilter << " range (" << lastFilter << "," << j << ")";
//						cout << "\n\t" << data->at(i).RealSeq().substr(lastFilter,j-lastFilter) << " : length " << j - lastFilter;
						for(int k = j; k > lastFilter; k--) {
							data->at(i).Remove[k] = true;
						}
						filter_count++;
					}
					lastFilter = j;
				}
			}
		}
		cout << " ... " << filter_count << " additional regions removed" << flush;
	}
	// Tidy the front and back
	if(options->RunBeforeInside() > 0 && options->Remove2Core()) {
		cout << "\n\tApplying front/back trimming for runs of " << options->RunBeforeInside();
		int seqTrimmed = 0;
		for (int i = 0; i < data->size(); i++) {
			if(options->IgnoreSequence(data->at(i).Name())) { continue; }
			// 1. Get the front ...
			bool DoOutside = false;
			for(int j = my_min(data->at(i).length(),options->RunBeforeInside())-1; j > 0; j--) {
				if(data->at(i).Filter(j)) { if(!DoOutside) { seqTrimmed++;} DoOutside = true;}
				if(DoOutside) { data->at(i).Inside[j] = false; data->at(i).Remove[j] = true; }
			}
			for(int j = 0; j < data->at(i).length(); j++) {
				if(!data->at(i).Filter(j)) { break; }
				data->at(i).Inside[j] = false;
			}
			// 2. And the back ...
			DoOutside = false;
			for(int j = my_max(0,data->at(i).length() - options->RunBeforeInside()); j < data->at(i).length(); j++)  {
				if(data->at(i).Filter(j)) { if(!DoOutside) { seqTrimmed++;} DoOutside = true; }
				if(DoOutside) { data->at(i).Inside[j] = false; data->at(i).Remove[j] = true; }
			}
			for(int j = data->at(i).length() - 1; j >=0; j--) {
				if(!data->at(i).Filter(j)) { break; }
				data->at(i).Inside[j] = false;
			}
		}
		cout << " resulting in " <<seqTrimmed << " sections removed" << flush;
	}
	// Get the summary statistics
	for(int i = 0; i < data->size(); i++) {
		data->at(i).CalculateSummary();
	}
	if(ignoreCount > 0) { cout << "\n\tThere were " << ignoreCount << " sequences ignored by filtering due to word/name lists"; }
}

// Extends the edges if required
void DoHMMfiltering() {

	cout << "\nThis is some experimental annotation and should not be used by anyone but SW! Go away!"; exit(-1);

	// 2-state HMM
	// Good : Bad
	// G:G = 1 - er [er = ErrorRate]
	// G:B = er
	// B:G = exit = (1 / l) [l = expected error length; geometric]
	// BB = 1 - (1 / l)
	for(int i = 0; i < data->size() ; i++) {
		GetHMMLabelling(options->ExpErrorProb(), 1.0 / (double) options->ExpErrorLength(), data->at(i) ,i);
	}
}


void GetHMMLabelling(double er, double leave, CSequence &seq, int seqNum) {
	int end = seq.length()-1;
	// Storage
	vector<double> fGood(seq.length(),0.0), fBad(seq.length(),0.0);	// Forward
	vector<double> bGood(seq.length(),0.0), bBad(seq.length(),0.0);	// Backward
	double logF, logB;
	// Eqm
	double fG = log(leave/(er+leave));	// pi_Good
	double fB = log(er/(er+leave));	// pi_Bad
	// transitions
	double GG = log(1.0 - er), GB = log(er), BG = log(leave), BB = log(1 - leave);

	// Rework the PPs so no weird edge cases
	double fix = 0.0001;
	for(int i = 0; i < seq.length(); i++) { if(PP[seqNum][i] < fix) { PP[seqNum][i] = fix; } if(PP[seqNum][i] > 1.0 - fix) { PP[seqNum][i] = 1.0 - fix; } }

	// Forward
	//////////////////////////////////////////
	// Initialise
	fGood[0] = fG + GG + log(PP[seqNum][0]);
	addLogProb(fB + BG + log(PP[seqNum][0]), fGood[0]);
	fBad[0] = fG + GB + log(1.0 - PP[seqNum][0]);
	addLogProb(fB + BB + log(1.0 - PP[seqNum][0]), fBad[0]);
	// Recurse
	for(int pos = 1 ; pos < seq.length(); pos++) {
		// Good
		fGood[pos] = fGood[pos-1] + GG + log(PP[seqNum][pos]);
		addLogProb(fBad[pos-1] + BG + log(PP[seqNum][pos]), fGood[pos]);
		// Bad
		fBad[pos] = fGood[pos-1] + GB + log(1.0 - PP[seqNum][pos]);
		addLogProb(fBad[pos-1] + BB + log(1.0 - PP[seqNum][pos]), fBad[pos]);
	}
	logF = fGood[end] + GG + fG;
	addLogProb(fGood[end] + GB + fB, logF);
	addLogProb(fBad[end] + BG + fG, logF);
	addLogProb(fBad[end] + BB + fB, logF);
	// **
	// Backward
	//////////////////////////////////////////
	// Initialise
/*	bGood[end] = fG + GG + log(PP[seqNum][end]);
	addLogProb(fB + BG + log(PP[seqNum][end]), bGood[end]);
	bBad[end] = fG + GB + log(1.0 - PP[seqNum][end]);
	addLogProb(fB + BB + log(1.0 - PP[seqNum][end]), bBad[end]);
	*/
	/*
	fGood[end] = fG + GG + log(PP[seqNum][end]);
	addLogProb(fB + BG + log(PP[seqNum][end]), fGood[end]);
	fBad[end] = fG + GB + log(1.0 - PP[seqNum][end]);
	addLogProb(fB + BB + log(1.0 - PP[seqNum][end]), fBad[end]);
	*/
	fGood[end] = fG + GG + log(PP[seqNum][end]);
	addLogProb(fB + BG + log(PP[seqNum][end]), fGood[end]);
	fBad[end] = fG + GB + log(1.0 - PP[seqNum][end]);
	addLogProb(fB + BB + log(1.0 - PP[seqNum][end]), fBad[end]);

	// Recurse
	for(int pos = seq.length()-2 ; pos >= 0 ; pos--) {
		// Good
		bGood[pos] = bGood[pos+1] + GG + log(PP[seqNum][pos]);
		addLogProb(bBad[pos+1] + BG + log(PP[seqNum][pos]), bGood[pos]);
		// Bad
		bBad[pos] = bGood[pos+1] + GB + log(1.0 - PP[seqNum][pos]);
		addLogProb(bBad[pos+1] + BB + log(1.0 - PP[seqNum][pos]), bBad[pos]);
	}
	logB = bGood[0] + GG + fG;
	addLogProb(bGood[0] + GB + fB, logB);
	addLogProb(bBad[0] + BG + fG,logB);
	addLogProb(bBad[0] + BB + fB,logB);

	cout << "\nHMM forward/backward not working quite right..."; exit(-1);
	cout << "\nFor sequence " << seqNum << " = " << seq.Name() << "\n\tForward: " << logF << "\tBackward: " << logB;
	for(int pos = 0; pos < seq.length() ; pos++) {
		cout << "\n\t["<<pos<<"] " << fGood[pos] << "\t" << fBad[pos] << "\t" << bGood[pos] << "\t" << bBad[pos]  << " from PP: " << PP[seqNum][pos] << " -> newPP(good): " << (fGood[pos] + bGood[pos] - logF) << " ->" << exp(fGood[pos] + bGood[pos] - logF);
	}

	// Write over the original PPs with these new ones
	for(int pos = 0; pos < seq.length() ; pos++) {
		PP[seqNum][pos] = exp(fGood[pos] + bGood[pos] - logF);
	}

	exit(-1);
}


double mean(vector <double> vec) {
	double ret = 0.0;
	for(int i = 0; i < vec.size(); i++) {
		ret += vec[i];
	}
	return ret / (double)vec.size();
}

double stdev(vector <double> vec) {
	double ret = 0, ave = mean(vec);
	for(int i = 0; i < vec.size(); i++) {
		ret += (vec[i] - ave) * (vec[i] - ave);
	}
	ret /= (double) vec.size();
	return sqrt(ret);
}
