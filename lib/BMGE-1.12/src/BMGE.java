/*
  BMGE (Block Mapping and Gathering with Entropy): selection of phylogenetic 
  informative regions from multiple sequence alignments.
  [Version 1.12]
  Copyright (C) 2010  Alexis Criscuolo 

  This file is part of BMGE.

  BMGE is free software;  you can redistribute it and/or modify it under the 
  terms of the GNU General Public License  as published by the Free Software 
  Foundation; either version 2 of the License, or (at your option) any later 
  version.
  
  BMGE is distributed  in the hope that  it will be useful,  but WITHOUT ANY 
  WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS 
  FOR A  PARTICULAR PURPOSE.  See the  GNU General  Public License  for more 
  details.
  
  You should have  received a copy  of the GNU  General Public License along 
  with this program; if not, write to the 
  Free Software Foundation, Inc., 
  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Contact: 
  Unite de Biologie Moleculaire du Gene chez les Extremophiles
  Departement de Microbiologie
  INSTITUT  PASTEUR
  25 rue du Dr Roux - 75015 Paris  (France)

  alexis.criscuolo@pasteur.fr
*/

import jap.*;
import Jama.*;
import java.io.*;
import java.util.*;
import java.text.*;

public class BMGE {

    //##### constants #####
    static final byte AA = 0; static final byte DNA = 1; static final byte COD = 2;  static final byte RY = 3; 
    static final byte PHYLIP = 0; static final byte FASTA = 1; static final byte PAUP = 2; static final byte HTML = 3;
    static final byte PHYLIP_TAX = 4; static final byte PAUP_TAX = 5; static final byte PHYLIP_ACCN = 6; static final byte PAUP_ACCN = 7; 
    static final byte NO = 0; static final byte WARNING = 1; static final byte YES = 2; static final byte FAST = 3;
    static final String BLANK = "                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ";
    static final int REPLICATE = 10;
    static final double REJECTION_PVALUE = 0.1;
    static final char[] NT = {'A','C','G','T'};
    static final char[] AM = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

    static final int FAST_STATIONARY_DIV = 1000;
    // Khi2 5% critical values
    static final double KHI2DOF3 = 7.815;
    static final double KHI2DOF19 = 30.144;
    // Khi2 10% critical values
    //static final double KHI2DOF3 = 6.251;
    //static final double KHI2DOF19 = 27.203;
    

    //##### options #####
    static final String[] optName = {"-i","-t","-m","-w","-g","-o","-c","-h","-b","-l","-s"};
    static String[] optChoice;
    static byte characterState;           // -t
    static double hmin, hmax;             // -h
    static Matrix similarity;             // -m
    static double rowGapRate, colGapRate; // -g
    static int slidingWindow;             // -w
    static int minBlockSize;              // -b
    static int lengthDisplay;             // -l
    static byte stationarity;             // -s

    //##### data storage #####
    static Alignment alignment, alignment_, alignmentCleaned, alignmentComplement;
    static BitSet keepRow, keepCol;
    static ArrayList<String> rwLabel;

    //##### entropy and gap rate computing #####
    static ArrayList<Double> entropy, smoothedEntropy, gapCol, gapRow;
    static double[] frequency, eigenvalue;
    static Matrix density;

    //##### stationarity estimation #####
    static ArrayList<ArrayList<Meter>> alignmentMeter, alignmentGamma; 
    static ArrayList<ArrayList<Double>> alignmentStuartValue;
    static double stuartMax, sigma, sigmaMax, sigmaCentile;
    static int b_sigmaMax;
    
    //##### i/o #####
    static File infile, outfile;          // -i -o -c
    static BufferedReader in;
    static BufferedWriter out;
    static ArrayList<String> outputOption;
    static ArrayList<File> outputFileName;
    static byte outputFormat, outputType;
    static int codonCode, maxLabelLength;

    //##### stuffs #####
    static int as, b, bb, i, j, c, n, cpt1, cpt2, left, right, lgth, s1, s2, r, best_character, min_stscore, max, inv;
    static double d, g, h, v, up, down, maxGapFrequency, pValue, score;
    static char ch, c1, c2;
    static String line, label, codon1, codon2;
    static StringBuffer sequence;
    static ArrayList<Integer> blockStart, blockEnd;
    static boolean ok, previousCol, recompute, pam0;
    static BitSet removed;
    static Meter meter;
    static ArrayList<Double> sigmaArray, sigmaArray2;

    //##### static method stuffs #####
    static int i_, l_;
    static double[] array_;
    static String _b;


    public static void main(String[] args) throws IOException {

	//################################
	//#####   reading options   ######
	//################################
	Arrays.sort(optName);
	optChoice = new String[optName.length];
	Arrays.fill(optChoice , "null");
	outputOption = new ArrayList<String>(0);
	outputFileName = new ArrayList<File>(0);
	if ( args.length == 0 ) {
	    System.out.println("   mandatory parameters: -i 'infile' -t 'type'");
	    System.out.println("   use option -? for a description of the arguments");
	    System.exit(0);
	}
	    
	try {
	    b = -1;
	    while ( ++b < args.length ) {
		if ( args[b].equals("-i")        // input file
		     || args[b].equals("-h")     // gap rate
		     || args[b].equals("-g")     // gap rate
		     || args[b].equals("-w")     // sliding window
		     || args[b].equals("-b")     // min block size
		     || args[b].equals("-t")     // character state type
		     || args[b].equals("-l")     // length display
		     || args[b].equals("-s")     // stationarity
		     || args[b].equals("-m") ) { // similarity matrix
		    line = args[b];
		    optChoice[ Arrays.binarySearch(optName , line) ] = args[++b];
		}
		else {
		    if ( args[b].equals("-?") ) {
			displayUserGuide(); 
			System.exit(0);
		    }
		    if ( args[b].startsWith("-o") || args[b].startsWith("-c") ) {
			outputOption.add(args[b]); outputFileName.add(new File(args[++b]));
		    }
		}
		if ( args[b].startsWith("-") ) {
		    System.out.println("   incorrect option " + args[b-1]);
		    System.exit(0);
		}
	    }
	} catch ( ArrayIndexOutOfBoundsException e ) {
	    System.out.println("   incorrect options");
	    System.exit(0);
	}
	//##### input file #####
	b = Arrays.binarySearch(optName , "-i");
	if ( optChoice[b].equals("null") ) {
	    System.out.println("   the input file name is not given (option -i)");
	    System.exit(0);
	}
	infile = new File(optChoice[b]);
	if ( ! infile.exists() ) {
	    System.out.println("   problem with the input file : " + infile.toString() + " does not exist");
	    System.exit(0);
	}
	//##### character state type #####
	b = Arrays.binarySearch(optName , "-t");
	if ( optChoice[b].equals("null") ) {
	    System.out.println("   the type of character states is not given (option -t [AA,DNA,RNA,CODON])");
	    System.exit(0);
	}
	characterState = -1;
	if ( optChoice[b].equals("AA") ) characterState = AA;
	if ( optChoice[b].equals("DNA") || line.equals("RNA") ) characterState = DNA;
	if ( optChoice[b].equals("CODON") ) characterState = COD;
	if ( characterState == -1 ) {
	    System.out.println("   invalid type of character states (option -t [AA,DNA,RNA,CODON])");
	    System.exit(0);
	}
	//##### sliding windows #####
	b = Arrays.binarySearch(optName , "-w");
	if ( optChoice[b].equals("null") ) slidingWindow = 3;
	else {
	    try { slidingWindow = Integer.parseInt(optChoice[b]); }
	    catch ( NumberFormatException e ) {
		System.out.println("   incorrect sliding window size (option -w)");
		System.exit(0);
	    }
	}
	if ( slidingWindow % 2 == 0 ) {
	    System.out.println("   the size of the sliding window must be odd (option -w)");
	    System.exit(0);
	}
	//##### minimum block size #####
	b = Arrays.binarySearch(optName , "-b");
	if ( optChoice[b].equals("null") ) minBlockSize = 5;
	else {
	    try { minBlockSize = Integer.parseInt(optChoice[b]); }
	    catch ( NumberFormatException e ) {
		System.out.println("   incorrect minimum block size (option -b)");
		System.exit(0);
	    }
	}
	if ( minBlockSize < 0 ) {
	    System.out.println("   the minimum block size must be greater than 0 (option -b)");
	    System.exit(0);
	}
	//##### gap rate #####
	b = Arrays.binarySearch(optName , "-g");
	rowGapRate = 1.0; colGapRate = 0.2;
	if ( ! optChoice[b].equals("null") ) {
	    line = optChoice[b];
	    c = line.indexOf(":");
	    if ( c != -1 ) {
		try {
		    rowGapRate = Double.parseDouble(line.substring(0 , c));
		    colGapRate = Double.parseDouble(line.substring(++c));
		} catch ( NumberFormatException e ) {
		    System.out.println("   the gap rates are not correct (option -g)");
		    System.exit(0);
		}
	    }
	    else {
		try {
		    colGapRate = Double.parseDouble(line);
		} catch ( NumberFormatException e ) {
		    System.out.println("   the gap rate is not correct (option -g)");
		    System.exit(0);
		}
	    }
	    if ( (Math.min(rowGapRate,colGapRate) < 0) || (Math.max(rowGapRate,colGapRate) > 1.0) ) {
		System.out.println("   the gap rates must range from 0 to 1");
		System.exit(0);
	    }
	}  
	//##### entropy cut-off #####
	b = Arrays.binarySearch(optName , "-h");
	hmin = 0; hmax = 0.5;
	if ( ! optChoice[b].equals("null") ) {
	    line = optChoice[b];
	    c = line.indexOf(":");
	    if ( c != -1 ) {
		try {
		    hmin = Double.parseDouble(line.substring(0 , c));
		    hmax = Double.parseDouble(line.substring(++c));
		} catch ( NumberFormatException e ) {
		    System.out.println("   the entropy cut-off values are not correct (option -h)");
		    System.exit(0);
		}
	    }
	    else {
		try {
		    hmax = Double.parseDouble(line);
		} catch ( NumberFormatException e ) {
		    System.out.println("   the entropy cut-off value is not correct (option -h)");
		    System.exit(0);
		}
	    }
	    if ( (Math.min(hmin,hmax) < 0) || (Math.max(hmin,hmax) > 1.0) ) {
		System.out.println("   the entropy cut-off values must range from 0 to 1");
		System.exit(0);
	    }
	}  
	//##### similarity matrix #####
	similarity = Model.getID(1); pam0 = false;
	b = Arrays.binarySearch(optName , "-m");
	if ( optChoice[b].equals("null") ) {
	    switch ( characterState ) {
	    case COD:
	    case AA: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM62); break;
	    case DNA: similarity = Model.getDNAPAM(100 , 2);
	    }
	}
	else {
	    if ( optChoice[b].startsWith("BLOSUM") ) {
		if ( (characterState != AA) && (characterState != COD) ) {
		    System.out.println("   the similarity matrix (option -m) is not compatible with the character states (option -t)");
		    System.exit(0);
		}
		try {c = Integer.parseInt(optChoice[b].substring(6));} catch ( NumberFormatException e ) {c = 0;}
		switch (c) {
		case 30: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM30); break;
		case 35: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM35); break;
		case 40: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM40); break;
		case 45: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM45); break;
		case 50: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM50); break;
		case 55: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM55); break;
		case 60: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM60); break;
		case 62: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM62); break;
		case 65: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM65); break;
		case 70: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM70); break;
		case 75: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM75); break;
		case 80: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM80); break;
		case 85: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM85); break;
		case 90: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM90); break;
		case 95: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM95); break;
		case 100: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM100); break;
		default: similarity = Model.getAminoAcidSimilarity(Model.BLOSUM62); break;
		}
	    }
	    if ( optChoice[b].equals("PAM0") || optChoice[b].equals("ID") ) {
		switch ( characterState ) {
		case COD:
		case AA: similarity = Model.getID(20); break;
		case DNA: similarity = Model.getID(4); break;
		} 
		pam0 = true;
	    }
	    if ( optChoice[b].startsWith("DNAPAM") ) {
		if ( characterState != DNA ) {
		    System.out.println("   the similarity matrix (option -m) is not compatible with the character states (option -t)");
		    System.exit(0);
		}
		line = optChoice[b].substring(6);
		c = line.indexOf(":");
		try {
		    if ( c == -1 ) i = Integer.parseInt(line);
		    else i = Integer.parseInt(line.substring(0 , c));
		} catch ( NumberFormatException e ) {
		    System.out.println("   the similarity matrix DNAPAM must be defined with an integer value, e.g. DNAPAM47 (option -m)");
		    System.exit(0);
		}
		try {
		    if ( c == -1 ) h = 1.0;
		    else h = Double.parseDouble(line.substring(++c));
		} catch ( NumberFormatException e ) {
		    System.out.println("   incorrect DNAPAM transition/tranversion ratio (option -m)");
		    System.exit(0);
		}
		similarity = Model.getDNAPAM(i , h);
	    }
	    if ( similarity.getRowDimension() * similarity.getColumnDimension() == 1 ) {
		System.out.println("   incorrect similarity matrix name (option -m)");
		System.exit(0);
	    }
	}
	//##### stationarity-based cleaning #####
	b = Arrays.binarySearch(optName , "-s");
	stationarity = NO;
	if ( optChoice[b].toUpperCase().equals("NO") )  stationarity = NO;
	if ( optChoice[b].toUpperCase().equals("YES") ) stationarity = YES;
	if ( optChoice[b].toUpperCase().equals("FAST") ) stationarity = FAST;
	//##### length display #####
	b = Arrays.binarySearch(optName , "-l");
	if ( optChoice[b].equals("null") ) lengthDisplay = Integer.MAX_VALUE;
	else {
	    try { lengthDisplay = Integer.parseInt(optChoice[b]); }
	    catch ( NumberFormatException e ) {
		System.out.println("   incorrect number of characters per line (option -l)");
		System.exit(0);
	    }
	}
	if ( lengthDisplay < 0 ) {
	    System.out.println("   the number of characters per line must be greater than 0 (option -l)");
	    System.exit(0);
	}


	if ( (stationarity != NO) && (characterState == COD) ) {
	    System.out.println("  the stationary-based character trimming (option -s) is not available for codon alignments (-t CODON)");
	    System.out.print("  set character states as nucleotides (-t DNA -s ");
	    if ( stationarity == YES ) System.out.print("YES"); else System.out.print("FAST"); 
	    System.out.println(") or convert them into amino acids (-t CODON -h 1 -g 1 -w 1 -oaa)");
	    System.out.println("  in order to use the stationary-based character trimming (option -s)");
	    System.exit(0);
	}



	//##################################
	//#####   reading alignment   ######
	//##################################
	in = new BufferedReader(new FileReader(infile));
	alignment_ = new Alignment();
	try {
	    label = "phylip";
	    line = "";
	    while ( line.length() == 0 ) line = in.readLine().trim();
	    if ( line.startsWith(">") ) {  // fasta format
		label = line.substring(1).trim();
		sequence = new StringBuffer("");
		line = in.readLine().trim();
		while ( ! line.startsWith(">") ) {
		    sequence = sequence.append(line);
		    line = in.readLine().trim();
		}
		alignment_.add(sequence.toString() , label);
		while ( true ) {
		    label = line.substring(1);
		    sequence = new StringBuffer("");
		    line = in.readLine().trim();
		    while ( ! line.startsWith(">") ) {
			sequence = sequence.append(line);
			line = in.readLine().trim();
		    }
		    alignment_.add(sequence.toString() , label);
		}
	    }
	    else {
		while ( true ) {
		    line = in.readLine().trim();
		    if ( line.length() == 0 ) continue;
		    b = line.indexOf(" ");
		    alignment_.add( line.substring(b).trim() , line.substring(0 , b).trim() );
		}
	    }
	} catch ( NullPointerException e ) { 
	    if ( ! label.equals("phylip") )
		alignment_.add(sequence.toString() , label); 
	} finally { in.close(); }
	if ( alignment_.size() * alignment_.length() == 0 ) {
	    System.out.println("   problem with the alignment inside the file " + infile.getName());
	    System.exit(0);
	}
	//##################################
	//#####   storing alignment   ######
	//##################################
	alignment = new Alignment();
	i = -1;
	switch (characterState) {
	case AA:  
	    System.out.print("   Amino acid sequence alignment "); 
	    alignment = new AminoAcidAlignment(); 
	    while ( ++i < alignment_.size() ) {
		if ( ! alignment.add(AminoAcidAlignment.filter(alignment_.getSequence(i)) , alignment_.getLabel(i)) ) {
		    System.out.println("");
		    System.out.println("   problem with sequence " + (i+1) + " : " + alignment_.getLabel(i));
		    System.exit(0);
		}
	    }
	    break;
	case DNA: 
	    System.out.print("   DNA sequence alignment ");
	    alignment = new NucleotideAlignment();
	    while ( ++i < alignment_.size() ) {
		if ( ! alignment.add(NucleotideAlignment.filter(alignment_.getSequence(i)) , alignment_.getLabel(i)) ) {
		    System.out.println("");
		    System.out.println("   problem with sequence " + (i+1) + " : " + alignment_.getLabel(i));
		    System.exit(0);
		}
	    }
	    break;
	case COD: 
	    System.out.print("   Codon sequence alignment ");
	    alignment = new CodonAlignment();
	    while ( ++i < alignment_.size() ) {
		if ( ! alignment.add(CodonAlignment.filter(alignment_.getSequence(i)) , alignment_.getLabel(i)) ) {
		    System.out.println("");
		    System.out.println("   problem with sequence " + (i+1) + " : " + alignment_.getLabel(i));
		    if ( alignment_.getSequence(i).length() % 3 != 0 )
			System.out.println("   this does not contain trinucleotide character states : length = 3 * " 
					   +  (alignment_.getSequence(i).length()/3) 
					   + " + " + (alignment_.getSequence(i).length()%3) );
		    System.exit(0);
		}
	    }
	    break;
	}
	System.out.println(infile.getName());
	alignment_ = null;
	System.out.println("   before : " + alignment.size() + " sequences / " + alignment.length() + " characters");




	keepRow = new BitSet(alignment.size()); i = -1; while ( ++i < alignment.size() ) keepRow.set(i); 
	do {
	    //##############################################################################
	    //#####   storing a copy of the initial alignment to compute parameters   ######
	    //##############################################################################
	    alignment_ = new Alignment();
	    i = -1;
	    switch (characterState) {
	    case AA:  alignment_ = alignment.toAminoAcidAlignment(); break;
	    case DNA: alignment_ = alignment.toNucleotideAlignment(); break;
	    case COD: alignment_ = alignment.toCodonAlignment(); break;
	    }
	    b = 0; i = -1; while ( ++i < keepRow.size() ) if ( ! keepRow.get(i) ) alignment_.removeRow(b); else b++;
	    //######################################################
	    //#####   computing entropy values and gap rates  ######
	    //######################################################
	    bb = alignment_.length();
	    entropy = new ArrayList<Double>(alignment_.length()); smoothedEntropy = new ArrayList<Double>(bb); 
	    j = bb; while ( --j >= 0 ) { entropy.add( new Double(0) ); smoothedEntropy.add( new Double(0) ); }
	    ok = (0 < hmin) || (hmax < 1.0) || ( stationarity != NO ) || outputOption.contains("-oh") || outputOption.contains("-ch");
	    if ( ok || (colGapRate < 1.0) ) {
		gapCol = new ArrayList<Double>(alignment_.length());
		g = alignment_.size(); maxGapFrequency = (g - 1.0) / g; as = alignment_.getAlphabetSize();
		j =-1;
		while ( ++j < bb ) {
		    if ( j % 100 == 0 ) System.out.print(" " + toString(100.0*((double)j)/((double)bb) , 2) + "% \r");
		    g = alignment_.getColGapRate(j); gapCol.add(new Double(g));
		    if ( (! ok) || (g >= maxGapFrequency) ) continue;
		    frequency = alignment_.getFrequencies(j); 
		    if ( pam0 ) {
			try { h = 0; c = as; while ( (d=frequency[--c]) != 1.0 ) h -= ( d == 0 ) ? 0 : d * log( as , d ); }
			catch ( ArrayIndexOutOfBoundsException ee ) { entropy.set(j , new Double(h)); smoothedEntropy.set(j , new Double(h)); }
			continue;
		    }
		    density = new Matrix(as , as);
		    try { c = as; while ( (d=frequency[--c]) != 1.0 ) density.set(c , c , d); }
		    catch ( ArrayIndexOutOfBoundsException e ) { 
			density = density.times(similarity); density = density.times(1.0 / density.trace());
			eigenvalue = round(density.eig().getRealEigenvalues());
			try { h = 0; c = as; while ( true ) h -= eigenvalue[--c] * log( as , eigenvalue[c] ); }
			catch ( ArrayIndexOutOfBoundsException ee ) { entropy.set(j , new Double(h)); smoothedEntropy.set(j , new Double(h)); }			
		    }
		}
	    }
	    alignment_ = null;
	    
	    //#########################
	    //#####   smoothing   #####
	    //#########################
	    if ( (slidingWindow > 1) && (hmax < 1.0) ) {
		smoothedEntropy = new ArrayList<Double>(alignment.length());
		j = -1;
		while ( ++j < alignment.length() ) {
		    up = 0; down = 0; c = Math.max(-1 , j - slidingWindow/2 - 1);
		    while ( ++c < Math.min(j + slidingWindow/2 + 1 , alignment.length()) ) {
			up += (1.0 - gapCol.get(c).doubleValue()) * entropy.get(c).doubleValue(); 
			down += (1.0 - gapCol.get(c).doubleValue());
		    }
		    if ( down == 0 ) down = 1e-5;
		    smoothedEntropy.add( new Double(up/down) );
		}
	    }
	    //################################################################################
	    //##### computing blocks by distinguishing among uncertainty and variability #####
	    //################################################################################
	    //##### looking for the conserved positions #####
	    keepCol = new BitSet(alignment.length());
	    j = -1; while ( ++j < alignment.length() ) if ( (smoothedEntropy.get(j).doubleValue() < hmax) /*&& (entropy.get(j).doubleValue() < hmax)*/ ) keepCol.set(j);
	    //##### looking for the left- and right- limit of each successive conserved and nonconserved blocks #####
	    if ( hmax < 1.0 ) {
		ok = true;
		while ( ! ok ) { 
		    ok = false;
		    blockStart = new ArrayList<Integer>(0); blockEnd = new ArrayList<Integer>(0);
		    j = -1; while ( ! keepCol.get(++j) ) {}               // the first block must be conserved
		    blockStart.add(new Integer(j)); 
		    while ( ++j < alignment_.length() ) 
			if ( keepCol.get(j-1) != keepCol.get(j) ) {
			    blockEnd.add(new Integer(j-1)); blockStart.add(new Integer(j));
			}
		    if ( ! keepCol.get(j-1) ) blockStart.remove(blockStart.size() - 1); // the last block must be conserved
		    else blockEnd.add(new Integer(j-1));
		    i = -1;
		    while ( (i += 2) < blockStart.size() - 1 ) {
			lgth = blockEnd.get(i).intValue() - blockStart.get(i).intValue();
			left = Math.max(0 , blockEnd.get(i-1).intValue() - lgth);
			right = Math.min(alignment_.length() , blockStart.get(i+1).intValue() + lgth) - 1;
			g = 0; j = left - 1; while ( ++j <= right ) g += gapCol.get(j).doubleValue();
			if ( g < 0.3 ) {
			    j = blockStart.get(i-1).intValue() - 1; up = 0; down = 0;   
			    while ( ++j <= blockEnd.get(i+1).intValue() ) {
				up += (1 - gapCol.get(j).doubleValue()) * entropy.get(j).doubleValue();
				down += (1 - gapCol.get(j).doubleValue());
			    }
			    h = up / down;
			    if ( h <= hmax ) {   // the two flanking conserved blocks and the nonconserved one 
				j = blockStart.get(i).intValue()-1;    //          have a good average entropy
				while ( ++j <= blockEnd.get(i).intValue() ) keepCol.set(j);
				ok = true;
			    }
			}
		    }
		}
	    }
	    //##########################################
	    //#####   removing gapped positions   ######
	    //##########################################
	    if ( colGapRate < 1.0 ) {
		j = -1; while ( ++j < alignment.length() ) if ( gapCol.get(j).doubleValue() > colGapRate ) keepCol.set(j , false);
	    }
	    //#################################
	    //##### removing short blocks #####
	    //#################################
	    if ( minBlockSize > 0 ) {
		blockStart = new ArrayList<Integer>(0); blockEnd = new ArrayList<Integer>(0);
		cpt1 = 0; j = -1;
		while ( ++j < alignment.length() ) {
		    if ( keepCol.get(j) && (cpt1 == 0) ) {blockStart.add(new Integer(j)); cpt1++;}
		    if ( (! keepCol.get(j)) && (cpt1 > 0) ) {blockEnd.add(j); cpt1 = 0;}
		}
		if ( blockStart.size() > blockEnd.size() ) blockEnd.add(new Integer(--j));
		i = -1;
		while ( ++i < blockStart.size() ) {
		    if ( blockEnd.get(i).intValue() - blockStart.get(i).intValue() < minBlockSize ) {
			j = blockStart.get(i).intValue() - 1; 
			while ( ++j < blockEnd.get(i).intValue() ) keepCol.set(j , false);
		    }
		}
	    }
	    //###########################################
	    //##### removing too constant positions #####
	    //###########################################
	    if ( hmin > 0 ) {
		j = -1; while ( ++j < alignment.length() ) if ( smoothedEntropy.get(j).doubleValue() <= hmin ) keepCol.set(j , false);
	    }
	    //##########################################
	    //#####   cleaning alignment          ######
	    //##########################################
	    alignment_ = new Alignment();
	    switch (characterState) {
	    case AA:  alignment_ = alignment.toAminoAcidAlignment(); break;
	    case DNA: alignment_ = alignment.toNucleotideAlignment(); break;
	    case COD: alignment_ = alignment.toCodonAlignment(); break;
	    }
	    j = -1; 
	    while ( ++j < alignment.length() ) 
		if ( ! keepCol.get(j) ) { i = -1; while ( ++i < alignment.size() ) alignment_.setCharAt(i , j , '-'); }
	    //##########################################
	    //#####   looking for stationarity    ######
	    //##########################################
	    if ( stationarity != NO ) {
		System.out.print("\r" + BLANK.substring(0 , 100) + "\r   after :  " 
				 + keepRow.cardinality() + " sequences / " 
				 + keepCol.cardinality() + " characters");
		// computing initial Stuart test X_ij values
		stuartMax = 0;
		switch (characterState) {
		case AA:  stuartMax = KHI2DOF19; break;
		case COD: 
		case DNA: stuartMax = KHI2DOF3; break;
		}
		down = alignment.size() * (alignment.size() - 1) / 2;
		alignmentMeter = new ArrayList<ArrayList<Meter>>(alignment.size());
		alignmentStuartValue = new ArrayList<ArrayList<Double>>(alignment.size()) ;
		ok = true; up = 0; i = -1;
		while ( ++i < alignment.size() ) {
		    alignmentMeter.add(new ArrayList<Meter>(i)); 
		    alignmentStuartValue.add(new ArrayList<Double>(i)); 
		    j = -1; 
		    while ( ++j < i ) {
			alignmentMeter.get(i).add( new Meter() ); 
			alignmentStuartValue.get(i).add( new Double(0) ); 
			if ( keepRow.get(i) && keepRow.get(j) ) {
			    switch (characterState) {
			    case AA:  alignmentMeter.get(i).set(j , new AminoAcidMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));  break;
			    case DNA: alignmentMeter.get(i).set(j , new NucleotideMeter(alignment_.getSequence(i) , alignment_.getSequence(j))); break;
			    case COD: alignmentMeter.get(i).set(j , new CodonMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));      break;
			    }
			    //System.out.println(i + " " + j);
			    //s1 = -1; while ( ++s1 < 20 ) { s2 = -1; while ( ++s2 < 20 ) System.out.print(alignmentMeter.get(i).get(j).get(s1,s2) + " "); } 
			    //System.out.println("");
			    alignmentStuartValue.get(i).set(j , new Double(alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
			    if ( alignmentStuartValue.get(i).get(j).doubleValue() > stuartMax ) ok = false; else up++;
			}
		    }
		}

		/*System.out.println("");
		  i = -1;
		  while ( ++i < alignment.size() ) {
		  j = -1; while ( ++j < i ) System.out.print(" " + String.format(Locale.ENGLISH , "%.5f" , new Double(alignmentStuartValue.get(i).get(j))));
		  System.out.println("");
		  }*/

		// performing stationary-based trimming
		if ( ! ok ) {
		    alignmentGamma = new ArrayList<ArrayList<Meter>>(alignment.size());
		    i = -1;
		    while ( ++i < alignment.size() ) {
			alignmentGamma.add(new ArrayList<Meter>(i));
			j = -1; 
			while ( ++j < i ) {
			    switch (characterState) {
			    case AA:  alignmentGamma.get(i).add( new AminoAcidMeter() );  break;
			    case DNA: alignmentGamma.get(i).add( new NucleotideMeter() ); break;
			    case COD: alignmentGamma.get(i).add( new CodonMeter() );      break;
			    }
			}
		    }

		    switch (stationarity) {
		    case YES:
			while ( ! ok ) {
			    System.out.print("\r" + BLANK.substring(0 , 100) + "\r   after :  " 
					     + keepRow.cardinality() + " sequences / " 
					     + keepCol.cardinality() + " characters"
					     + "   [" + ((int)(100*up/down)) + "%]");
			    
			    // computing gamma
			    i = -1;
			    while ( ++i < alignment.size() ) {
				if ( keepRow.get(i) ) {
				    j = -1; 
				    while ( ++j < i ) {
					if ( keepRow.get(j) ) {
					    //b = alignmentGamma.get(i).get(j).size();
					    switch (characterState) {
					    case DNA:  
						for ( char nt1 : NT ) {
						    if ( (s1=getIndex(DNA , nt1)) == -1 ) continue;
						    for ( char nt2 : NT ) {
							if ( nt1 == nt2 ) continue;
							if ( (s2=getIndex(DNA , nt2)) == -1 ) continue;
							if ( alignmentMeter.get(i).get(j).get(nt1,nt2) <= 0 ) continue;

							alignmentMeter.get(i).get(j).decrement(nt1 , nt2);
							alignmentGamma.get(i).get(j).set(s1 , s2 , alignmentStuartValue.get(i).get(j).doubleValue() * (alignmentStuartValue.get(i).get(j).doubleValue() - alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
							alignmentMeter.get(i).get(j).increment(nt1 , nt2);
						    }
						}
						break;
					    case AA:  
						for ( char nt1 : AM ) {
						    if ( (s1=getIndex(AA , nt1)) == -1 ) continue; 
						    for ( char nt2 : AM ) {
							if ( nt1 == nt2 ) continue;
							if ( (s2=getIndex(AA , nt2)) == -1 ) continue;
							if ( alignmentMeter.get(i).get(j).get(nt1,nt2) <= 0 ) continue;

							alignmentMeter.get(i).get(j).decrement(nt1 , nt2);
							alignmentGamma.get(i).get(j).set(s1 , s2 , alignmentStuartValue.get(i).get(j).doubleValue() * (alignmentStuartValue.get(i).get(j).doubleValue() - alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
							alignmentMeter.get(i).get(j).increment(nt1 , nt2);
						    }
						}
						break;
					    }
					}
				    }
				}
			    }
			    
			    
			    // computing sigma
			    sigmaMax = 0; b_sigmaMax = -1; b = -1;
			    while ( ++b < alignment.length() ) {
				if ( keepCol.get(b) && (entropy.get(b).doubleValue() > 0) ) {
				    sigma = 0; i = -1;
				    while ( ++i < alignment.size() ) {
					if ( keepRow.get(i) && ((c1 = alignment.charAt(i,b)) != '?') && (c1 != '-') ) {
					    j = -1;
					    while ( ++j < i ) 
						if ( keepRow.get(j) && ((c2 = alignment.charAt(j,b)) != '?') 
						     && (c2 != '-') && (c2 != c1) ) 
						    sigma += alignmentGamma.get(i).get(j).get(c1 , c2);
					}
				    }
				    /*if ( b % 10000 == 0 )
				      System.out.print(" [" + b + "," + b_sigmaMax
				      + "," + String.format(Locale.ENGLISH , "%.5f" , new Double(sigma))
				      + "," + String.format(Locale.ENGLISH , "%.5f" , new Double(sigmaMax))
				      + "]");*/
				    if ( sigma > sigmaMax ) { sigmaMax = sigma; b_sigmaMax = b; }
				}
			    }
			    
			    // removing the worst character
			    b = b_sigmaMax; keepCol.set(b , false);
			    
			    // updating
			    /*ok = true; recompute = false; i = -1;
			      while ( ++i < alignment.size() ) {
			      if ( keepRow.get(i) && ((c1 = alignment.charAt(i,b)) != '?') && (c1 != '-') ) {
			      j = -1;
			      while ( ++j < i ) 
			      if ( keepRow.get(j) && ((c2 = alignment.charAt(j,b)) != '?') 
			      && (c2 != '-') && (c1 != c2) ) {
			      alignmentMeter.get(i).get(j).decrement(c1 , c2);
			      //d = alignmentStuartValue.get(i).get(j).doubleValue() / Math.exp(alignmentGamma.get(i).get(j).get(c1 , c2));
			      d = alignmentStuartValue.get(i).get(j).doubleValue() - alignmentGamma.get(i).get(j).get(c1 , c2) / alignmentStuartValue.get(i).get(j).doubleValue();
			      //d = alignmentStuartValue.get(i).get(j).doubleValue();
			      //d -= Math.exp( Math.log(alignmentGamma.get(i).get(j).get(c1 , c2)) - Math.log(d) );
			      alignmentStuartValue.get(i).set(j , new Double(d));
			      if ( d > stuartMax ) ok = false;
			      if ( (d < 0) || Double.isNaN(d) ) recompute = true;
			      }
			      
			      }			
			      }*/
			    i = -1; while ( ++i < alignment.size() ) alignment_.setCharAt(i , b , '-'); 
			    ok = true; d = 0; i = -1;
			    while ( ++i < alignment.size() ) {
				if ( keepRow.get(i) ) {
				    j = -1; 
				    while ( ++j < i ) {
					if ( keepRow.get(j) ) {
					    switch (characterState) {
					    case AA:  alignmentMeter.get(i).set(j , new AminoAcidMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));  break;
					    case DNA: alignmentMeter.get(i).set(j , new NucleotideMeter(alignment_.getSequence(i) , alignment_.getSequence(j))); break;
					    case COD: alignmentMeter.get(i).set(j , new CodonMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));      break;
					    }
					    alignmentStuartValue.get(i).set(j , new Double(alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
					    if ( alignmentStuartValue.get(i).get(j).doubleValue() > stuartMax ) ok = false; 
					    else d++;					
					}
				    }
				}
			    }
			    up = Math.max(d , up);
			}
			break;

		    case FAST:
			inv = 0; b = alignment.length(); while ( --b >= 0 ) inv += (keepCol.get(b) && entropy.get(b).doubleValue() == 0) ? 1 : 0;
			while ( ! ok ) {
			    System.out.print("\r" + BLANK.substring(0 , 100) + "\r   after :  " 
					     + keepRow.cardinality() + " sequences / " 
					     + keepCol.cardinality() + " characters"
					     + "   [" + ((int)(100*up/down)) + "%]");
			    // computing gamma
			    i = alignment.size();
			    while ( --i >= 0 ) {
				if ( ! keepRow.get(i) ) continue;
				//System.out.print(i + " " );
				j = i; 
				while ( --j >= 0 ) {
				    if ( ! keepRow.get(j) ) continue;
				    d = alignmentStuartValue.get(i).get(j).doubleValue();
				    switch (characterState) {
				    case DNA:  
					for ( char nt1 : NT ) {
					    if ( (s1=getIndex(DNA , nt1)) == -1 ) continue;
					    for ( char nt2 : NT ) {
						if ( (nt1 == nt2)
						     || (alignmentMeter.get(i).get(j).get(nt1,nt2) <= 0)
						     || ((s2=getIndex(DNA , nt2)) == -1) ) continue;
						alignmentMeter.get(i).get(j).decrement(nt1 , nt2);
						alignmentGamma.get(i).get(j).set(s1 , s2 , d * (d - alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
						alignmentMeter.get(i).get(j).increment(nt1 , nt2);
					    }
					}
					break;
				    case AA:  
					for ( char nt1 : AM ) {
					    for ( char nt2 : AM ) {
						if ( alignmentMeter.get(i).get(j).get(nt1,nt2) <= 0 ) continue;
						if ( nt1 == nt2 ) continue;
						
						alignmentMeter.get(i).get(j).decrement(nt1 , nt2);
						alignmentGamma.get(i).get(j).set(getIndex(AA , nt1) , getIndex(AA , nt2) , d * (d - alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
						alignmentMeter.get(i).get(j).increment(nt1 , nt2);
					    }
					}
					break;
				    }
				}
			    }
			    // computing sigma
			    sigmaArray = new ArrayList<Double>(alignment.length());
			    sigmaArray2 = new ArrayList<Double>(alignment.length());
			    b = -1;
			    while ( ++b < alignment.length() ) {
				sigmaArray.add( new Double(0) );
				if ( ! (keepCol.get(b) && (entropy.get(b).doubleValue() > 0)) ) continue;
				sigma = 0; i = alignment.size();
				while ( --i >= 0 ) {
				    if ( keepRow.get(i) && ((c1 = alignment.charAt(i,b)) != '?') && (c1 != '-') ) {
					j = -1;
					while ( ++j < i ) 
					    sigma += ( keepRow.get(j) 
						       && ((c2 = alignment.charAt(j,b)) != '?') 
						       && (c2 != '-') 
						       && (c2 != c1) ) ? alignmentGamma.get(i).get(j).get(c1 , c2) : 0;
				    }
				}
				sigmaArray.set(b , new Double(sigma)); sigmaArray2.add(new Double(sigma));
			    }
			    // removing the worst characters
			    Collections.sort(sigmaArray2); 
			    cpt1 = Math.max(2 , (int)((keepCol.cardinality()/*-inv*/)*(down-up) / (down*FAST_STATIONARY_DIV)) ); cpt2 = 0; //System.out.println("" + cpt1);
			    b = sigmaArray2.size();
			    while ( --b >= 0 ) {
				sigmaMax = sigmaArray2.get(b).doubleValue(); j = -1;
				while ( (++j < alignment.length()) && (cpt1 > cpt2) ) {
				    if ( keepCol.get(j) && (sigmaArray.get(j).doubleValue() == sigmaMax) ) {
					keepCol.set(j , false);
					i = alignment.size(); while ( --i >= 0 ) alignment_.setCharAt(i , j , '-'); 
					++cpt2;
				    }
				}
			    }
			    // updating
			    ok = true; d = 0; i = -1;
			    while ( ++i < alignment.size() ) {
				if ( ! keepRow.get(i) ) continue;
				j = -1; 
				while ( ++j < i ) {
				    if ( ! keepRow.get(j) ) continue;
				    switch (characterState) {
				    case AA:  alignmentMeter.get(i).set(j , new AminoAcidMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));  break;
				    case DNA: alignmentMeter.get(i).set(j , new NucleotideMeter(alignment_.getSequence(i) , alignment_.getSequence(j))); break;
				    case COD: alignmentMeter.get(i).set(j , new CodonMeter(alignment_.getSequence(i) , alignment_.getSequence(j)));      break;
				    }
				    alignmentStuartValue.get(i).set(j , new Double(alignmentMeter.get(i).get(j).getStuartMarginalSymmetryTestValue()));
				    if ( alignmentStuartValue.get(i).get(j).doubleValue() > stuartMax ) { ok = false; /*System.out.println("[" + i + "-" + j + "]" + alignmentStuartValue.get(i).get(j).doubleValue());*/ }
				    else d++;					
				}
			    }
			    up = Math.max(d , up);
			}
			
			break;
		    }
		}

	    }
	
	    ok = true; i = -1; 
	    while ( ++i < keepRow.size() ) 
		if ( keepRow.get(i) && (alignment_.getRowGapRate(i) >= rowGapRate) ) {
		    keepRow.set(i , false); ok = false;
		}
	    if ( keepRow.cardinality() == 0 ) ok = true;
	} while ( ! ok ); // if ok=false here, then one must restart the whole procedure 
	                  // because at least one sequence was removed















	System.out.print("\r" + BLANK.substring(0 , 100) + "\r   after :  " 
			 + keepRow.cardinality() + " sequences / " 
			 + keepCol.cardinality() + " characters");
	
	//##########################################
	//#####   cleaning alignment          ######
	//##########################################
	alignmentCleaned = new Alignment(); alignmentComplement = new Alignment();
	switch (characterState) {
	case AA: 
	    alignmentCleaned = alignment.toAminoAcidAlignment(); if ( outputOption.contains("-c") ) alignmentComplement = alignment.toAminoAcidAlignment(); 
	    break;
	case DNA: 
	    alignmentCleaned = alignment.toNucleotideAlignment(); if ( outputOption.contains("-c") ) alignmentComplement = alignment.toNucleotideAlignment(); 
	    break;
	case COD: 
	    alignmentCleaned = alignment.toCodonAlignment(); if ( outputOption.contains("-c") ) alignmentComplement = alignment.toCodonAlignment(); 
	    break;
	}
	j = keepCol.size(); while ( --j >= 0 ) if ( ! keepCol.get(j) ) alignmentCleaned.removeColumn(j);
	i = keepRow.size(); while ( --i >= 0 ) if ( ! keepRow.get(i) ) alignmentCleaned.removeRow(i);
	if ( outputOption.contains("-c") ) {
	    j = keepCol.size(); while ( --j >= 0 ) if ( keepCol.get(j) ) alignmentComplement.removeColumn(j);
	    i = keepRow.size(); while ( --i >= 0 ) if ( ! keepRow.get(i) ) alignmentComplement.removeRow(i);
	}
	/*
	  b =0; j = -1; 
	  while ( ++j < keepCol.size() ) 
	  if ( ! keepCol.get(j) ) alignmentCleaned.removeColumn(b);
	  else b++;
	  b = 0; i = -1;
	  while ( ++i < keepRow.size() ) 
	  if ( ! keepRow.get(i) ) alignmentCleaned.removeRow(b);
	  else b++;
	  b =0; j = -1; 
	  while ( ++j < keepCol.size() ) 
	  if ( keepCol.get(j) ) alignmentComplement.removeColumn(b);
	  else b++;
	  if ( outputOption.contains("-c") ) {
	  b = 0; i = -1;
	  while ( ++i < keepRow.size() ) 
	  if ( ! keepRow.get(i) ) alignmentComplement.removeRow(b);
	  else b++;
	  }
	*/
	if ( alignmentCleaned.length() != -1 )
	    System.out.println("\r                                          \r   after :  " 
			       + alignmentCleaned.size() + " sequences / " + alignmentCleaned.length() + " characters");
	else
	    System.out.println("\r                                          \r   after :  " 
			       + alignmentCleaned.size() + " sequences / 0 characters");
	
	


	//##################################
	//#####   writing alignment   ######
	//##################################
	bb = -1;
	while ( ++bb < outputOption.size() ) {
	    line = outputOption.get(bb).toLowerCase();
	    alignment_ = new Alignment();
	    if ( line.startsWith("-o") ) {
		switch (characterState) {
		case AA:  alignment_ = alignmentCleaned.toAminoAcidAlignment(); break;
		case DNA: alignment_ = alignmentCleaned.toNucleotideAlignment(); break;
		case COD: alignment_ = alignmentCleaned.toCodonAlignment(); break;
		}
	    }
	    if ( line.startsWith("-c") ) {
		switch (characterState) {
		case AA:  alignment_ = alignmentComplement.toAminoAcidAlignment(); break;
		case DNA: alignment_ = alignmentComplement.toNucleotideAlignment(); break;
		case COD: alignment_ = alignmentComplement.toCodonAlignment(); break;
		}
	    }
	    //##### output format #####
	    outputFormat = PHYLIP_ACCN;
	    if ( line.indexOf("p") != -1 )   outputFormat = PHYLIP;
	    if ( line.indexOf("pp") != -1 )  outputFormat = PHYLIP_TAX;
	    if ( line.indexOf("ppp") != -1 ) outputFormat = PHYLIP_ACCN;
	    if ( line.indexOf("f") != -1 )   outputFormat = FASTA;
	    if ( line.indexOf("n") != -1 )   outputFormat = PAUP;
	    if ( line.indexOf("nn") != -1 )  outputFormat = PAUP_TAX;
	    if ( line.indexOf("nnn") != -1 ) outputFormat = PAUP_ACCN;
	    if ( line.indexOf("h") != -1 )   outputFormat = HTML;
	    //##### codon position(s) #####
	    codonCode = 0;
	    if ( line.indexOf("1") != -1 ) codonCode += 1;
	    if ( line.indexOf("2") != -1 ) codonCode += 2;
	    if ( line.indexOf("3") != -1 ) codonCode += 4;
	    if ( codonCode == 0 ) codonCode = 7;
	    //##### translation into another character type #####
	    outputType = characterState;
	    if ( line.indexOf("aa") != -1 ) {
		outputType = AA;
		if ( characterState == DNA ) {
		    System.out.println("   impossible to convert into amino acid sequence alignment"); System.exit(0);
		}
		if ( characterState != AA ) alignment_ = alignment_.toAminoAcidAlignment(); 
	    }
	    else {
		if ( line.indexOf("co") != -1 ) {
		    if ( characterState == DNA ) {
			System.out.println("   impossible to convert into codon sequence alignment"); System.exit(0);
		    }
		    if ( codonCode == 7 ) { outputType = COD; if ( characterState != COD ) alignment_ = alignment_.toCodonAlignment(); }
		    else {
			if ( characterState != COD ) alignment_ = alignment_.toCodonAlignment();
			outputType = DNA; alignment_ = alignment_.toNucleotideAlignment(codonCode);
		    }
		}
		else {
		    if ( line.indexOf("dna") != -1 ) {
			outputType = DNA; 
			if ( codonCode == 7 ) alignment_ = alignment_.toNucleotideAlignment();
			else {
			    if ( characterState == DNA ) {
				System.out.println("   impossible to convert into codon sequence alignment"); System.exit(0);
			    }
			    if ( characterState != COD ) alignment_ = alignment_.toCodonAlignment();
			    alignment_ = alignment_.toNucleotideAlignment(codonCode);
			}
		    }
		    else {
			if ( line.indexOf("ry") != -1 ) {
			    outputType = DNA;
			    if ( codonCode == 7 ) alignment_ = alignment_.toRYcodingAlignment();
			    else {
				if ( characterState == DNA ) {
				    System.out.println("   impossible to convert into codon sequence alignment"); System.exit(0);
				}
				if ( characterState != COD ) alignment_ = alignment_.toCodonAlignment();
				alignment_ = alignment_.toNucleotideAlignment(codonCode);
				alignment_ = alignment_.toRYcodingAlignment();
			    }
			}
		    }
		}
	    }
	    //##### rewriting taxon names for PHYLIP, PHYLIP_RW, PAUP or PAUP_RW output files #####
	    rwLabel = new ArrayList<String>(0);
	    if ( (outputFormat == PHYLIP) || (outputFormat == PAUP) ) {
		rwLabel = new ArrayList<String>(alignment_.size());
		i = -1; while ( ++i < alignment_.size() ) rwLabel.add(noBlank(alignment_.getLabel(i).trim()));
	    }
	    if ( (outputFormat == PHYLIP_TAX) || (outputFormat == PAUP_TAX) ) {
		rwLabel = new ArrayList<String>(alignment_.size());
		i = -1;
		while ( ++i < alignment_.size() ) {
		    rwLabel.add(""); line = alignment_.getLabel(i).trim();
		    if ( line.endsWith("]") ) {
			b = line.lastIndexOf("["); c = line.lastIndexOf("]"); 
			if ( (b != -1) && (b < c) ) rwLabel.set(i , noBlank(line.substring(b+1 , c)));
			else
			    rwLabel.set(i , noBlank(line));
		    }
		    else
			rwLabel.set(i , noBlank(line));
		}

	    }
	    if ( (outputFormat == PHYLIP_ACCN) || (outputFormat == PAUP_ACCN) ) {
		rwLabel = new ArrayList<String>(alignment_.size());
		i = -1;
		while ( ++i < alignment_.size() ) {
		    rwLabel.add(""); line = alignment_.getLabel(i).trim();
		    if ( line.endsWith("]") ) {
			b = line.lastIndexOf("["); c = line.lastIndexOf("]"); 
			if ( (b != -1) && (b < c) ) {
			    rwLabel.set(i , noBlank(line.substring(b+1 , c)));
			    if ( (c = line.lastIndexOf("|")) != -1 ) {
				while ( line.charAt(--c) == '|' ) {}
				b = ++c; while ( --b >= 0 ) if ( line.charAt(b) == '|' ) break;
				if ( (b != -1) && (b+1 != c) ) rwLabel.set(i , rwLabel.get(i) + "_____" + noBlank(line.substring(b+1 , c)));
				else rwLabel.set(i , rwLabel.get(i) + "_____" + noBlank(line.substring(0 , line.lastIndexOf("[")).trim()));
			    }
			    else
				rwLabel.set(i , rwLabel.get(i) + "_____" + noBlank(line.substring(0 , line.lastIndexOf("[")).trim()));
			}
			else
			    rwLabel.set(i , noBlank(line.trim()));
		    }
		    else
			rwLabel.set(i , noBlank(line.trim()));
		}
	    }
	    if ( rwLabel.size() > 0 ) {
		i = -1;
		while ( ++i < rwLabel.size() ) {
		    line = rwLabel.get(i);
		    //System.out.println(line);
		    if ( i != rwLabel.lastIndexOf(line) ) {
			j = i-1;
			while ( ++j < rwLabel.size() ) {
			    label = rwLabel.get(j);
			    if ( line.equals(label) ) {
				label = alignment_.getLabel(j).trim();
				b = label.lastIndexOf("["); c = label.lastIndexOf("]"); 
				if ( (b != -1) && (c != -1) && (b < c) 
				     && ((c = label.lastIndexOf("|")) != -1) ) {
				    //System.out.println(label + " " + c);
				    while ( label.charAt(--c) == '|' ) {}
				    b = ++c; while ( --b >= 0 ) if ( label.charAt(b) == '|' ) break;
				    if ( (b != -1) && (b+1 != c) ) rwLabel.set(j , rwLabel.get(j) 
									       + "_____" + noBlank(label.substring(b+1 , c)));
				    else 
					rwLabel.set(j , rwLabel.get(j) 
						    + "_____" + noBlank(label.substring(0 , label.lastIndexOf("[")).trim()));
				}
				else rwLabel.set(j , rwLabel.get(j) + "_____" + noBlank(label.trim()));
			    }
			}
		    }
		}
		i = -1;
		while ( ++i < rwLabel.size() ) {
		    line = rwLabel.get(i);
		    if ( i != rwLabel.lastIndexOf(line) ) {
			b = 0; j = i-1;
			while ( ++j < rwLabel.size() ) if ( line.equals( rwLabel.get(j) ) ) rwLabel.set(j , rwLabel.get(j) + "___" + (++b));
		    }
		}
	    }


	    /*rwLabel = new ArrayList<String>(alignment_.size());
	      i = -1;
	      while ( ++i < alignment_.size() ) {
	      rwLabel.add(""); line = alignment_.getLabel(i); b = line.lastIndexOf("["); c = line.lastIndexOf("]");
	      if ( (b != -1) && (c != -1) && (b < c) ) rwLabel.set(i , noBlank(line.substring(++b , c)));
	      }
	      i = -1; 
	      while ( ++i < rwLabel.size() ) {
	      cpt1 = 1; j = i;
	      while ( ++j < rwLabel.size() ) {
	      if ( rwLabel.get(i).equals( rwLabel.get(j) ) ) {
	      cpt1++;
	      line = alignment_.getLabel(j); b = line.indexOf("|"); c = line.indexOf("|", ++b);
	      if ( (b != -1) && (c != -1) ) {while ( line.indexOf("|", c + 1) != -1 ) {b = ++c; c = line.indexOf("|", b);}}
	      if ( (b != -1) && (c != -1) && (b + 1 < c) ) 
	      rwLabel.set(j , rwLabel.get(j) + "_____" + noBlank(line.substring(b , c)));
	      else {
	      if ( line.indexOf("[") != -1 ) 
	      rwLabel.set(j , rwLabel.get(j) + "_____" + noBlank(line.substring(0 , line.indexOf("["))));
	      else
	      rwLabel.set(j , noBlank(line));
	      }
	      }
	      }
	      if ( cpt1 > 1 ) {
	      line = alignment_.getLabel(i); b = line.indexOf("|"); c = line.indexOf("|", ++b); 
	      if ( (b != -1) && (c != -1) ) {while ( line.indexOf("|", c + 1) != -1 ) {b = ++c; c = line.indexOf("|", b);}}
	      if ( (b != -1) && (c != -1) && (b + 1 < c) ) 
	      rwLabel.set(i , rwLabel.get(i) + "_____" + noBlank(line.substring(b , c)));
	      else {
	      if ( line.indexOf("[") != -1 ) 
	      rwLabel.set(i , rwLabel.get(i) + "_____" + noBlank(line.substring(0 , line.indexOf("["))));
	      else
	      rwLabel.set(i , noBlank(line));
	      }
	      }
	      }
	      i = -1; 
	      while ( ++i < rwLabel.size() ) {
	      if ( rwLabel.get(i).length() > 100 ) rwLabel.set(i , rwLabel.get(i).substring(0 , 100)); /////////
	      cpt1 = 1; j = i;
	      while ( ++j < rwLabel.size() ) 
	      if ( rwLabel.get(i).equals( rwLabel.get(j) ) ) rwLabel.set(j , rwLabel.get(j) + "___" + (++cpt1));
	      if ( cpt1 > 1 ) rwLabel.set(i , rwLabel.get(i) + "___1");
	      }*/
	    //##### writing file #####
	    outfile = outputFileName.get(bb);
	    out = new BufferedWriter(new FileWriter(outfile));
	    switch (outputFormat) {
	    case FASTA:
		i = -1;
		while ( ++i < alignment_.size() ) {
		    out.write(">" + alignment_.getLabel(i)); out.newLine();
		    b = 0;
		    do {
			if ( characterState != COD )
			    out.write(alignment_.getSequence(i).substring(b , Math.min(b+lengthDisplay,alignment_.length()))); 
			else
			    out.write(alignment_.getSequence(i).substring(b , Math.min(b+lengthDisplay,alignment_.length()))); 
			out.newLine();
			b += lengthDisplay;
		    } while ( Math.min(b,alignment_.length()) < alignment_.length() );
		}
		break;

	    case PHYLIP:
	    case PHYLIP_TAX:
	    case PHYLIP_ACCN:
		maxLabelLength = -1; i = -1; 
		while ( ++i < rwLabel.size() ) 
		    maxLabelLength = (maxLabelLength > rwLabel.get(i).length()) ? maxLabelLength : rwLabel.get(i).length();
		out.write(" " + alignment_.size());
		if ( outputType == COD ) out.write(" " + (3*alignment_.length()));
		else out.write(" " + alignment_.length());		
		out.newLine();
		i = -1;
		while ( ++i < alignment_.size() ) {
		    out.write( (rwLabel.get(i) + BLANK).substring(0 , maxLabelLength + 1) );
		    out.write( alignment_.getSequence(i) );
		    out.newLine();
		}
		break;
		
	    case PAUP:
	    case PAUP_TAX:
	    case PAUP_ACCN:
		maxLabelLength = -1; i = -1; 
		while ( ++i < rwLabel.size() ) 
		    maxLabelLength = (maxLabelLength > rwLabel.get(i).length()) ? maxLabelLength : rwLabel.get(i).length();
		out.write("#NEXUS"); out.newLine(); out.newLine();
		out.write("begin data;"); out.newLine();
		out.write("   dimensions ntax=" + alignment_.size());
		if ( outputType == COD ) out.write(" nchar=" + (3*alignment_.length()) + ";");
		else out.write(" nchar=" + alignment_.length() + ";");
		out.newLine();
		if ( outputType == AA ) out.write("   format datatype=PROTEIN;");
		else out.write("   format datatype=NUCLEOTIDE;");
		out.newLine();
		out.write("   matrix"); out.newLine();
		i = -1;
		while ( ++i < alignment_.size() ) {
		    out.write( "     " + (rwLabel.get(i) + BLANK).substring(0 , maxLabelLength + 1) );
		    out.write( alignment_.getSequence(i) );
		    out.newLine();
		}
		out.write("   ;"); out.newLine();
		out.write("end;"); out.newLine();
		break;

	    case HTML: 
		maxLabelLength = -1; i = -1;
		while ( ++i < alignment.size() ) 
		    maxLabelLength = (maxLabelLength > alignment.getLabel(i).length()) ? maxLabelLength : alignment.getLabel(i).length();
		maxLabelLength = (maxLabelLength > 5) ? maxLabelLength : 5;
		out.write("<html>"); out.newLine(); 
		out.write("<body bgcolor=\"#FFFFFF\" text=\"#D3D4D4\">"); out.newLine();
		out.write("<pre>"); out.newLine();
		out.newLine(); out.newLine(); 
		sequence = new StringBuffer("");
		while ( sequence.length() < alignment.length() + 10 ) sequence = sequence.append(BLANK);
		line = sequence.substring(0 , alignment.length() + 10);
		out.write(BLANK.substring(0 , maxLabelLength + 1)); out.write(line); out.newLine();
		//##### histograms #####		
		i = 11;
		while ( --i >= 0 ) {
		    switch (i) {
		    case 10: out.write(BLANK.substring(0 , maxLabelLength - 3) + htmlEntry("1.0_" , "" , "#000000")); break;
		    case 5:  out.write(BLANK.substring(0 , maxLabelLength - 3) + htmlEntry("0.5_" , "" , "#000000")); break;
		    case 0:  out.write(BLANK.substring(0 , maxLabelLength - 3) + htmlEntry("0.0_" , "" , "#000000")); break;
		    default: out.write(BLANK.substring(0 , maxLabelLength + 1)); break;
		    }
		    j = -1;
		    while ( ++j < alignment.length() ) {
			if ( smoothedEntropy.get(j).doubleValue() < ((double)i)/10.0 ) {
			    if ( gapCol.get(j).doubleValue() <= ((double)i)/10.0 ) {
				if ( characterState == COD ) out.write("   ");
				else out.write(" ");
			    }
			    else {
				if ( characterState == COD ) out.write("===");
				else out.write("=");
			    }
			}
			else {
			    if ( smoothedEntropy.get(j).doubleValue() > ((double)(i+0.5))/10.0 )
				if ( characterState == COD ) out.write(htmlEntry("<b>:::</b>" , "" , "#000000"));
				else out.write(htmlEntry("<b>:</b>" , "" , "#000000"));
			    else
				if ( characterState == COD ) out.write(htmlEntry("<b>...</b>" , "" , "#000000"));
				else out.write(htmlEntry("<b>.</b>" , "" , "#000000"));
			}
		    }
		    out.newLine();
		}
		//##### bars #####
		out.newLine(); 
		sequence = new StringBuffer(line);
		lgth = alignment.length(); 
		if ( characterState == COD ) { lgth *= 3; sequence = sequence.append(line).append(line); }
		j = -1; 
		while ( ++j <= lgth ) {
		    if ( j % 5 == 0 ) {
			if ( j % 10 == 0 ) sequence.setCharAt(j , '|');
			else sequence.setCharAt(j , '-');
		    }
		    else sequence.setCharAt(j , '=');
		}
		sequence = sequence.deleteCharAt(0);
		out.write(BLANK.substring(0 , maxLabelLength + 1));
		out.write(htmlEntry("" + sequence.substring(0 , lgth) + "" , "#F0F0F0" , "#000000"));
		out.newLine(); 
		sequence = new StringBuffer(line); if ( characterState == COD ) sequence = sequence.append(line).append(line);
		j = 0;
		while ( ++j < lgth ) {
		    if ( j % 10 == 0 ) {
			label = (new Integer(j)).toString();
			sequence = sequence.replace(--j , (j = --j + label.length()) , label);
		    }
			
		}
		out.write(BLANK.substring(0 , maxLabelLength + 1));
		out.write(htmlEntry(sequence.substring(0 , lgth) , "#F0F0F0" , "#000000"));
		out.newLine(); 
		sequence = new StringBuffer(line); if ( characterState == COD ) sequence = sequence.append(line).append(line);
		j = -1; 
		while ( ++j <= lgth ) {
		    if ( j % 5 == 0 ) {
			if ( j % 10 == 0 ) sequence.setCharAt(j , '|');
			else sequence.setCharAt(j , '-');
		    }
		    else sequence.setCharAt(j , '=');
		}
		sequence = sequence.deleteCharAt(0);
		out.write(BLANK.substring(0 , maxLabelLength + 1));
		out.write(htmlEntry(sequence.substring(0 , lgth) , "#F0F0F0" , "#000000"));
		out.newLine(); 
		sequence = null; label = null; line = null;
		out.newLine(); 
		//##### alignment #####
		i = -1;
		while ( ++i < alignment.size() ) {
		    out.write(htmlEntry((alignment.getLabel(i) + BLANK).substring(0 , maxLabelLength + 1) , "" , "#000000")); 
		    j = -1;
		    while ( ++j < alignment.length() ) {
			if ( characterState == COD ) line = alignment.getCodonAt(i , j);
			else line = alignment.getCharAt(i , j);
			if ( (! keepCol.get(j)) || (! keepRow.get(i)) ) out.write(line);
			else {
			    if ( ((characterState == COD) 
				  && alignment.getMajorityCharacter(j).contains(CodonAlignment.toAminoAcidCharacter(line) + ""))
				 || ((characterState != COD) && alignment.getMajorityCharacter(j).contains(line)) ) 
				out.write(htmlEntry(line , "#000000" , "#FFFFFF")); 
			    else out.write(htmlEntry(line , "" , "#000000"));
			}
		    }
		    out.newLine();
		}
		out.write("<span style=\"color:#000000\">"); out.newLine();
		if ( stationarity != NO ) {
		    out.newLine(); out.newLine(); 
		    out.write("Stuart's (1955) test p-values"); out.newLine(); out.newLine(); 
		    i = -1;
		    while ( ++i < alignment.size() ) {
			out.write((alignment.getLabel(i) + BLANK).substring(0 , maxLabelLength + 1)); 
			j = -1;
			while ( ++j < i ) {
			    switch (characterState) {
			    case AA:  meter = new AminoAcidMeter(alignment.getSequence(i) , alignment.getSequence(j));  break;
			    case DNA: meter = new NucleotideMeter(alignment.getSequence(i) , alignment.getSequence(j)); break;
			    case COD: meter = new NucleotideMeter(alignment.getSequence(i) , alignment.getSequence(j)); break;
			    }
			    pValue = meter.getStuartMarginalSymmetryTestPvalue();
			    out.write("  " + String.format(Locale.ENGLISH , "%.6f" , new Double(pValue)));
			}
			out.newLine();
		    }
		    out.newLine(); out.newLine(); 	    
		    i = -1;
		    while ( ++i < alignment.size() ) {
			out.write((alignment.getLabel(i) + BLANK).substring(0 , maxLabelLength + 1)); 
			j = -1;
			while ( ++j < i ) {
			    switch (characterState) {
			    case AA:  meter = new AminoAcidMeter(alignmentCleaned.getSequence(i) , alignmentCleaned.getSequence(j));  break;
			    case DNA: meter = new NucleotideMeter(alignmentCleaned.getSequence(i) , alignmentCleaned.getSequence(j)); break;
			    case COD: meter = new CodonMeter(alignmentCleaned.getSequence(i) , alignmentCleaned.getSequence(j)); break;
			    }
			    pValue = meter.getStuartMarginalSymmetryTestPvalue();
			    out.write("  " + String.format(Locale.ENGLISH , "%.6f" , new Double(pValue)));
			}
			out.newLine();
		    }
		}
		out.newLine(); out.newLine();
		out.write("Characters : " + keepCol.cardinality() + " selected  " 
			  + (alignment.length()-keepCol.cardinality()) + " removed"); out.newLine();
		out.write("  selected: "); 
		b = -1; while ( ++b < alignment.length() ) if ( keepCol.get(b) ) out.write(" " + (b+1));
		/*
		  left = -1; b = -1;
		  while ( ++b < alignment.length() ) {
		  if ( keepCol.get(b) && (left == -1) ) { left = b+1; out.write(" " + left); }
		  if ( (! keepCol.get(b)) && (left != -1) ) { if ( b-1 > left ) out.write("-" + b); left = -1; }
		  }
		  if ( keepCol.get(--b) && (left != -1) ) if ( b-1 > left ) out.write("-" + b); 
		*/
		out.newLine(); out.write("  removed:  "); 
		b = -1; while ( ++b < alignment.length() ) if ( ! keepCol.get(b) ) out.write(" " + (b+1));
		/*left = -1; b = -1;
		  while ( ++b < alignment.length() ) {
		  if ( (! keepCol.get(b)) && (left == -1) ) { left = b+1; out.write(" " + left); }
		  if ( keepCol.get(b) && (left != -1) ) { if ( b-1 > left ) out.write("-" + b); left = -1; }
		  }
		  if ( (! keepCol.get(--b)) && (left != -1) )  if ( b-1 > left ) out.write("-" + b);*/
		out.newLine(); out.newLine();
		lgth = (int) Math.rint(Math.log10(alignment.length())) + 1;
		out.write(("ch." + BLANK).substring(0 , lgth));
		out.write(" entropy     smooth. entr.     gap rate");
		out.newLine();
		j = -1;
		while ( ++j < alignment.length() ) {
		    b = j+1;
		    line = (b + BLANK).substring(0 , lgth);
		    out.write(line);
		    out.write(" " + String.format(Locale.ENGLISH , "%.6f" , entropy.get(j))); 
		    out.write("       " + String.format(Locale.ENGLISH , "%.6f" , smoothedEntropy.get(j)));
		    out.write("       " + String.format(Locale.ENGLISH , "%.6f" , gapCol.get(j))); 
		    out.newLine();
		}
		



		out.write("</span>"); out.newLine(); 
		out.write("</pre>"); out.newLine(); 
		out.write("</body>"); out.newLine(); 
		out.write("</html>"); out.newLine(); 
		break;
	    }
	    
	    out.close();
	    


	    
	}

	




    }


    public static double[] round(double[] array) {
	array_ = new double[(l_=array.length)]; i_ = l_;
	while ( --i_ >= 0 ) array_[i_] = ( Math.floor(10000000000.0*Math.abs(array[i_])) != 0 ) ? array[i_] : 0;
	return array_;
	/*
	  l_ = array.length; array_ = new double[l_]; Arrays.fill(array_ , 0);
	  i_ = l_;
	  while ( --i_ >= 0 ) if ( Math.floor(10000000000.0*Math.abs(array[i_])) != 0 ) array_[i_] = array[i_];
	  return array_;
	*/
    }

    public static double log(double base , double a) {
	if ( Math.abs(a) == 0 ) return 0;
	return ( base == 4.0 ) ? Math.log(a) / 1.386294361 : ( base == 20.0 ) ? Math.log(a) / 2.995732274 : Math.log(a) / Math.log(base);
    }

    public static String toString(double d , int limit) {
	//return String.format(Locale.ENGLISH , "%." + limit + "f" , new Double(d));
	NumberFormat df = NumberFormat.getInstance(Locale.US); df.setMaximumFractionDigits(limit); df.setMinimumFractionDigits(limit);
	return df.format(d);

    }

    public static String htmlEntry(String entry , String backgroundColor , String textColor) {
	if ( (! backgroundColor.equals("")) && (! textColor.equals("")) )
	    return "<span style=\"background-color:" + backgroundColor + ";color:" + textColor + "\">" + entry + "</span>";
	if ( ! backgroundColor.equals("") )
	    return "<span style=\"background-color:" + backgroundColor + "\">" + entry + "</span>";
	if ( ! textColor.equals("") )
	    return "<span style=\"color:" + textColor + "\">" + entry + "</span>";
	return entry;
    }

    public static String noBlank(String line) { // as in PAUP => ( ) [ ] { } / \ , ; : = *'"`+- < >
	return line.replace(' ' , '_').replace('(' , '_').replace(')' , '_').replace('[' , '_').replace(']' , '_').replace('{' , '_').replace('}' , '_').replace('/' , '_').replace('\\' , '_').replace(',' , '_').replace(';' , '_').replace(':' , '_').replace('=' , '_').replace('*' , '_').replace('\"' , '_').replace('\'' , '_').replace('+' , '_').replace('-' , '_').replace('<' , '_').replace('>' , '_').replace('|' , '_').replace('.' , '_');
    }


    public static boolean neq(String str1 , String str2) {
	return ! str1.equals(str2);
    }

    public static int getIndex(byte whichCharState , char charState) {
	return ( whichCharState == DNA ) 
	    ? ( charState == 'A' ) ? NucleotideMeter.A 
	    : ( charState == 'C' ) ? NucleotideMeter.C 
	    : ( charState == 'G' ) ? NucleotideMeter.G 
	    : ( charState == 'T' ) ? NucleotideMeter.T 
	    : -1
	    : ( whichCharState == AA )
	    ? ( charState == 'A' ) ? AminoAcidMeter.A
	    : ( charState == 'R' ) ? AminoAcidMeter.R
	    : ( charState == 'N' ) ? AminoAcidMeter.N
	    : ( charState == 'D' ) ? AminoAcidMeter.D
	    : ( charState == 'C' ) ? AminoAcidMeter.C
	    : ( charState == 'Q' ) ? AminoAcidMeter.Q
	    : ( charState == 'E' ) ? AminoAcidMeter.E
	    : ( charState == 'G' ) ? AminoAcidMeter.G
	    : ( charState == 'H' ) ? AminoAcidMeter.H
	    : ( charState == 'I' ) ? AminoAcidMeter.I
	    : ( charState == 'L' ) ? AminoAcidMeter.L
	    : ( charState == 'K' ) ? AminoAcidMeter.K
	    : ( charState == 'M' ) ? AminoAcidMeter.M
	    : ( charState == 'F' ) ? AminoAcidMeter.F
	    : ( charState == 'P' ) ? AminoAcidMeter.P
	    : ( charState == 'S' ) ? AminoAcidMeter.S
	    : ( charState == 'T' ) ? AminoAcidMeter.T
	    : ( charState == 'W' ) ? AminoAcidMeter.W
	    : ( charState == 'Y' ) ? AminoAcidMeter.Y
	    : ( charState == 'V' ) ? AminoAcidMeter.V
	    : -1
	    : -1;
	
	/*switch ( whichCharState ) {
	case DNA:
	    switch ( charState ) {
	    case 'A': return NucleotideMeter.A;
	    case 'C': return NucleotideMeter.C;
	    case 'G': return NucleotideMeter.G;
	    case 'T': return NucleotideMeter.T;
	    }
	    return -1;
	case AA:
	    switch ( charState ) {
	    case 'A': return AminoAcidMeter.A;
	    case 'R': return AminoAcidMeter.C;
	    case 'N': return AminoAcidMeter.G;
	    case 'D': return AminoAcidMeter.D;
	    case 'C': return AminoAcidMeter.C;
	    case 'Q': return AminoAcidMeter.Q;
	    case 'E': return AminoAcidMeter.E;
	    case 'G': return AminoAcidMeter.G;
	    case 'H': return AminoAcidMeter.H;
	    case 'I': return AminoAcidMeter.I;
	    case 'L': return AminoAcidMeter.L;
	    case 'K': return AminoAcidMeter.K;
	    case 'M': return AminoAcidMeter.M;
	    case 'F': return AminoAcidMeter.F;
	    case 'P': return AminoAcidMeter.P;
	    case 'S': return AminoAcidMeter.S;
	    case 'T': return AminoAcidMeter.T;
	    case 'W': return AminoAcidMeter.W;
	    case 'Y': return AminoAcidMeter.Y;
	    case 'V': return AminoAcidMeter.V;
	    }
	    return -1;
	}
	return -1;*/
    }


    public static void displayUserGuide() {
	System.out.println("");
	System.out.println(" BMGE comes with ABSOLUTELY NO WARRANTY.  This is free software, and you are welcome to");
	System.out.println(" redistribute it under certain conditions.  See the file COPYING.txt for details.");
	System.out.println("");
	System.out.println(" BMGE (version 1.12) arguments :");
	//System.out.println("");
	System.out.println("   -i <infile> : input file in fasta or phylip sequential format");
	//System.out.println("");
	System.out.println("   -t [AA,DNA,CODON] : sequence coding in the input file (Amino Acid-, DNA-, RNA-, or");
	System.out.println("                       CODON-coding sequences, respectively)");
	//System.out.println("");
	System.out.println("   -m BLOSUM<n> :   for Amino Acid or CODON sequence alignment; name of the BLOSUM");
	System.out.println("                    matrix used to estimate the entropy-like value for each character");
	System.out.println("                    (n = 30, 35, 40, ..., 60, 62, 65, ..., 90, 95; default: BLOSUM62)");
	System.out.println("   -m DNAPAM<n:r> : for DNA or RNA sequence alignment; name of the PAM matrix (n ranges");
	System.out.println("                    from 1 to 10,000) and transition/transvertion ratio r value (r ranges");
	System.out.println("                    from 0 to 10,000) used to estimate the entropy-like value for each");
	System.out.println("                    character (default: DNAPAM100:2)");
	System.out.println("   -m DNAPAM<n> :   same as previous option but with r = 1");
	System.out.println("   -m [ID,PAM0] :   for all sequence coding; identity matrix used to estimate entropy-");
	System.out.println("                    like values for each characters");
	//System.out.println("");
	System.out.println("   -g <rate_max> :          real number corresponding to the maximum gap rate allowed");
	System.out.println("                            per character (ranges from 0 to 1; default: 0.2)");
	System.out.println("   -g <col_rate:row_rate> : real numbers corresponding to the maximum gap rates allowed");
	System.out.println("                            per sequence and character, respectively (range from 0 to 1;");
	System.out.println("                            default: 0:0.2)");
	//System.out.println("");
	System.out.println("   -h <thr_max> :         real number corresponding to the maximum entropy threshold");
	System.out.println("                          (ranges from 0 to 1; default: 0.5)");
	System.out.println("   -h <thr_min:thr_max> : real numbers corresponding to the minimum and maximum entropy");
	System.out.println("                          threshold, respectively (range from 0 to 1; default: 0:0.5)");
	//System.out.println("");
	System.out.println("   -b <min_size> : integer number corresponding to the minimum length of selected");
	System.out.println("                   region(s) (ranges from 1 to alignment length; default: 5)");
	//System.out.println("");
	System.out.println("   -w <size> : sliding window size (must be odd; ranges from 1 to alignment length; if");
	System.out.println("               set to 1, then entropy-like values are not smoothed; default: 3)");
	//System.out.println("");
	System.out.println("   -s [NO,YES] : if set to YES, performs a stationarity-based trimming of the multiple");
	System.out.println("                 sequence alignement (default: NO)");
	//System.out.println("");
	System.out.println("   -o<x> <outfile> : output file in phylip sequential (-o, -op, -opp, -oppp), fasta (-of),");
	System.out.println("                     nexus (-on, -onn, -onnn) or html (-oh) format; for phylip and nexus");
	System.out.println("                     format, options -opp and -onn allow NCBI-formatted sequence names to");
	System.out.println("                     be renamed onto their taxon name only; options -oppp and -onnn allow");
	System.out.println("                     renaming onto 'taxon name'_____'accession number' (default: -oppp)");
	System.out.println("   -c<x> <outfile> : same as previous option but for the complementary alignment, except ");
	System.out.println("                     for html output (i.e. -ch do not exist)");
	System.out.println("   -o<y> <outfile> : converts the trimmed alignment in Amino Acid- (-oaa), DNA- (-odna),");
	System.out.println("                     codon- (-oco) or RY- (-ory) coding sequences (can be combined with ");
	System.out.println("                     -o<x> options; default: no conversion)");
	System.out.println("                     (see documentation for more details and more output options)");
	System.out.println("");
    }







}
