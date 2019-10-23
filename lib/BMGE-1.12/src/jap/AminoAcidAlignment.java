/*
  BMGE (Block Mapping and Gathering with Entropy): selection of phylogenetic informative regions from multiple sequence alignments
  Copyright (C) 2010  Alexis Criscuolo 

  This file is part of BMGE.

  BMGE is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  BMGE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Contact: 
  Unité de Biologie Moléculaire du Gène chez les Extrêmophiles
  Département de Microbiologie
  INSTITUT  PASTEUR
  25 rue du Dr Roux - 75015 Paris  (France)

  alexis.criscuolo@pasteur.fr
*/

package jap;
import java.util.*;

/**
 *  The <code>AminoAcidAlignment</code> class extends the <code>Alignment</code> class to store amino acid sequence alignment.
 *  @version 1.0
 **/

public class AminoAcidAlignment extends Alignment {

    private final char[] character = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X','?','-'};
    //                A  Alanine                R  Arginine               N  Asparagine             D  Aspartic acid
    private final int A = 0;  private final int R = 1;  private final int N = 2;  private final int D = 3;  
    //                C  Cysteine               Q  Glutamine              E  Glutamic acid          G  Glycine
    private final int C = 4;  private final int Q = 5;  private final int E = 6;  private final int G = 7;  
    //                H  Histidine              I  Isoleucine             L  Leucine                K  Lysine
    private final int H = 8;  private final int I = 9;  private final int L = 10; private final int K = 11; 
    //                M  Methionine             F  Phenylalanine          P  Proline                S  Serine
    private final int M = 12; private final int F = 13; private final int P = 14; private final int S = 15; 
    //                T  Threonine              W  Tryptophan             Y  Tyrosine               V  Valine
    private final int T = 16; private final int W = 17; private final int Y = 18; private final int V = 19;
    //                B  D,N                    Z  Q,E                    X,? unknown               -  gap                      
    private static double[] freq;
    private static double total;
    private static TreeSet<Character> alphabet;
    private static int i, j, c;
    private static StringBuffer sb;
    private static String codon;
    private static AminoAcidAlignment aaAlig;
    private static CodonAlignment coAlig, coAlig_;
    private static NucleotideAlignment dnaAlig;

    /**
     *  Constructs ans empty <code>AminoAcidAlignment</code>.
     */
    public AminoAcidAlignment() {
	super();
	alphabet = new TreeSet<Character>();
	c = -1; while ( ++c < character.length ) alphabet.add(new Character(character[c]));
    }

    /**
     *  Returns the specified sequence but all non amino acid characters replaced by 'X'.
     *  @param sequence the <code>String</code> sequence to be filtered 
     *  @return the filtered <code>String</code> sequence
     */
    public static String filter(String sequence) {
	sb = new StringBuffer(sequence);
	j = -1; while ( ++j < sb.length() ) if ( ! alphabet.contains(new Character(sb.charAt(j))) ) sb.setCharAt(j , new Character('X'));
	return sb.toString();
    }

    public boolean setCharAt(int row , int col , char ch) {
	if ( ! alphabet.contains(new Character(ch)) ) return false;
	return super.setCharAt(row, col, ch);
    }

    public int getAlphabetSize() {
	return 20;
    }

    public double[] getFrequencies(int col) {
	freq = new double[getAlphabetSize()];  // [fA,fR,fN,fD,fC,fQ,fE,fG,fH,fI,fL,fK,fM,fF,fP,fS,fT,fW,fY,fV]
	Arrays.fill(freq , 0);
	total = 0; i = -1;
	while ( ++i < super.size() ) {
	    total++;
	    switch ( super.getSequence(i).charAt(col) ) {
	    case 'A': freq[A]++; break;
	    case 'R': freq[R]++; break;
	    case 'N': freq[N]++; break;
	    case 'D': freq[D]++; break;
	    case 'C': freq[C]++; break;
	    case 'Q': freq[Q]++; break;
	    case 'E': freq[E]++; break;
	    case 'G': freq[G]++; break;
	    case 'H': freq[H]++; break;
	    case 'I': freq[I]++; break;
	    case 'L': freq[L]++; break;
	    case 'K': freq[K]++; break;
	    case 'M': freq[M]++; break;
	    case 'F': freq[F]++; break;
	    case 'P': freq[P]++; break;
	    case 'S': freq[S]++; break;
	    case 'T': freq[T]++; break;
	    case 'W': freq[W]++; break;
	    case 'Y': freq[Y]++; break;
	    case 'V': freq[V]++; break;
	    case 'B': freq[D] += 0.5; freq[N] += 0.5; break;
	    case 'Z': freq[Q] += 0.5; freq[E] += 0.5; break;     
	    case 'X': c = -1; while ( ++c < this.getAlphabetSize() ) freq[c] += 0.05; break;
	    default: total--; break;
	    }
	}
	c = -1; while ( ++c < this.getAlphabetSize() ) freq[c] /= total;
	return freq;
    }

    public AminoAcidAlignment toAminoAcidAlignment() {
	aaAlig = new AminoAcidAlignment();
	i = -1; while ( ++i < super.size() ) aaAlig.add(super.getSequence(i) , super.getLabel(i));
	return aaAlig;
    }
		
    public NucleotideAlignment toNucleotideAlignment() {
	coAlig_ = this.toCodonAlignment();
	return coAlig_.toNucleotideAlignment(7);
    }

    public CodonAlignment toCodonAlignment() {
	coAlig = new CodonAlignment();
	i = -1;
	while ( ++i < super.size() ) {
	    sb = new StringBuffer("");
	    j = -1;
	    while ( ++j < super.length() ) {
		switch(super.charAt(i , j)) {
		case 'A': codon = "GCX"; break;
		case 'R': codon = "MGX"; break;
		case 'N': codon = "AAY"; break;
		case 'D': codon = "GAY"; break;
		case 'C': codon = "TGY"; break;
		case 'Q': codon = "CAR"; break;
		case 'E': codon = "GAR"; break;
		case 'G': codon = "GGX"; break;
		case 'H': codon = "CAY"; break;
		case 'I': codon = "ATH"; break;
		case 'L': codon = "YTX"; break;
		case 'K': codon = "AAR"; break;
		case 'M': codon = "ATG"; break;
		case 'F': codon = "TTY"; break;
		case 'P': codon = "CCX"; break;
		case 'S': codon = "WSX"; break;
		case 'T': codon = "ACX"; break;
		case 'W': codon = "TGG"; break;
		case 'Y': codon = "TAY"; break;
		case 'V': codon = "GTX"; break;
		case 'X': codon = "XXX"; break;
		case '-': codon = "---"; break;
		default: codon="???"; break;
		}
		sb = sb.append(codon);
	    }
	    coAlig.add(sb.toString() , super.getLabel(i));
	}
	return coAlig;
    }

    public NucleotideAlignment toRYcodingAlignment() {
	dnaAlig = new NucleotideAlignment();
	i = -1;
	while ( ++i < super.size() ) {
	    sb = new StringBuffer("");
	    j = -1;
	    while ( ++j < super.length() ) {
		switch(super.charAt(i , j)) {
		case 'A': codon = "RYX"; break;
		case 'R': codon = "XRX"; break;
		case 'N': codon = "RRY"; break;
		case 'D': codon = "RRY"; break;
		case 'C': codon = "YRY"; break;
		case 'Q': codon = "YRR"; break;
		case 'E': codon = "RRR"; break;
		case 'G': codon = "RRX"; break;
		case 'H': codon = "YRY"; break;
		case 'I': codon = "RYX"; break;
		case 'L': codon = "YYX"; break;
		case 'K': codon = "RRR"; break;
		case 'M': codon = "RYR"; break;
		case 'F': codon = "YYY"; break;
		case 'P': codon = "YYX"; break;
		case 'S': codon = "XXX"; break;
		case 'T': codon = "RYX"; break;
		case 'W': codon = "YRR"; break;
		case 'Y': codon = "YRY"; break;
		case 'V': codon = "RYX"; break;
		case 'X': codon = "XXX"; break;
		case '-': codon = "---"; break;
		default: codon="???"; break;
		}
		sb = sb.append(codon);
	    }
	    dnaAlig.add(sb.toString() , super.getLabel(i));
	}
	return dnaAlig;
    }



	
}

	

