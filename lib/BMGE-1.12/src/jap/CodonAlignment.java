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
 *  The <code>CodonAlignment</code> class extends the <code>Alignment</code> class to store codon sequence alignment.
 *  @version 1.0
 **/

public class CodonAlignment extends Alignment {

    private final String[] codonCharacter = {/*A*/ "GCT","GCC","GCA","GCG", "GCU", "GCX", "GCM","GCR","GCW","GCS","GCY","GCK","GCB","GCD","GCH","GCV",
					     /*R*/ "CGT","CGC","CGA","CGG","AGA","AGG", "CGU", "CGX","AGR","MGX", "CGM","CGR","CGW","CGS","CGY","CGK","CGB","CGD","CGH","CGV",
					     /*N*/ "AAT","AAC", "AAU", "AAY",
					     /*D*/ "GAT","GAC", "GAU", "GAY",
					     /*C*/ "TGT","TGC", "UGU","UGC", "TGY","UGY",
					     /*Q*/ "CAA","CAG", "CAR",
					     /*E*/ "GAA","GAG", "GAR",
					     /*G*/ "GGT","GGC","GGA","GGG", "GGU", "GGX", "GGM","GGR","GGW","GGS","GGY","GGK","GGB","GGD","GGH","GGV",
					     /*H*/ "CAT","CAC", "CAU", "CAY",
					     /*I*/ "ATT","ATC","ATA", "AUU","AUC","AUA", "ATH","AUH",
					     /*L*/ "TTA","TTG","CTT","CTC","CTA","CTG", "UUA","UUG","CUU","CUC","CUA","CUG", "TTR","CTX","YTX","UUR","CUX","YUX", "CTM","CTR","CTW","CTS","CTY","CTK","CTB","CTD","CTH","CTV", "CUM","CUR","CUW","CUS","CUY","CUK","CUB","CUD","CUH","CUV", "YTM","YTR","YTW","YTS","YTY","YTK","YTB","YTD","YTH","YTV", "YUM","YUR","YUW","YUS","YUY","YUK","YUB","YUD","YUH","YUV", 
					     /*K*/ "AAA","AAG", "AAR",
					     /*M*/ "ATG", "AUG",
					     /*F*/ "TTT","TTC", "UUU","UUC", "TTY","UUY",
					     /*P*/ "CCT","CCC","CCA","CCG", "CCU", "CCX", "CCM","CCR","CCW","CCS","CCY","CCK","CCB","CCD","CCH","CCV",
					     /*S*/ "TCT","TCC","TCA","TCG","AGT","AGC", "UCU","UCC","UCA","UCG","ACU", "TCX","UCX","AGY","WSX", "TCM","TCR","TCW","TCS","TCY","TCK","TCB","TCD","TCH","TCV", "UCM","UCR","UCW","UCS","UCY","UCK","UCB","UCD","UCH","UCV", "WSM","WSR","WSW","WSS","WSY","WSK","WSB","WSD","WSH","WSV", 
					     /*T*/ "ACT","ACC","ACA","ACG", "ACU", "ACX", "ACM","ACR","ACW","ACS","ACY","ACK","ACB","ACD","ACH","ACV",
					     /*W*/ "TGG", "UGG",
					     /*Y*/ "TAT","TAC", "UAU","UAC", "TAY","UAY",
					     /*V*/ "GTT","GTC","GTA","GTG", "GUU","GUC","GUA","GUG", "GTX","GUX", "GTM","GTR","GTW","GTS","GTY","GTK","GTB","GTD","GTH","GTV", "GUM","GUR","GUW","GUS","GUY","GUK","GUB","GUD","GUH","GUV",
					     /*X*/ "XXX",
					     /*?*/ "???",
					     /*-*/ "---"};
					     
    private final char[] aaCharacter = {'A','A','A','A', 'A', 'A', 'A','A','A','A','A','A','A','A','A','A',
					'R','R','R','R','R','R', 'R', 'R','R','R', 'R','R','R','R','R','R','R','R','R','R',
					'N','N', 'N', 'N',
					'D','D', 'D', 'D',
					'C','C', 'C','C', 'C','C',
					'Q','Q', 'Q',
					'E','E', 'E',
					'G','G','G','G', 'G', 'G', 'G','G','G','G','G','G','G','G','G','G',
					'H','H', 'H', 'H',
					'I','I','I', 'I','I','I', 'I','I',
					'L','L','L','L','L','L', 'L','L','L','L','L','L', 'L','L','L','L','L','L', 'L','L','L','L','L','L','L','L','L','L', 'L','L','L','L','L','L','L','L','L','L', 'L','L','L','L','L','L','L','L','L','L', 'L','L','L','L','L','L','L','L','L','L',
					'K','K', 'K',
					'M', 'M',
					'F','F', 'F','F', 'F','F',
					'P','P','P','P', 'P', 'P', 'P','P','P','P','P','P','P','P','P','P',
					'S','S','S','S','S','S', 'S','S','S','S','S', 'S','S','S','S', 'S','S','S','S','S','S','S','S','S','S', 'S','S','S','S','S','S','S','S','S','S', 'S','S','S','S','S','S','S','S','S','S',
					'T','T','T','T', 'T', 'T', 'T','T','T','T','T','T','T','T','T','T',
					'W', 'W',
					'Y','Y', 'Y','Y', 'Y','Y',
					'V','V','V','V', 'V','V','V','V', 'V','V', 'V','V','V','V','V','V','V','V','V','V', 'V','V','V','V','V','V','V','V','V','V',
					'X',
					'?',
					'-'};
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

    private NucleotideAlignment dnaAlignment;
    private static ArrayList<String> codonAlphabet;
    private static ArrayList<Character> aminoAcidAlphabet;

    private static int i, j, c, cc, cs;
    private static String codon;
    private static char ch;
    private static StringBuffer sb;
    private static double[] freq;
    private static double total;
    private static boolean rna;
    private static AminoAcidAlignment aaAlig;
    private static CodonAlignment coAlig;
    private static NucleotideAlignment dnaAlig;

    public CodonAlignment() {
	super();
	codonAlphabet = new ArrayList<String>();
	aminoAcidAlphabet = new ArrayList<Character>();
	i = -1; 
	while ( ++i < codonCharacter.length ) {
	    codonAlphabet.add(codonCharacter[i]);
	    aminoAcidAlphabet.add(new Character(aaCharacter[i]));
	}	    
	dnaAlignment = new NucleotideAlignment();
    }

    public static char toAminoAcidCharacter(String codon) {
	if ( (cc = codonAlphabet.indexOf(codon)) != -1 ) return aminoAcidAlphabet.get(cc).charValue();
	return 'X';
    }

    public static String toAminoAcidSequence(String codonSequence) {
	sb = new StringBuffer(codonSequence.substring(0 , codonSequence.length()/3));
	cs = 0; while ( cs < sb.length() ) sb.setCharAt(cs , toAminoAcidCharacter(codonSequence.substring(3*cs , 3*(++cs))));
	return sb.toString();
    }

    /**
     *  Returns the specified sequence but all non codon characters replaced by 'XXX'.
     *  @param sequence the <code>String</code> sequence to be filtered 
     *  @return the filtered <code>String</code> sequence
     */
    public static String filter(String sequence) {
	sb = new StringBuffer(sequence);
	i = 0; 
	while ( i < sb.length()/3 ) 
	    if ( ! codonAlphabet.contains(sequence.substring(3*i , 3*(++i))) ) sb = sb.replace(3*(--i) , 3*(++i) , "XXX");
	return sb.toString();
    }

    public boolean add(String sequence) {
	if ( sequence.length() % 3 == 0 ) 
	    if ( super.add(toAminoAcidSequence(sequence)) ) return dnaAlignment.add(sequence , super.getLabel(super.size()-1));
	return false;
    }

    public boolean add(String sequence , String label) {
	if ( sequence.length() % 3 == 0 ) 
	    if ( super.add(toAminoAcidSequence(sequence) , label) ) return dnaAlignment.add(sequence , label);
	return false;
    }

    public boolean add(int row , String sequence) {
	if ( sequence.length() % 3 == 0 )
	    if ( super.add(row , toAminoAcidSequence(sequence)) ) return dnaAlignment.add(row , sequence , super.getLabel(super.size()-1));
	return false;
    }

    public boolean add(int row , String sequence , String label) {
	if ( sequence.length() % 3 == 0 )
	    if ( super.add(row , toAminoAcidSequence(sequence) , label) ) return dnaAlignment.add(row , sequence , label);
	return false;
    }

    public boolean set(int row , String sequence) {
	if ( sequence.length() % 3 == 0 ) 
	    if ( super.set(row , toAminoAcidSequence(sequence)) ) return dnaAlignment.set(row , sequence);
	return false;
    }

    public boolean set(int row , String sequence , String label) {
	if ( sequence.length() % 3 == 0 )
	    if ( super.set(row , toAminoAcidSequence(sequence) , label) ) return dnaAlignment.set(row , sequence , label);
	return false;
    }

    public boolean removeColumn(int col) {
	j = 3 * col;
	if ( super.removeColumn(col) ) return dnaAlignment.removeColumn(j) && dnaAlignment.removeColumn(j) && dnaAlignment.removeColumn(j);  
	return false;
    }

    public String getSequence(int row) {
	return dnaAlignment.getSequence(row);
    }

    public boolean setCharAt(int row , int col , char ch) {
	if ( ! aminoAcidAlphabet.contains(new Character(ch)) ) return false;
	if ( super.setCharAt(row , col , ch) ) {
	    rna = dnaAlignment.isRNA();
	    switch (ch) {
	    case 'A': codon = "GCX"; break;
	    case 'R': codon = "MGX"; break;
	    case 'N': codon = "AAY"; break;
	    case 'D': codon = "GAY"; break;
	    case 'C': if ( rna ) codon = "UGY"; else codon = "TGY"; break;
	    case 'Q': codon = "CAR"; break;
	    case 'E': codon = "GAR"; break;
	    case 'G': codon = "GGX"; break;
	    case 'H': codon = "CAY"; break;
	    case 'I': if ( rna ) codon = "AUH"; else codon = "ATH"; break;
	    case 'L': if ( rna ) codon = "YUX"; else codon = "YTX"; break;
	    case 'K': codon = "AAR"; break;
	    case 'M': if ( rna ) codon = "AUG"; else codon = "ATG"; break;
	    case 'F': if ( rna ) codon = "UUY"; else codon = "TTY"; break;
	    case 'P': codon = "CCX"; break;
	    case 'S': codon = "WSX"; break;
	    case 'T': codon = "ACX"; break;
	    case 'W': if ( rna ) codon = "UGG"; else codon = "TGG"; break;
	    case 'Y': if ( rna ) codon = "UAY"; else codon = "TAY"; break;
	    case 'V': if ( rna ) codon = "GUX"; else codon = "GTX"; break;
	    case 'X': codon = "XXX"; break;
	    case '-': codon = "---"; break;
	    default: codon="???"; break;
	    }
	    return dnaAlignment.setCharAt(row , 3*col , codon.charAt(0)) 
		&& dnaAlignment.setCharAt(row , 3*col+1 , codon.charAt(1))
		&& dnaAlignment.setCharAt(row , 3*col+2 , codon.charAt(2));
	}
	return false;
    }

    public String getCodonAt(int row , int col) {
	if ( (super.size() > row) && (row >= 0) && (super.length() > col) && (col >= 0) ) 
	    return dnaAlignment.getSequence(row).substring(3*col , 3*(col+1));
	return null;
    }

    public boolean setCodonAt(int row , int col , String codon) {
	return codonAlphabet.contains(codon) 
	    && super.setCharAt(row , col , toAminoAcidCharacter(codon))
	    && dnaAlignment.setCharAt(row , 3*col , codon.charAt(0))
	    && dnaAlignment.setCharAt(row , 3*col+1 , codon.charAt(1))
	    && dnaAlignment.setCharAt(row , 3*col+2 , codon.charAt(2));
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
	return toNucleotideAlignment(7);
    }

    public NucleotideAlignment toNucleotideAlignment(int codonCode) {
	dnaAlig = new NucleotideAlignment();
	if ( (codonCode <= 0) || (codonCode > 7) ) return dnaAlig; // no codon position
	i = -1; while ( ++i < super.size() ) dnaAlig.add(dnaAlignment.getSequence(i) , super.getLabel(i));
	if ( codonCode == 7 ) return dnaAlig; // codon positions 123
	switch (codonCode) {
	case 4:               // codon position 3
	case 6: c = 0; break; // codon positions 23 => remove position 1
	case 5: c = 1; break; // codon positions 13 => remove position 2
	case 1:               // codon position 1
	case 2:               // codon position 2
	case 3: c = 2; break; // codon positions 12 => remove position 3
	}
	j = super.length();
	while ( --j >= 0 ) dnaAlig.removeColumn(3*j + c);
	if ( (codonCode == 3) || (codonCode == 5) || (codonCode == 6) ) return dnaAlig;
	switch (codonCode) {
	case 4: c = 1; break; // codon position 3 => remove position 2
	case 2: c = 0; break; // codon position 2 => remove position 1
	case 1: c = 1; break; // codon position 1 => remove position 2
	}
	j = dnaAlig.length();
	while ( --j >= 0 ) dnaAlig.removeColumn(2*j + c);
	return dnaAlig;
    }

    public CodonAlignment toCodonAlignment() {
	coAlig = new CodonAlignment();
	i = -1; while ( ++i < super.size() ) coAlig.add(dnaAlignment.getSequence(i) , super.getLabel(i));
	return coAlig;
    }

    public NucleotideAlignment toRYcodingAlignment() {
	dnaAlig = new NucleotideAlignment();
	i = -1; 
	while ( ++i < super.size() ) {
	    dnaAlig.add(dnaAlignment.getSequence(i) , super.getLabel(i));
	    //System.out.println( dnaAlig.getSequence(i) );
	    j = -1;
	    while ( ++j < dnaAlig.length() ) {
		ch = dnaAlig.charAt(i , j);
		switch(ch) {
		case 'A':
		case 'G': dnaAlig.setCharAt(i , j , 'R'); break;
		case 'C':
		case 'T': dnaAlig.setCharAt(i , j , 'Y'); break;
		case '?': dnaAlig.setCharAt(i , j , '?'); break;
		case '-': dnaAlig.setCharAt(i , j , '-'); break;
		default: dnaAlig.setCharAt(i , j , 'X'); break;
		}
	    }
	    //System.out.println( dnaAlig.getSequence(i) );
	}
	return dnaAlig;
    }


}

	

