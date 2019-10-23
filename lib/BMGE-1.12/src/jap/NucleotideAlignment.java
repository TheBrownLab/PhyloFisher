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
 *  The <code>NucleotideAlignment</code> class extends the <code>Alignment</code> class to store DNA sequence alignment.
 *  @version 1.0
 **/

public class NucleotideAlignment extends Alignment {

    private final char[] character = {'A','G','C','T','U','M','R','W','S','Y','K','B','D','H','V','N','X','?','-'};
    //                A Adenosine              G Guanine                C Cytosine               T Thymine / U  Uracil
    private final int A = 0; private final int G = 1; private final int C = 2; private final int T = 3;
    //                M A,C     R A,G     W A,T     S C,G     Y C,T     K G,T
    //                B C,G,T   D A,G,T   H A,C,T   V A,C,G
    //                N,X,? unknown       -  gap     
    private final double THD = 1.0 / 3.0;
    private final double PVAL = 0.1;

    private static double[] freq;
    private static double total;
    private static TreeSet<Character> alphabet;
    private static int i, j, b;
    private static char c;
    private static StringBuffer sb;
    private static NucleotideAlignment dnaAlig;

    /**
     *  Constructs an empty <code>NucleotideAlignment</code>.
     */
    public NucleotideAlignment() {
	super();
	alphabet = new TreeSet<Character>();
	i = character.length; while ( --i >= 0 ) alphabet.add(new Character(character[i]));
    }

    /**
     *  Returns the specified sequence but all non DNA characters replaced by 'X'.
     *  @param sequence the <code>String</code> sequence to be filtered 
     *  @return the filtered <code>String</code> sequence
     */
    public static String filter(String sequence) {
	sb = new StringBuffer(sequence);
	j = -1; while ( ++j < sb.length() ) if ( ! alphabet.contains(new Character(sb.charAt(j))) ) sb.setCharAt(j , new Character('X'));
	return sb.toString();
    }

    public boolean isRNA() {
	i = -1; while ( ++i < super.size() ) if ( super.getSequence(i).indexOf("U") != -1 ) return true;
	return false;
    }

    public boolean setCharAt(int row , int col , char ch) {
	if ( ! alphabet.contains(new Character(ch)) ) return false;
	return super.setCharAt( row , col , ch );
    }

    public int getAlphabetSize() {
	return 4;
    }

    public double[] getFrequencies(int col) {
	freq = new double[getAlphabetSize()];  // [f_A , f_G , f_C , f_T]
	Arrays.fill(freq , 0);
	total = 0; i = super.size();
	while ( --i >= 0 ) {
	    total++;
	    switch ( super.getSequence(i).charAt(col) ) {
	    case 'A': freq[A]++; break;
	    case 'G': freq[G]++; break;
	    case 'C': freq[C]++; break;
	    case 'T': 
	    case 'U': freq[T]++; break;
	    case 'M': freq[A] += 0.5; freq[C] += 0.5; break;
	    case 'R': freq[A] += 0.5; freq[G] += 0.5; break;
	    case 'W': freq[A] += 0.5; freq[T] += 0.5; break;
	    case 'S': freq[C] += 0.5; freq[G] += 0.5; break;
	    case 'Y': freq[C] += 0.5; freq[T] += 0.5; break;
	    case 'K': freq[G] += 0.5; freq[T] += 0.5; break;
	    case 'B': freq[C] += THD; freq[G] += THD; freq[T] += THD; break;
	    case 'D': freq[A] += THD; freq[G] += THD; freq[T] += THD; break;
	    case 'H': freq[A] += THD; freq[C] += THD; freq[T] += THD; break;
	    case 'V': freq[A] += THD; freq[C] += THD; freq[G] += THD; break;
	    case 'N':
	    case 'X': freq[A] += 0.25; freq[G] += 0.25; freq[C] += 0.25; freq[T] += 0.25; break;
	    default: total--; break;
	    }
	}
	freq[A] /= total; freq[G] /= total; freq[C] /= total; freq[T] /= total;
	return freq;
    }

    public NucleotideAlignment toNucleotideAlignment() {
	dnaAlig = new NucleotideAlignment();
	i = -1; while ( ++i < super.size() ) dnaAlig.add(super.getSequence(i) , super.getLabel(i));
	return dnaAlig;
    }

    public NucleotideAlignment toRYcodingAlignment() {
	dnaAlig = new NucleotideAlignment();
	i = -1;
	while ( ++i < super.size() ) {
	    sb = new StringBuffer(super.getSequence(i));
	    j = -1;
	    while ( ++j < super.length() ) {
		switch(sb.charAt(j)) {
		case 'A':
		case 'G': sb.setCharAt(j , 'R'); break;
		case 'C':
		case 'T': sb.setCharAt(j , 'Y'); break;
		case '?': sb.setCharAt(j , '?'); break;
		case '-': sb.setCharAt(j , '-'); break;
		default: sb.setCharAt(j , 'X'); break;
		}
	    }
	    dnaAlig.add(sb.toString() , super.getLabel(i));
	}
	return dnaAlig;
    }






}

	

