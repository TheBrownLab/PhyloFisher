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
 *  The <code>RyAlignment</code> class extends the <code>Alignment</code> class to store RY-coding sequence alignment.
 *  @version 1.0
 **/

public class RYcodingAlignment extends Alignment {

    private final char[] character = {'R','Y','N','X','?','-'};
    /*  R  A,G        Y  C,T        N,X,?  unknown        -  gap     */

    private final int R = 0; private final int Y = 1;
    private double[] freq;
    private double total;

    private TreeSet<Character> alphabet;
    private int i;
    private StringBuffer sb;

    public RYcodingAlignment() {
	super();
	alphabet = new TreeSet<Character>();
	i = -1; while ( ++i < character.length ) alphabet.add(new Character(character[i]));
    }

    private String Filter(String sequence) {
	sb = new StringBuffer(sequence);
	i = -1;
	while ( ++i < sb.length() ) 
	    if ( ! alphabet.contains( new Character(sb.charAt(i)) ) ) sb.setCharAt(i , new Character('X'));
	return sb.toString();
    }

    public boolean add(String sequence) {
	return super.add( Filter(sequence) );
    }

    public boolean add(String sequence , String label) {
	return super.add( Filter(sequence) , label );
    }

    public boolean add(int row , String sequence) {
	return super.add( row , Filter(sequence) );
    }

    public boolean add(int row , String sequence , String label) {
	return super.add( row , Filter(sequence) , label );
    }

    public boolean set(int row , String sequence) {
	return super.set( row , Filter(sequence) );
    }

    public boolean set(int row , String sequence , String label) {
	return super.set( row , Filter(sequence) , label );
    }

    public boolean setCharAt(int row , int col , char ch) {
	if ( ! alphabet.contains(new Character(ch)) ) return false;
	return super.setCharAt( row , col , ch );
    }

    public int getAlphabetSize() {
	return 2;
    }

    public double[] getFrequencies(int col) {
	freq = new double[getAlphabetSize()];  // [f_R , f_Y]
	Arrays.fill(freq , 0);
	total = 0;
	i = -1;
	while ( ++i < super.size() ) {
	    total++;
	    switch ( super.getSequence(i).charAt(col) ) {
	    case 'R': freq[R]++; break;
	    case 'Y': freq[Y]++; break;
	    case '?':
	    case 'N':
	    case 'X': freq[R] += 0.5; freq[Y] += 0.5; break;
	    default: total--; break;
	    }
	}
	freq[R] /= total; freq[Y] /= total;
	return freq;
    }
}

	

