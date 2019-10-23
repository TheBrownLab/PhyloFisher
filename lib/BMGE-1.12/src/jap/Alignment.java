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
 *  The <code>Alignment</code> class can be used to store every type of sequence alignment.
 *  @version 1.0
 **/

public class Alignment {

    private ArrayList<StringBuffer> sequence;
    private ArrayList<String> label;

    private ArrayList<Character> alphabet;
    private double[] freq;
    private double total;
    private int i, j, c;
    private char ch;
    private String line;
    private ArrayList<String> charState;
    private ArrayList<Integer> charStateCount;

    public Alignment() {
	sequence = new ArrayList<StringBuffer>(0);
	label = new ArrayList<String>(0);
	alphabet = new ArrayList<Character>(0);
    }

    public int size() {
	return sequence.size();
    }

    public int length() {
	return ( this.size() == 0 ) ? -1 : sequence.get(0).length();
	/*
	  if ( this.size() == 0 ) return -1;
	  return sequence.get(0).length();
	*/
    }

    /*private String Filter(String sequence) {
      sb = new StringBuffer( sequence );
      i = -1;
      while ( ++i < sb.length() ) {
      ch = sb.charAt(i);
      if ( ch == ' ' ) {
      sb.deleteCharAt(i--);
      continue;
      }
      if ( (ch != '-') && (! alphabet.contains(new Character(ch))) ) alphabet.add( new Character(ch) );
      }
      return sb.toString();
      }*/

    public boolean add(String sequence) {
	i = 0; line = "label" + i; while ( ! this.label.contains(line) ) line = "label" + (++i);
	return this.add(sequence , line);
    }

    public boolean add(String sequence , String label) {
	if ( (this.size() == 0) || (sequence.length() == this.length()) ) {
	    this.sequence.add(new StringBuffer(sequence.toUpperCase())); this.label.add(label); return true;
	}
	return false;
    }

    public boolean add(int row , String sequence) {
	i = 0; line = "label" + i;
	while ( ! this.label.contains(line) ) line = "label" + (++i);
	return this.add(row , sequence , line);
    }

    public boolean add(int row , String sequence , String label) {
	if ( (this.size() > row) && (row >= 0) && (sequence.length() == this.length()) ) {
	    this.sequence.add(row , new StringBuffer(sequence.toUpperCase())); this.label.add(row , label); return true;
	}
	return false;
    }

    public boolean set(int row , String sequence) {
	return this.set(row , sequence , this.label.get(row));
    }

    public boolean set(int row , String sequence , String label) {
	if ( (this.size() > row) && (row >= 0) && (sequence.length() == this.length()) ) {
	    this.sequence.set(row , new StringBuffer(sequence.toUpperCase())); this.label.set(row , label); return true;
	}
	return false;
    }

    public boolean removeRow(int row) {
	if ( (this.size() > row) && (row >= 0) ) {sequence.remove(row); label.remove(row); return true;}
	return false;
    }
	
    public boolean removeColumn(int col) {
	if ( (this.length() > col) && (col >= 0) ) {
	    i = -1; while ( ++i < this.size() ) sequence.set(i , sequence.get(i).deleteCharAt(col));
	    return true;
	}
	return false;
    }
	
    public String getSequence(int row) {
	return ( (this.size() > row) && (row >= 0) ) ? sequence.get(row).toString() : null;
	/*
	  if ( (this.size() > row) && (row >= 0) ) return sequence.get(row).toString();
	  return null;
	*/
    }
	    
    public String getLabel(int row) {
	return ( (this.size() > row) && (row >= 0) ) ? label.get(row) : null;
	/*
	  if ( (this.size() > row) && (row >= 0) ) return label.get(row);
	  return null;
	*/
    }
	    
    public char charAt(int row , int col) {
	return ( (this.size() > row) && (row >= 0) && (this.length() > col) && (col >= 0) ) 
	    ? sequence.get(row).charAt(col) 
	    : (char) 0;
	/*
	  if ( (this.size() > row) && (row >= 0) && (this.length() > col) && (col >= 0) ) return sequence.get(row).charAt(col);
	  return (char)0;
	*/
    }

    public String getCharAt(int row , int col) {
	return Character.toString(this.charAt(row , col));
    }

    public String getCodonAt(int row , int col) {
	return null;
    }

    public boolean setCodonAt(int row , int col , String cod) {
	return false;
    }

    public boolean setCharAt(int row , int col , char ch) {
	if ( (this.size() > row) && (row >= 0) && (this.length() > col) && (col >= 0) ) {
	    sequence.get(row).setCharAt(col , ch); return true;
	}
	return false;
    }

    public int getAlphabetSize() {
	return alphabet.size();
    }

    public double[] getFrequencies(int col) {
	freq = new double[this.getAlphabetSize()];
	Arrays.fill(freq , 0);
	return freq;
    }

    public double getColGapRate(int col) {
	if ( (this.length() > col) && (col >= 0) ) {
	    c = 0; i = this.size();
	    while ( --i >= 0 ) c += ( sequence.get(i).charAt(col) == '-' ) ? 1 : 0;
	    return ((double)c)/this.size();
	}
	return 0;
	/*
	  if ( (this.length() > col) && (col >= 0) ) {
	  c = 0; total = 0; i = -1;
	  while ( ++i < this.size() ) {total++; if ( sequence.get(i).charAt(col) == '-' ) c++;}
	  return ((double)c)/total;
	  }
	  return 0;
	*/
    }	

    public double getRowGapRate(int row) {
	if ( (this.size() > row) && (row >= 0) ) {
	    c = 0; j = this.length();
	    while ( --j >= 0 ) c += ( sequence.get(row).charAt(j) == '-' ) ? 1 : 0;
	    return ((double)c)/this.length();
	}
	return 0;
	/*
	  if ( (this.size() > row) && (row >= 0) ) {
	  c = 0; total = 0; j = -1;
	  while ( ++j < this.length() ) {total++; if ( sequence.get(row).charAt(j) == '-' ) c++;}
	  return ((double)c)/total;
	  }
	  return 0;
	*/
    }	

    public ArrayList<String> getMajorityCharacter(int col) {
	charState = new ArrayList<String>(0);
	if ( (this.length() > col) && (col >= 0) ) {
	    charStateCount = new ArrayList<Integer>(0);
	    total = 0; i = -1;
	    while ( ++i < size() ) {
		line = this.getCharAt(i , col);
		if ( line.equals("-") || line.equals("X") || line.equals("?") ) continue;
		if ( ! charState.contains(line) ) {
		    charState.add(line); charStateCount.add(new Integer(1));
		    if ( total < 1 ) total = 1;
		}
		else {
		    j = charState.indexOf(line); c = charStateCount.get(j).intValue();
		    charStateCount.set(j , new Integer(++c));
		    if ( total < c ) total = c;
		}
	    }
	    c = -1;
	    while ( ++c < charState.size() ) 
		if ( charStateCount.get(c).intValue() < total ) {charState.remove(c); charStateCount.remove(c--);}
	}
	return charState;
    }
	

    public AminoAcidAlignment toAminoAcidAlignment() {
	return new AminoAcidAlignment();
    }

    public NucleotideAlignment toNucleotideAlignment() {
	return new NucleotideAlignment();
    }

    public NucleotideAlignment toNucleotideAlignment(int codonCode) {
	return new NucleotideAlignment();
    }

    public CodonAlignment toCodonAlignment() {
	return new CodonAlignment();
    }

    public NucleotideAlignment toRYcodingAlignment() {
	return new NucleotideAlignment();
    }

    public BitSet stationaryTrimming( ArrayList<Double> entropy , BitSet keepCol , BitSet keepRow ) {
	return new BitSet();
    }


}

	

