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
import Jama.*;

public class CodonMeter extends Meter {

    private final int SIZE = 4; 

    private NucleotideMeter nm;
    private static int i, j, s;
    private static String codon1, codon2;
    private static CodonMeter com_copy;

    public CodonMeter() {
	nm = new NucleotideMeter();
    }

    public CodonMeter(String codonSequence1 , String codonSequence2) {
	nm = new NucleotideMeter(codonSequence1 , codonSequence2);
    }   
		
    public void increment( String codonState1 , String codonState2 ) {
	if ( (codonState1.length() == 3) && (codonState2.length() == 3) ) {
	    nm.increment(codonState1.charAt(0) , codonState2.charAt(0));
	    nm.increment(codonState1.charAt(1) , codonState2.charAt(1));
	    nm.increment(codonState1.charAt(2) , codonState2.charAt(2));
	}
    }

    public void increment( char charState1 , char charState2 ) {
	codon1 = ""; 
	switch (charState1) {
	case 'A': codon1 = "GCX"; break;
	case 'R': codon1 = "MGX"; break;
	case 'N': codon1 = "AAY"; break;
	case 'D': codon1 = "GAY"; break;
	case 'C': codon1 = "TGY"; break;
	case 'Q': codon1 = "CAR"; break;
	case 'E': codon1 = "GAR"; break;
	case 'G': codon1 = "GGX"; break;
	case 'H': codon1 = "CAY"; break;
	case 'I': codon1 = "ATH"; break;
	case 'L': codon1 = "YTX"; break;
	case 'K': codon1 = "AAR"; break;
	case 'M': codon1 = "ATG"; break;
	case 'F': codon1 = "TTY"; break;
	case 'P': codon1 = "CCX"; break;
	case 'S': codon1 = "WSX"; break;
	case 'T': codon1 = "ACX"; break;
	case 'W': codon1 = "TGG"; break;
	case 'Y': codon1 = "TAY"; break;
	case 'V': codon1 = "GTX"; break;
	case 'X': codon1 = "XXX"; break;
	case '-': codon1 = "---"; break;
	default:  codon1 = "???"; break;
	}
	codon2 = ""; 
	switch (charState2) {
	case 'A': codon2 = "GCX"; break;
	case 'R': codon2 = "MGX"; break;
	case 'N': codon2 = "AAY"; break;
	case 'D': codon2 = "GAY"; break;
	case 'C': codon2 = "TGY"; break;
	case 'Q': codon2 = "CAR"; break;
	case 'E': codon2 = "GAR"; break;
	case 'G': codon2 = "GGX"; break;
	case 'H': codon2 = "CAY"; break;
	case 'I': codon2 = "ATH"; break;
	case 'L': codon2 = "YTX"; break;
	case 'K': codon2 = "AAR"; break;
	case 'M': codon2 = "ATG"; break;
	case 'F': codon2 = "TTY"; break;
	case 'P': codon2 = "CCX"; break;
	case 'S': codon2 = "WSX"; break;
	case 'T': codon2 = "ACX"; break;
	case 'W': codon2 = "TGG"; break;
	case 'Y': codon2 = "TAY"; break;
	case 'V': codon2 = "GTX"; break;
	case 'X': codon2 = "XXX"; break;
	case '-': codon2 = "---"; break;
	default:  codon2 = "???"; break;
	}
	nm.increment(codon1 , codon2);
    }

    public void decrement( String codonState1 , String codonState2 ) {
	if ( (codonState1.length() == 3) && (codonState2.length() == 3) ) {
	    nm.decrement(codonState1.charAt(0) , codonState2.charAt(0));
	    nm.decrement(codonState1.charAt(1) , codonState2.charAt(1));
	    nm.decrement(codonState1.charAt(2) , codonState2.charAt(2));
	}
    }

    public void decrement( char charState1 , char charState2 ) {
	codon1 = ""; 
	switch (charState1) {
	case 'A': codon1 = "GCX"; break;
	case 'R': codon1 = "MGX"; break;
	case 'N': codon1 = "AAY"; break;
	case 'D': codon1 = "GAY"; break;
	case 'C': codon1 = "TGY"; break;
	case 'Q': codon1 = "CAR"; break;
	case 'E': codon1 = "GAR"; break;
	case 'G': codon1 = "GGX"; break;
	case 'H': codon1 = "CAY"; break;
	case 'I': codon1 = "ATH"; break;
	case 'L': codon1 = "YTX"; break;
	case 'K': codon1 = "AAR"; break;
	case 'M': codon1 = "ATG"; break;
	case 'F': codon1 = "TTY"; break;
	case 'P': codon1 = "CCX"; break;
	case 'S': codon1 = "WSX"; break;
	case 'T': codon1 = "ACX"; break;
	case 'W': codon1 = "TGG"; break;
	case 'Y': codon1 = "TAY"; break;
	case 'V': codon1 = "GTX"; break;
	case 'X': codon1 = "XXX"; break;
	case '-': codon1 = "---"; break;
	default:  codon1 = "???"; break;
	}
	codon2 = ""; 
	switch (charState2) {
	case 'A': codon2 = "GCX"; break;
	case 'R': codon2 = "MGX"; break;
	case 'N': codon2 = "AAY"; break;
	case 'D': codon2 = "GAY"; break;
	case 'C': codon2 = "TGY"; break;
	case 'Q': codon2 = "CAR"; break;
	case 'E': codon2 = "GAR"; break;
	case 'G': codon2 = "GGX"; break;
	case 'H': codon2 = "CAY"; break;
	case 'I': codon2 = "ATH"; break;
	case 'L': codon2 = "YTX"; break;
	case 'K': codon2 = "AAR"; break;
	case 'M': codon2 = "ATG"; break;
	case 'F': codon2 = "TTY"; break;
	case 'P': codon2 = "CCX"; break;
	case 'S': codon2 = "WSX"; break;
	case 'T': codon2 = "ACX"; break;
	case 'W': codon2 = "TGG"; break;
	case 'Y': codon2 = "TAY"; break;
	case 'V': codon2 = "GTX"; break;
	case 'X': codon2 = "XXX"; break;
	case '-': codon2 = "---"; break;
	default:  codon2 = "???"; break;
	}
	nm.decrement(codon1 , codon2);
    }

    public double get(int row , int col) {
	return nm.get(row , col);
    }

    public void set(int row , int col , double x) {
	nm.set(row , col , x);
    }

    public int size() { 
	return SIZE; 
    }

    public CodonMeter copy() {
	com_copy = new CodonMeter();
	i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) com_copy.set(i , j , nm.get(i , j)); }
	return com_copy;
    }

    public double getBowkerSymmetryTestPvalue() {
	return nm.getBowkerSymmetryTestPvalue();
    }

    public double getEvansHoenigSymmetryTestPvalue() {
	return nm.getEvansHoenigSymmetryTestPvalue();
    }

    public double getStuartMarginalSymmetryTestPvalue() {
	return nm.getStuartMarginalSymmetryTestPvalue();
    }

}
