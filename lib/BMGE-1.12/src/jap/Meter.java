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


/**
 *  The <code>Meter</code> class can be used to store every type of sequence alignment.
 *  @version 1.0
 **/
public class Meter {

    private double[][] ctgc; 

    public Meter() {
    	ctgc = new double[0][0];
    }
    
    public Meter(String sequence1 , String sequence2) {
	this();
    }
    
    public void increment( char charState1 , char charState2 ) {}
    public void increment( String codonState1 , String codonState2 ) {}
    public void increment( int row , int col ) {}
    public void decrement( char charState1 , char charState2 ) {}
    public void decrement( String codonState1 , String codonState2 ) {}
    public void decrement( int row , int col ) {}
    public double get(int row , int col) { return -1.0; }
    public double get( char charState1 , char charState2 ) { return -1.0; }
    public void set(int row , int col , double x) {}
    public int size() { return 0; }
    public Meter copy() { return new Meter(); }
    public double getBowkerSymmetryTestPvalue() { return 0; }
    public double getEvansHoenigSymmetryTestPvalue() { return 0; }
    public double getStuartMarginalSymmetryTestValue() { return 0; }
    public double getStuartMarginalSymmetryTestPvalue() { return 0; }
    public char getChar(int charIndex) { return '-'; }


}



