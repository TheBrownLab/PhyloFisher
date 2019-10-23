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

public class AminoAcidMeter extends Meter {

    private final char[] character = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X','?','-'};
    //                      A  Alanine                      R  Arginine                     N  Asparagine           
    public final static int A = 0;  public final static int R = 1;  public final static int N = 2;  
    //                      D  Aspartic acid                C  Cysteine                     Q  Glutamine 
    public final static int D = 3;  public final static int C = 4;  public final static int Q = 5;  
    //                      E  Glutamic acid                G  Glycine                      H  Histidine 
    public final static int E = 6;  public final static int G = 7;  public final static int H = 8;  
    //                      I  Isoleucine                   L  Leucine                      K  Lysine
    public final static int I = 9;  public final static int L = 10; public final static int K = 11; 
    //                      M  Methionine                   F  Phenylalanine                P  Proline 
    public final static int M = 12; public final static int F = 13; public final static int P = 14; 
    //                      S  Serine                       T  Threonine                    W  Tryptophan
    public final static int S = 15; public final static int T = 16; public final static int W = 17; 
    //                      Y  Tyrosine                     V  Valine
    public final static int Y = 18; public final static int V = 19;
    //               B  D,N                   Z  Q,E                   X,? unknown              -  gap                      
    public final int B = 20; public final int Z = 21;
    private final int SIZE = 20;

    private double[][] ctgc; 
    private static int i, j, s, m, index1, index2, dof, p;
    private static double chi2, up, down, sum_r, sum_c;
    private static char c1, c2;
    private static AminoAcidMeter aam_copy;
    private static Matrix sumRow, sumCol, vector_U, vector_UT, matrix_V, matrix;
    private static int[] index;

    public AminoAcidMeter() {
	ctgc = new double[SIZE][SIZE];
	i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) ctgc[i][j] = 0; }
    }

    public AminoAcidMeter(String aminoAcidSequence1 , String aminoAcidSequence2) {
	this();
	m = Math.min(aminoAcidSequence1.length() , aminoAcidSequence2.length());
	s = -1;	while ( ++s < m ) increment( aminoAcidSequence1.charAt(s) , aminoAcidSequence2.charAt(s) );
	/*i = -1;
	  while ( ++i < SIZE ) {
	  System.out.print(" " + character[i] + "  ");
	  j = -1; while ( ++j < SIZE ) System.out.print( ctgc[i][j] + " " );
	  System.out.println("");
	  }*/
    }   
		
    public void increment( char charState1 , char charState2 ) {
	index1 = ( charState1 == 'A' ) ? A : ( charState1 == 'R' ) ? R : ( charState1 == 'N' ) ? N 
	    : ( charState1 == 'D' ) ? D : ( charState1 == 'C' ) ? C : ( charState1 == 'Q' ) ? Q 
	    : ( charState1 == 'E' ) ? E : ( charState1 == 'G' ) ? G : ( charState1 == 'H' ) ? H 
	    : ( charState1 == 'I' ) ? I : ( charState1 == 'L' ) ? L : ( charState1 == 'K' ) ? K
	    : ( charState1 == 'M' ) ? M : ( charState1 == 'F' ) ? F : ( charState1 == 'P' ) ? P 
	    : ( charState1 == 'S' ) ? S : ( charState1 == 'T' ) ? T : ( charState1 == 'W' ) ? W 
	    : ( charState1 == 'Y' ) ? Y : ( charState1 == 'V' ) ? V : -1;
	index2 = ( charState2 == 'A' ) ? A : ( charState2 == 'R' ) ? R : ( charState2 == 'N' ) ? N 
	    : ( charState2 == 'D' ) ? D : ( charState2 == 'C' ) ? C : ( charState2 == 'Q' ) ? Q 
	    : ( charState2 == 'E' ) ? E : ( charState2 == 'G' ) ? G : ( charState2 == 'H' ) ? H 
	    : ( charState2 == 'I' ) ? I : ( charState2 == 'L' ) ? L : ( charState2 == 'K' ) ? K
	    : ( charState2 == 'M' ) ? M : ( charState2 == 'F' ) ? F : ( charState2 == 'P' ) ? P 
	    : ( charState2 == 'S' ) ? S : ( charState2 == 'T' ) ? T : ( charState2 == 'W' ) ? W 
	    : ( charState2 == 'Y' ) ? Y : ( charState2 == 'V' ) ? V : -1;
	/*
	  index1 = -1;
	  switch (charState1) {
	  case 'A': index1 = A; break;
	  case 'R': index1 = R; break;
	  case 'N': index1 = N; break;
	  case 'D': index1 = D; break;
	  case 'C': index1 = C; break;
	  case 'Q': index1 = Q; break;
	  case 'E': index1 = E; break;
	  case 'G': index1 = G; break;
	  case 'H': index1 = H; break;
	  case 'I': index1 = I; break;
	  case 'L': index1 = L; break;
	  case 'K': index1 = K; break;
	  case 'M': index1 = M; break;
	  case 'F': index1 = F; break;
	  case 'P': index1 = P; break;
	  case 'S': index1 = S; break;
	  case 'T': index1 = T; break;
	  case 'W': index1 = W; break;
	  case 'Y': index1 = Y; break;
	  case 'V': index1 = V; break;
	  }
	  index2 = -1;
	  switch (charState2) {
	  case 'A': index2 = A; break;
	  case 'R': index2 = R; break;
	  case 'N': index2 = N; break;
	  case 'D': index2 = D; break;
	  case 'C': index2 = C; break;
	  case 'Q': index2 = Q; break;
	  case 'E': index2 = E; break;
	  case 'G': index2 = G; break;
	  case 'H': index2 = H; break;
	  case 'I': index2 = I; break;
	  case 'L': index2 = L; break;
	  case 'K': index2 = K; break;
	  case 'M': index2 = M; break;
	  case 'F': index2 = F; break;
	  case 'P': index2 = P; break;
	  case 'S': index2 = S; break;
	  case 'T': index2 = T; break;
	  case 'W': index2 = W; break;
	  case 'Y': index2 = Y; break;
	  case 'V': index2 = V; break;
	  }
	*/
	if ( (index1 != -1) && (index2 != -1) ) ctgc[index1][index2]++;
	else {
	    if ( (index1 != -1) && (index2 == -1) ) {
		switch (charState2) {
		case 'B': ctgc[index1][D] += 0.5; ctgc[index1][N] += 0.5; break;
		case 'Z': ctgc[index1][Q] += 0.5; ctgc[index1][E] += 0.5; break;
		case 'X': j = -1; while ( ++j < SIZE ) ctgc[index1][j] += 0.05; break;
		}
	    }
	    if ( (index1 == -1) && (index2 != -1) ) {
		switch (charState1) {
		case 'B': ctgc[D][index2] += 0.5; ctgc[N][index2] += 0.5; break;
		case 'Z': ctgc[Q][index2] += 0.5; ctgc[E][index2] += 0.5; break;
		case 'X': i = -1; while ( ++i < SIZE ) ctgc[i][index2] += 0.05; break;
		}
	    }
	    if ( (index1 == -1) && (index2 == -1) ) {
		switch (charState1) {
		case 'B':
		    switch (charState2) {
		    case 'B': ctgc[D][D] += 0.25; ctgc[D][N] += 0.25; ctgc[N][D] += 0.25; ctgc[N][N] += 0.25; break;
		    case 'Z': ctgc[D][Q] += 0.25; ctgc[D][E] += 0.25; ctgc[N][Q] += 0.25; ctgc[N][E] += 0.25; break;
		    case 'X': j = -1; while ( ++j < SIZE ) { ctgc[D][j] += 0.025; ctgc[N][j] += 0.025; } break;
		    }
		case 'Z':
		    switch (charState2) {
		    case 'B': ctgc[Q][D] += 0.25; ctgc[Q][N] += 0.25; ctgc[E][D] += 0.25; ctgc[E][N] += 0.25; break;
		    case 'Z': ctgc[Q][Q] += 0.25; ctgc[Q][E] += 0.25; ctgc[E][Q] += 0.25; ctgc[E][E] += 0.25; break;
		    case 'X': j = -1; while ( ++j < SIZE ) { ctgc[Q][j] += 0.025; ctgc[E][j] += 0.025; } break;
		    }
		case 'X' :
		    switch (charState2) {
		    case 'B': i = -1; while ( ++i < SIZE ) { ctgc[i][D] += 0.025; ctgc[i][N] += 0.025; } break;
		    case 'Z': i = -1; while ( ++i < SIZE ) { ctgc[i][Q] += 0.025; ctgc[i][E] += 0.025; } break;
		    case 'X': i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) ctgc[i][j] += 0.0025; } break;
		    }
		}
	    }
	}
    }
    

    public void decrement( char charState1 , char charState2 ) {
	index1 = ( charState1 == 'A' ) ? A : ( charState1 == 'R' ) ? R : ( charState1 == 'N' ) ? N 
	    : ( charState1 == 'D' ) ? D : ( charState1 == 'C' ) ? C : ( charState1 == 'Q' ) ? Q 
	    : ( charState1 == 'E' ) ? E : ( charState1 == 'G' ) ? G : ( charState1 == 'H' ) ? H 
	    : ( charState1 == 'I' ) ? I : ( charState1 == 'L' ) ? L : ( charState1 == 'K' ) ? K
	    : ( charState1 == 'M' ) ? M : ( charState1 == 'F' ) ? F : ( charState1 == 'P' ) ? P 
	    : ( charState1 == 'S' ) ? S : ( charState1 == 'T' ) ? T : ( charState1 == 'W' ) ? W 
	    : ( charState1 == 'Y' ) ? Y : ( charState1 == 'V' ) ? V : -1;
	index2 = ( charState2 == 'A' ) ? A : ( charState2 == 'R' ) ? R : ( charState2 == 'N' ) ? N 
	    : ( charState2 == 'D' ) ? D : ( charState2 == 'C' ) ? C : ( charState2 == 'Q' ) ? Q 
	    : ( charState2 == 'E' ) ? E : ( charState2 == 'G' ) ? G : ( charState2 == 'H' ) ? H 
	    : ( charState2 == 'I' ) ? I : ( charState2 == 'L' ) ? L : ( charState2 == 'K' ) ? K
	    : ( charState2 == 'M' ) ? M : ( charState2 == 'F' ) ? F : ( charState2 == 'P' ) ? P 
	    : ( charState2 == 'S' ) ? S : ( charState2 == 'T' ) ? T : ( charState2 == 'W' ) ? W 
	    : ( charState2 == 'Y' ) ? Y : ( charState2 == 'V' ) ? V : -1;
	/*
	  index1 = -1;
	  switch (charState1) {
	  case 'A': index1 = A; break;
	  case 'R': index1 = R; break;
	  case 'N': index1 = N; break;
	  case 'D': index1 = D; break;
	  case 'C': index1 = C; break;
	  case 'Q': index1 = Q; break;
	  case 'E': index1 = E; break;
	  case 'G': index1 = G; break;
	  case 'H': index1 = H; break;
	  case 'I': index1 = I; break;
	  case 'L': index1 = L; break;
	  case 'K': index1 = K; break;
	  case 'M': index1 = M; break;
	  case 'F': index1 = F; break;
	  case 'P': index1 = P; break;
	  case 'S': index1 = S; break;
	  case 'T': index1 = T; break;
	  case 'W': index1 = W; break;
	  case 'Y': index1 = Y; break;
	  case 'V': index1 = V; break;
	  }
	  index2 = -1;
	  switch (charState2) {
	  case 'A': index2 = A; break;
	  case 'R': index2 = R; break;
	  case 'N': index2 = N; break;
	  case 'D': index2 = D; break;
	  case 'C': index2 = C; break;
	  case 'Q': index2 = Q; break;
	  case 'E': index2 = E; break;
	  case 'G': index2 = G; break;
	  case 'H': index2 = H; break;
	  case 'I': index2 = I; break;
	  case 'L': index2 = L; break;
	  case 'K': index2 = K; break;
	  case 'M': index2 = M; break;
	  case 'F': index2 = F; break;
	  case 'P': index2 = P; break;
	  case 'S': index2 = S; break;
	  case 'T': index2 = T; break;
	  case 'W': index2 = W; break;
	  case 'Y': index2 = Y; break;
	  case 'V': index2 = V; break;
	  }
	*/
	if ( (index1 != -1) && (index2 != -1) ) ctgc[index1][index2]--;
	else {
	    if ( (index1 != -1) && (index2 == -1) ) {
		switch (charState2) {
		case 'B' : ctgc[index1][D] -= 0.5; ctgc[index1][N] -= 0.5; break;
		case 'Z' : ctgc[index1][Q] -= 0.5; ctgc[index1][E] -= 0.5; break;
		case 'X' : j = -1; while ( ++j < SIZE ) ctgc[index1][j] -= 0.05; break;
		}
	    }
	    if ( (index1 == -1) && (index2 != -1) ) {
		switch (charState1) {
		case 'B' : ctgc[D][index2] -= 0.5; ctgc[N][index2] -= 0.5; break;
		case 'Z' : ctgc[Q][index2] -= 0.5; ctgc[E][index2] -= 0.5; break;
		case 'X' : i = -1; while ( ++i < SIZE ) ctgc[i][index2] -= 0.05; break;
		}
	    }
	    if ( (index1 == -1) && (index2 == -1) ) {
		switch (charState1) {
		case 'B' :
		    switch (charState2) {
		    case 'B' : ctgc[D][D] -= 0.25; ctgc[D][N] -= 0.25; ctgc[N][D] -= 0.25; ctgc[N][N] -= 0.25; break;
		    case 'Z' : ctgc[D][Q] -= 0.25; ctgc[D][E] -= 0.25; ctgc[N][Q] -= 0.25; ctgc[N][E] -= 0.25; break;
		    case 'X' : j = -1; while ( ++j < SIZE ) { ctgc[D][j] -= 0.025; ctgc[N][j] -= 0.025; } break;
		    }
		case 'Z' :
		    switch (charState2) {
		    case 'B' : ctgc[Q][D] -= 0.25; ctgc[Q][N] -= 0.25; ctgc[E][D] -= 0.25; ctgc[E][N] -= 0.25; break;
		    case 'Z' : ctgc[Q][Q] -= 0.25; ctgc[Q][E] -= 0.25; ctgc[E][Q] -= 0.25; ctgc[E][E] -= 0.25; break;
		    case 'X' : j = -1; while ( ++j < SIZE ) { ctgc[Q][j] -= 0.025; ctgc[E][j] -= 0.025; } break;
		    }
		case 'X' :
		    switch (charState2) {
		    case 'B' : i = -1; while ( ++i < SIZE ) { ctgc[i][D] -= 0.025; ctgc[i][N] -= 0.025; } break;
		    case 'Z' : i = -1; while ( ++i < SIZE ) { ctgc[i][Q] -= 0.025; ctgc[i][E] -= 0.025; } break;
		    case 'X' : i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) ctgc[i][j] -= 0.0025; } break;
		    }
		}
	    }
	}
    }

    public double get( char charState1 , char charState2 ) {
	index1 = ( charState1 == 'A' ) ? A : ( charState1 == 'R' ) ? R : ( charState1 == 'N' ) ? N 
	    : ( charState1 == 'D' ) ? D : ( charState1 == 'C' ) ? C : ( charState1 == 'Q' ) ? Q 
	    : ( charState1 == 'E' ) ? E : ( charState1 == 'G' ) ? G : ( charState1 == 'H' ) ? H 
	    : ( charState1 == 'I' ) ? I : ( charState1 == 'L' ) ? L : ( charState1 == 'K' ) ? K
	    : ( charState1 == 'M' ) ? M : ( charState1 == 'F' ) ? F : ( charState1 == 'P' ) ? P 
	    : ( charState1 == 'S' ) ? S : ( charState1 == 'T' ) ? T : ( charState1 == 'W' ) ? W 
	    : ( charState1 == 'Y' ) ? Y : ( charState1 == 'V' ) ? V : -1;
	index2 = ( charState2 == 'A' ) ? A : ( charState2 == 'R' ) ? R : ( charState2 == 'N' ) ? N 
	    : ( charState2 == 'D' ) ? D : ( charState2 == 'C' ) ? C : ( charState2 == 'Q' ) ? Q 
	    : ( charState2 == 'E' ) ? E : ( charState2 == 'G' ) ? G : ( charState2 == 'H' ) ? H 
	    : ( charState2 == 'I' ) ? I : ( charState2 == 'L' ) ? L : ( charState2 == 'K' ) ? K
	    : ( charState2 == 'M' ) ? M : ( charState2 == 'F' ) ? F : ( charState2 == 'P' ) ? P 
	    : ( charState2 == 'S' ) ? S : ( charState2 == 'T' ) ? T : ( charState2 == 'W' ) ? W 
	    : ( charState2 == 'Y' ) ? Y : ( charState2 == 'V' ) ? V : -1;
	/*
	  index1 = -1;
	  switch (charState1) {
	  case 'A': index1 = A; break;
	  case 'R': index1 = R; break;
	  case 'N': index1 = N; break;
	  case 'D': index1 = D; break;
	  case 'C': index1 = C; break;
	  case 'Q': index1 = Q; break;
	  case 'E': index1 = E; break;
	  case 'G': index1 = G; break;
	  case 'H': index1 = H; break;
	  case 'I': index1 = I; break;
	  case 'L': index1 = L; break;
	  case 'K': index1 = K; break;
	  case 'M': index1 = M; break;
	  case 'F': index1 = F; break;
	  case 'P': index1 = P; break;
	  case 'S': index1 = S; break;
	  case 'T': index1 = T; break;
	  case 'W': index1 = W; break;
	  case 'Y': index1 = Y; break;
	  case 'V': index1 = V; break;
	  }
	  index2 = -1;
	  switch (charState2) {
	  case 'A': index2 = A; break;
	  case 'R': index2 = R; break;
	  case 'N': index2 = N; break;
	  case 'D': index2 = D; break;
	  case 'C': index2 = C; break;
	  case 'Q': index2 = Q; break;
	  case 'E': index2 = E; break;
	  case 'G': index2 = G; break;
	  case 'H': index2 = H; break;
	  case 'I': index2 = I; break;
	  case 'L': index2 = L; break;
	  case 'K': index2 = K; break;
	  case 'M': index2 = M; break;
	  case 'F': index2 = F; break;
	  case 'P': index2 = P; break;
	  case 'S': index2 = S; break;
	  case 'T': index2 = T; break;
	  case 'W': index2 = W; break;
	  case 'Y': index2 = Y; break;
	  case 'V': index2 = V; break;
	  }
	*/
	if ( (index1 != -1) && (index2 != -1) ) return ctgc[index1][index2];
	else {
	    double val = 0;
	    if ( (index1 != -1) && (index2 == -1) ) {
		switch (charState2) {
		case 'B' : return 0.5 * (ctgc[index1][D] + ctgc[index1][N]);
		case 'Z' : return 0.5 * (ctgc[index1][Q] + ctgc[index1][E]);
		case 'X' : j = -1; while ( ++j < SIZE ) val += ctgc[index1][j]; return 0.05 * val;
		}
	    }
	    if ( (index1 == -1) && (index2 != -1) ) {
		switch (charState1) {
		case 'B' : return 0.5 * (ctgc[D][index2] + ctgc[N][index2]);
		case 'Z' : return 0.5 * (ctgc[Q][index2] + ctgc[E][index2]);
		case 'X' : i = -1; while ( ++i < SIZE ) val += ctgc[i][index2]; return 0.05 * val;
		}
	    }
	    if ( (index1 == -1) && (index2 == -1) ) {
		switch (charState1) {
		case 'B' :
		    switch (charState2) {
		    case 'B' : return 0.25 * (ctgc[D][D] + ctgc[D][N] + ctgc[N][D] + ctgc[N][N]);
		    case 'Z' : return 0.25 * (ctgc[D][Q] + ctgc[D][E] + ctgc[N][Q] + ctgc[N][E]);
		    case 'X' : j = -1; while ( ++j < SIZE ) val += (ctgc[D][j]+ctgc[N][j]); return 0.025 * val;
		    }
		case 'Z' :
		    switch (charState2) {
		    case 'B' : return 0.25 * (ctgc[Q][D] + ctgc[Q][N] + ctgc[E][D] + ctgc[E][N]);
		    case 'Z' : return 0.25 * (ctgc[Q][Q] + ctgc[Q][E] + ctgc[E][Q] + ctgc[E][E]); 
		    case 'X' : j = -1; while ( ++j < SIZE ) val += (ctgc[Q][j]+ctgc[E][j]); return 0.025 * val;
		    }
		case 'X' :
		    switch (charState2) {
		    case 'B' : i = -1; while ( ++i < SIZE ) val += (ctgc[i][D]+ctgc[i][N]); return 0.025 * val;
		    case 'Z' : i = -1; while ( ++i < SIZE ) val += (ctgc[i][Q]+ctgc[i][E]); return 0.025 * val; 
		    case 'X' : i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) val += ctgc[i][j]; } return 0.0025 * val;
		    }
		}
	    }
	}
	//System.out.println(charState1 + " " + charState2);
    	return -1.0;
    }

    public double get(int row , int col) {
	if ( (row >= 0) && (row < SIZE) && (col >= 0) && (col < SIZE) ) return ctgc[row][col];
	return -1.0;
    }

    public void set(int row , int col , double x) {
	/*if ( (row >= 0) && (row < SIZE) && (col >= 0) && (col < SIZE) && (x >= 0) )*/ ctgc[row][col] = x;
    }

    public int size() { 
	return SIZE; 
    }

    public AminoAcidMeter copy() {
	aam_copy = new AminoAcidMeter();
	i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) aam_copy.set(i , j , ctgc[i][j]); }
	return aam_copy;
    }

    public double getBowkerSymmetryTestPvalue() {
	chi2 = 0; dof = 0;
	i = -1;
	while ( ++i < SIZE ) {
	    j = -1;
	    while ( ++j < i ) {
		if ( ctgc[i][j] + ctgc[j][i] > 0 ) { 
		    chi2 += (ctgc[i][j] - ctgc[j][i]) * (ctgc[i][j] - ctgc[j][i]) / (ctgc[i][j] + ctgc[j][i]); dof++;
		}
	    }
	}
	return StatisticalTest.getKhi2Pvalue(chi2 , dof);
    }

    public double getEvansHoenigSymmetryTestPvalue() {
	chi2 = 0; dof = 0;
	p = -1;
	while ( ++p < SIZE - 1 ) {
	    up = 0; down = 0;
	    i = -1; while ( ++i < SIZE - p ) { up += ctgc[p+i][i] - ctgc[i][p+i]; down += ctgc[p+i][i] + ctgc[i][p+i]; }
	    if ( down > 0 ) { chi2 += up * up / down; dof++; }
	}
	return StatisticalTest.getKhi2Pvalue(chi2 , dof);
    }

    public double getStuartMarginalSymmetryTestValue() {
	dof = SIZE-1;
	matrix = new Matrix(dof , dof); sumRow = new Matrix(dof , 1); sumCol = new Matrix(1 , dof); 
	i = -1;
	while ( ++i < dof ) {
	    sum_r = 0; sum_c = 0;
	    j = -1; 
	    while ( ++j < SIZE ) { 
		sum_r += (int) ctgc[i][j]; sum_c += (int) ctgc[j][i]; 
		if ( j < dof ) matrix.set(i , j , (int) ctgc[i][j]); 
	    }
	    sumRow.set(i , 0 , sum_r); sumCol.set(0 , i , sum_c);
	}
	
	while ( true ) {
	    s = -1; while ( ++s < dof ) if ( sumRow.get(s , 0) + sumCol.get(0 , s) - 2 * matrix.get(s , s) == 0 ) break;
	    if ( s == dof ) break;
	    else {  // row and col s are full of non-diagonal 0 => one must remove this dof
		dof--; index = new int[dof];
		i = -1; j = -1; while ( ++i < dof+1 ) { if ( i != s ) index[++j] = i; }
		matrix = matrix.getMatrix(index , index); 
		sumRow = sumRow.getMatrix(index , 0 , 0);
		sumCol = sumCol.getMatrix(0 , 0 , index);
	    }
	}

	vector_U = new Matrix(dof , 1); matrix_V = new Matrix(dof , dof);
	i = -1;
	while ( ++i < dof ) {
	    vector_U.set(i , 0 , sumRow.get(i , 0) - sumCol.get(0 , i)); 
	    if ( vector_U.get(i , 0) == 0 ) vector_U.set(i , 0 , 0);
	    j = -1;
	    while ( ++j < dof ) {
		if ( i != j ) matrix_V.set(i , j , - matrix.get(i , j) - matrix.get(j , i));
		else matrix_V.set(i , i , sumRow.get(i , 0) + sumCol.get(0 , i) - 2 * matrix.get(i , i));
		if ( matrix_V.get(i , j) == 0 ) matrix_V.set(i , j , 0);
	    }
	}

	//while ( Math.abs(matrix_V.det()) < 1e-9 ) { 	// verifying the non-singularity of matrix_V
	while ( ! (new LUDecomposition(matrix_V)).isNonsingular() ) {
	    dof--; matrix_V = matrix_V.getMatrix(0 , dof-1 , 0 , dof-1); vector_U = vector_U.getMatrix(0 , dof-1 , 0 , 0);
	}

	while ( true ) {
	    vector_UT = vector_U.transpose();
	    matrix = matrix_V.inverse();
	    matrix = vector_UT.times(matrix);
	    matrix = matrix.times(vector_U);
	    chi2 = matrix.get(0 , 0);

	    if ( chi2 >= 0 ) break;
	    else {
		/*System.out.println(""); i = -1; while ( ++i < dof ) { j = -1; while ( ++j < dof ) 
		  System.out.print( (matrix_V.get(i , j) + "     ").substring(0 , 6) ); System.out.println(""); }
		  System.out.println(" " + matrix_V.det());*/
		do {
		    dof--; matrix_V = matrix_V.getMatrix(0 , dof-1 , 0 , dof-1); vector_U = vector_U.getMatrix(0 , dof-1 , 0 , 0);
		} while ( ! (new LUDecomposition(matrix_V)).isNonsingular() );
	    }
	}
	vector_U = null; vector_UT = null; matrix_V = null; matrix = null; sumRow = null; sumCol = null;
	return chi2;
    }
	
    public double getStuartMarginalSymmetryTestPvalue() {
	return StatisticalTest.getKhi2Pvalue(getStuartMarginalSymmetryTestValue() , SIZE-1);
    }
	




}
