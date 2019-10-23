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

public class NucleotideMeter extends Meter {

    private final char[] character = {'A','G','C','T','U','M','R','W','S','Y','K','B','D','H','V','N','X','?','-'};
    //  A Adenosine             G Guanine               C Cytosine        T Thymine / U  Uracil
    public static final int A = 0; public static final int G = 1; 
    public static final int C = 2; public static final int T = 3;
    //                  M A,C                            R A,G                            W A,T
    private final int[] M = {A , C}; private final int[] R = {A , G}; private final int[] W = {A , T}; 
    //                  S C,G                            Y C,T                            K G,T
    private final int[] S = {C , G}; private final int[] Y = {C , T}; private final int[] K = {G , T}; 
    //                  B C,G,T                              D A,G,T
    private final int[] B = {C , G , T}; private final int[] D = {A , G , T}; 
    //                  H A,C,T                              V A,C,G 
    private final int[] H = {A , C , T}; private final int[] V = {A , C , G}; 
    //                  N,X,? unknown       -  gap     
    private final int SIZE = 4; 
    private final double THD = 1.0 / 3.0; private final double STH = 1.0 / 6.0; 
    private final double NTH = 1.0 / 9.0; private final double TTH = 1.0 / 12.0; 

    private double[][] ctgc; 
    private int i, j, s, m, index1, index2, dof, p;
    private double chi2, up, down, sum_r, sum_c, val;
    private char c1, c2;
    private NucleotideMeter num_copy;
    private Matrix sumRow, sumCol, vector_U, vector_UT, matrix_V, matrix;
    private int[] index;

    public NucleotideMeter() {
	ctgc = new double[SIZE][SIZE];
	i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) ctgc[i][j] = 0; }
    }

    public NucleotideMeter(String nucleotideSequence1 , String nucleotideSequence2) {
	this();
	m = Math.min(nucleotideSequence1.length() , nucleotideSequence2.length());
	s = -1;	while ( ++s < m ) increment( nucleotideSequence1.charAt(s) , nucleotideSequence2.charAt(s) );
    }   
		
    public void increment( char charState1 , char charState2 ) {
	if ( (charState1 != '-') && (charState2 != '-') && (charState1 != '?') && (charState2 != '?') ) {
	    index1 = -1;
	    switch (charState1) {
	    case 'A': index1 = A; break;
	    case 'C': index1 = C; break;
	    case 'G': index1 = G; break;
	    case 'T': index1 = T; break;
	    }
	    index2 = -1;
	    switch (charState2) {
	    case 'A': index2 = A; break;
	    case 'C': index2 = C; break;
	    case 'G': index2 = G; break;
	    case 'T': index2 = T; break;
	    }
	    if ( (index1 != -1) && (index2 != -1) ) ctgc[index1][index2]++;
	    else {
		if ( (index1 != -1) && (index2 == -1) ) {
		    switch (charState2) {
		    case 'M': j = -1; while ( ++j < M.length ) ctgc[index1][M[j]] += 0.5; break;
		    case 'R': j = -1; while ( ++j < R.length ) ctgc[index1][R[j]] += 0.5; break;
		    case 'W': j = -1; while ( ++j < W.length ) ctgc[index1][W[j]] += 0.5; break;
		    case 'S': j = -1; while ( ++j < S.length ) ctgc[index1][S[j]] += 0.5; break;
		    case 'Y': j = -1; while ( ++j < Y.length ) ctgc[index1][Y[j]] += 0.5; break;
		    case 'K': j = -1; while ( ++j < K.length ) ctgc[index1][K[j]] += 0.5; break;
		    case 'B': j = -1; while ( ++j < B.length ) ctgc[index1][B[j]] += THD; break;
		    case 'D': j = -1; while ( ++j < D.length ) ctgc[index1][D[j]] += THD; break;
		    case 'H': j = -1; while ( ++j < H.length ) ctgc[index1][H[j]] += THD; break;
		    case 'V': j = -1; while ( ++j < V.length ) ctgc[index1][V[j]] += THD; break;
		    case 'X': j = -1; while ( ++j < SIZE ) ctgc[index1][j] += 0.25; break;
		    }
		}
		if ( (index1 == -1) && (index2 != -1) ) {
		    j = -1;
		    switch (charState1) {
		    case 'M': while ( ++j < M.length ) ctgc[M[j]][index2] += 0.5; break;
		    case 'R': while ( ++j < R.length ) ctgc[R[j]][index2] += 0.5; break;
		    case 'W': while ( ++j < W.length ) ctgc[W[j]][index2] += 0.5; break;
		    case 'S': while ( ++j < S.length ) ctgc[S[j]][index2] += 0.5; break;
		    case 'Y': while ( ++j < Y.length ) ctgc[Y[j]][index2] += 0.5; break;
		    case 'K': while ( ++j < K.length ) ctgc[K[j]][index2] += 0.5; break;
		    case 'B': while ( ++j < B.length ) ctgc[B[j]][index2] += THD; break;
		    case 'D': while ( ++j < D.length ) ctgc[D[j]][index2] += THD; break;
		    case 'H': while ( ++j < H.length ) ctgc[H[j]][index2] += THD; break;
		    case 'V': while ( ++j < V.length ) ctgc[V[j]][index2] += THD; break;
		    case 'X': while ( ++j < SIZE ) ctgc[j][index2] += 0.25; break;
		    }
		}
		if ( (index1 == -1) && (index2 == -1) ) {
		    i = -1; 
		    switch (charState1) {
		    case 'M':
			while( ++i < M.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[M[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[M[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[M[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[M[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[M[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[M[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[M[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[M[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[M[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[M[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[M[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'R':
			while( ++i < R.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[R[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[R[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[R[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[R[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[R[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[R[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[R[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[R[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[R[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[R[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[R[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'W':
			while( ++i < W.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[W[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[W[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[W[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[W[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[W[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[W[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[W[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[W[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[W[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[W[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[W[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'S':
			while( ++i < S.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[S[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[S[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[S[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[S[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[S[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[S[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[S[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[S[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[S[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[S[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[S[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'Y':
			while( ++i < Y.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[Y[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[Y[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[Y[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[Y[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[Y[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[Y[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[Y[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[Y[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[Y[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[Y[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[Y[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'K':
			while( ++i < K.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[K[i]][M[j]] += 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[K[i]][R[j]] += 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[K[i]][W[j]] += 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[K[i]][S[j]] += 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[K[i]][Y[j]] += 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[K[i]][K[j]] += 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[K[i]][B[j]] += STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[K[i]][D[j]] += STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[K[i]][H[j]] += STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[K[i]][V[j]] += STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[K[i]][j] += 0.125; break;
			    }
			}
			break;
		    case 'B':
			while( ++i < B.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[B[i]][M[j]] += STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[B[i]][R[j]] += STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[B[i]][W[j]] += STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[B[i]][S[j]] += STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[B[i]][Y[j]] += STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[B[i]][K[j]] += STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[B[i]][B[j]] += NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[B[i]][D[j]] += NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[B[i]][H[j]] += NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[B[i]][V[j]] += NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[B[i]][j] += TTH; break;
			    }
			}
			break;
		    case 'D':
			while( ++i < D.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[D[i]][M[j]] += STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[D[i]][R[j]] += STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[D[i]][W[j]] += STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[D[i]][S[j]] += STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[D[i]][Y[j]] += STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[D[i]][K[j]] += STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[D[i]][B[j]] += NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[D[i]][D[j]] += NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[D[i]][H[j]] += NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[D[i]][V[j]] += NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[D[i]][j] += TTH; break;
			    }
			}
			break;
		    case 'H':
			while( ++i < H.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[H[i]][M[j]] += STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[H[i]][R[j]] += STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[H[i]][W[j]] += STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[H[i]][S[j]] += STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[H[i]][Y[j]] += STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[H[i]][K[j]] += STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[H[i]][B[j]] += NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[H[i]][D[j]] += NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[H[i]][H[j]] += NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[H[i]][V[j]] += NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[H[i]][j] += TTH; break;
			    }
			}
			break;
		    case 'V':
			while( ++i < V.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[V[i]][M[j]] += STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[V[i]][R[j]] += STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[V[i]][W[j]] += STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[V[i]][S[j]] += STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[V[i]][Y[j]] += STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[V[i]][K[j]] += STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[V[i]][B[j]] += NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[V[i]][D[j]] += NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[V[i]][H[j]] += NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[V[i]][V[j]] += NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[V[i]][j] += TTH; break;
			    }
			}
			break;
		    case 'X' :
			while( ++i < SIZE ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[i][M[j]] += 0.125; break;
			    case 'R': while ( ++j < R.length ) ctgc[i][R[j]] += 0.125; break;
			    case 'W': while ( ++j < W.length ) ctgc[i][W[j]] += 0.125; break;
			    case 'S': while ( ++j < S.length ) ctgc[i][S[j]] += 0.125; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[i][Y[j]] += 0.125; break;
			    case 'K': while ( ++j < K.length ) ctgc[i][K[j]] += 0.125; break;
			    case 'B': while ( ++j < K.length ) ctgc[i][B[j]] += TTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[i][D[j]] += TTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[i][H[j]] += TTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[i][V[j]] += TTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[i][j] += 0.0625; break;
			    }
			}
			break;
		    }
		}
	    }
	}
    }

    public void decrement( char charState1 , char charState2 ) {
	if ( (charState1 != '-') && (charState2 != '-') && (charState1 != '?') && (charState2 != '?') ) {
	    index1 = -1;
	    switch (charState1) {
	    case 'A': index1 = A; break;
	    case 'C': index1 = C; break;
	    case 'G': index1 = G; break;
	    case 'T': index1 = T; break;
	    }
	    index2 = -1;
	    switch (charState2) {
	    case 'A': index2 = A; break;
	    case 'C': index2 = C; break;
	    case 'G': index2 = G; break;
	    case 'T': index2 = T; break;
	    }
	    if ( (index1 != -1) && (index2 != -1) ) ctgc[index1][index2]--;
	    else {
		if ( (index1 != -1) && (index2 == -1) ) {
		    switch (charState2) {
		    case 'M': j = -1; while ( ++j < M.length ) ctgc[index1][M[j]] -= 0.5; break;
		    case 'R': j = -1; while ( ++j < R.length ) ctgc[index1][R[j]] -= 0.5; break;
		    case 'W': j = -1; while ( ++j < W.length ) ctgc[index1][W[j]] -= 0.5; break;
		    case 'S': j = -1; while ( ++j < S.length ) ctgc[index1][S[j]] -= 0.5; break;
		    case 'Y': j = -1; while ( ++j < Y.length ) ctgc[index1][Y[j]] -= 0.5; break;
		    case 'K': j = -1; while ( ++j < K.length ) ctgc[index1][K[j]] -= 0.5; break;
		    case 'B': j = -1; while ( ++j < B.length ) ctgc[index1][B[j]] -= THD; break;
		    case 'D': j = -1; while ( ++j < D.length ) ctgc[index1][D[j]] -= THD; break;
		    case 'H': j = -1; while ( ++j < H.length ) ctgc[index1][H[j]] -= THD; break;
		    case 'V': j = -1; while ( ++j < V.length ) ctgc[index1][V[j]] -= THD; break;
		    case 'X': j = -1; while ( ++j < SIZE ) ctgc[index1][j] -= 0.25; break;
		    }
		}
		if ( (index1 == -1) && (index2 != -1) ) {
		    j = -1;
		    switch (charState1) {
		    case 'M': while ( ++j < M.length ) ctgc[M[j]][index2] -= 0.5; break;
		    case 'R': while ( ++j < R.length ) ctgc[R[j]][index2] -= 0.5; break;
		    case 'W': while ( ++j < W.length ) ctgc[W[j]][index2] -= 0.5; break;
		    case 'S': while ( ++j < S.length ) ctgc[S[j]][index2] -= 0.5; break;
		    case 'Y': while ( ++j < Y.length ) ctgc[Y[j]][index2] -= 0.5; break;
		    case 'K': while ( ++j < K.length ) ctgc[K[j]][index2] -= 0.5; break;
		    case 'B': while ( ++j < B.length ) ctgc[B[j]][index2] -= THD; break;
		    case 'D': while ( ++j < D.length ) ctgc[D[j]][index2] -= THD; break;
		    case 'H': while ( ++j < H.length ) ctgc[H[j]][index2] -= THD; break;
		    case 'V': while ( ++j < V.length ) ctgc[V[j]][index2] -= THD; break;
		    case 'X': while ( ++j < SIZE ) ctgc[j][index2] -= 0.25; break;
		    }
		}
		if ( (index1 == -1) && (index2 == -1) ) {
		    i = -1; 
		    switch (charState1) {
		    case 'M':
			while( ++i < M.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[M[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[M[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[M[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[M[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[M[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[M[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[M[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[M[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[M[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[M[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[M[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'R':
			while( ++i < R.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[R[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[R[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[R[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[R[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[R[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[R[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[R[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[R[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[R[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[R[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[R[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'W':
			while( ++i < W.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[W[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[W[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[W[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[W[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[W[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[W[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[W[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[W[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[W[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[W[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[W[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'S':
			while( ++i < S.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[S[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[S[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[S[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[S[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[S[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[S[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[S[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[S[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[S[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[S[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[S[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'Y':
			while( ++i < Y.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[Y[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[Y[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[Y[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[Y[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[Y[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[Y[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[Y[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[Y[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[Y[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[Y[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[Y[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'K':
			while( ++i < K.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[K[i]][M[j]] -= 0.25; break;
			    case 'R': while ( ++j < R.length ) ctgc[K[i]][R[j]] -= 0.25; break;
			    case 'W': while ( ++j < W.length ) ctgc[K[i]][W[j]] -= 0.25; break;
			    case 'S': while ( ++j < S.length ) ctgc[K[i]][S[j]] -= 0.25; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[K[i]][Y[j]] -= 0.25; break;
			    case 'K': while ( ++j < K.length ) ctgc[K[i]][K[j]] -= 0.25; break;
			    case 'B': while ( ++j < K.length ) ctgc[K[i]][B[j]] -= STH; break;
			    case 'D': while ( ++j < K.length ) ctgc[K[i]][D[j]] -= STH; break;
			    case 'H': while ( ++j < K.length ) ctgc[K[i]][H[j]] -= STH; break;
			    case 'V': while ( ++j < K.length ) ctgc[K[i]][V[j]] -= STH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[K[i]][j] -= 0.125; break;
			    }
			}
			break;
		    case 'B':
			while( ++i < B.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[B[i]][M[j]] -= STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[B[i]][R[j]] -= STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[B[i]][W[j]] -= STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[B[i]][S[j]] -= STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[B[i]][Y[j]] -= STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[B[i]][K[j]] -= STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[B[i]][B[j]] -= NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[B[i]][D[j]] -= NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[B[i]][H[j]] -= NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[B[i]][V[j]] -= NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[B[i]][j] -= TTH; break;
			    }
			}
			break;
		    case 'D':
			while( ++i < D.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[D[i]][M[j]] -= STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[D[i]][R[j]] -= STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[D[i]][W[j]] -= STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[D[i]][S[j]] -= STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[D[i]][Y[j]] -= STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[D[i]][K[j]] -= STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[D[i]][B[j]] -= NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[D[i]][D[j]] -= NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[D[i]][H[j]] -= NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[D[i]][V[j]] -= NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[D[i]][j] -= TTH; break;
			    }
			}
			break;
		    case 'H':
			while( ++i < H.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[H[i]][M[j]] -= STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[H[i]][R[j]] -= STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[H[i]][W[j]] -= STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[H[i]][S[j]] -= STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[H[i]][Y[j]] -= STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[H[i]][K[j]] -= STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[H[i]][B[j]] -= NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[H[i]][D[j]] -= NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[H[i]][H[j]] -= NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[H[i]][V[j]] -= NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[H[i]][j] -= TTH; break;
			    }
			}
			break;
		    case 'V':
			while( ++i < V.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[V[i]][M[j]] -= STH; break;
			    case 'R': while ( ++j < R.length ) ctgc[V[i]][R[j]] -= STH; break;
			    case 'W': while ( ++j < W.length ) ctgc[V[i]][W[j]] -= STH; break;
			    case 'S': while ( ++j < S.length ) ctgc[V[i]][S[j]] -= STH; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[V[i]][Y[j]] -= STH; break;
			    case 'K': while ( ++j < K.length ) ctgc[V[i]][K[j]] -= STH; break;
			    case 'B': while ( ++j < K.length ) ctgc[V[i]][B[j]] -= NTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[V[i]][D[j]] -= NTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[V[i]][H[j]] -= NTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[V[i]][V[j]] -= NTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[V[i]][j] -= TTH; break;
			    }
			}
			break;
		    case 'X' :
			while( ++i < SIZE ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) ctgc[i][M[j]] -= 0.125; break;
			    case 'R': while ( ++j < R.length ) ctgc[i][R[j]] -= 0.125; break;
			    case 'W': while ( ++j < W.length ) ctgc[i][W[j]] -= 0.125; break;
			    case 'S': while ( ++j < S.length ) ctgc[i][S[j]] -= 0.125; break;
			    case 'Y': while ( ++j < Y.length ) ctgc[i][Y[j]] -= 0.125; break;
			    case 'K': while ( ++j < K.length ) ctgc[i][K[j]] -= 0.125; break;
			    case 'B': while ( ++j < K.length ) ctgc[i][B[j]] -= TTH; break;
			    case 'D': while ( ++j < K.length ) ctgc[i][D[j]] -= TTH; break;
			    case 'H': while ( ++j < K.length ) ctgc[i][H[j]] -= TTH; break;
			    case 'V': while ( ++j < K.length ) ctgc[i][V[j]] -= TTH; break;
			    case 'X': while ( ++j < SIZE ) ctgc[i][j] -= 0.0625; break;
			    }
			}
			break;
		    }
		}
	    }
	}
    }

    public void increment( int row , int col ) {
	switch (row) {
	case A: 
	    switch (col) {
	    case A: increment('A','A'); break;
	    case C: increment('A','C'); break;
	    case G: increment('A','G'); break;
	    case T: increment('A','T'); break;
	    }
	case C: 
	    switch (col) {
	    case A: increment('C','A'); break;
	    case C: increment('C','C'); break;
	    case G: increment('C','G'); break;
	    case T: increment('C','T'); break;
	    }
	case G: 
	    switch (col) {
	    case A: increment('G','A'); break;
	    case C: increment('G','C'); break;
	    case G: increment('G','G'); break;
	    case T: increment('G','T'); break;
	    }
	case T: 
	    switch (col) {
	    case A: increment('T','A'); break;
	    case C: increment('T','C'); break;
	    case G: increment('T','G'); break;
	    case T: increment('T','T'); break;
	    }
	}
    }

    public void decrement( int row , int col ) {
	switch (row) {
	case A: 
	    switch (col) {
	    case A: decrement('A','A'); break;
	    case C: decrement('A','C'); break;
	    case G: decrement('A','G'); break;
	    case T: decrement('A','T'); break;
	    }
	case C: 
	    switch (col) {
	    case A: decrement('C','A'); break;
	    case C: decrement('C','C'); break;
	    case G: decrement('C','G'); break;
	    case T: decrement('C','T'); break;
	    }
	case G: 
	    switch (col) {
	    case A: decrement('G','A'); break;
	    case C: decrement('G','C'); break;
	    case G: decrement('G','G'); break;
	    case T: decrement('G','T'); break;
	    }
	case T: 
	    switch (col) {
	    case A: decrement('T','A'); break;
	    case C: decrement('T','C'); break;
	    case G: decrement('T','G'); break;
	    case T: decrement('T','T'); break;
	    }
	}
    }

    public double get( char charState1 , char charState2 ) {
	if ( (charState1 != '-') && (charState2 != '-') && (charState1 != '?') && (charState2 != '?') ) {
	    index1 = -1;
	    switch (charState1) {
	    case 'A': index1 = A; break;
	    case 'C': index1 = C; break;
	    case 'G': index1 = G; break;
	    case 'T': index1 = T; break;
	    }
	    index2 = -1;
	    switch (charState2) {
	    case 'A': index2 = A; break;
	    case 'C': index2 = C; break;
	    case 'G': index2 = G; break;
	    case 'T': index2 = T; break;
	    }
	    if ( (index1 != -1) && (index2 != -1) ) return ctgc[index1][index2];
	    else {
		val = 0;
		if ( (index1 != -1) && (index2 == -1) ) {
		    j = -1;
		    switch (charState2) {
		    case 'M': while ( ++j < M.length ) val += 0.5 * ctgc[index1][M[j]]; return val; 
		    case 'R': while ( ++j < R.length ) val += 0.5 * ctgc[index1][R[j]]; return val; 
		    case 'W': while ( ++j < W.length ) val += 0.5 * ctgc[index1][W[j]]; return val;
		    case 'S': while ( ++j < S.length ) val += 0.5 * ctgc[index1][S[j]]; return val;
		    case 'Y': while ( ++j < Y.length ) val += 0.5 * ctgc[index1][Y[j]]; return val;
		    case 'K': while ( ++j < K.length ) val += 0.5 * ctgc[index1][K[j]]; return val;
		    case 'B': while ( ++j < B.length ) val += THD * ctgc[index1][B[j]]; return val;
		    case 'D': while ( ++j < D.length ) val += THD * ctgc[index1][D[j]]; return val;
		    case 'H': while ( ++j < H.length ) val += THD * ctgc[index1][H[j]]; return val;
		    case 'V': while ( ++j < V.length ) val += THD * ctgc[index1][V[j]]; return val;
		    case 'X': while ( ++j < SIZE ) val += 0.25 * ctgc[index1][j]; return val;
		    }
		}
		if ( (index1 == -1) && (index2 != -1) ) {
		    j = -1;
		    switch (charState1) {
		    case 'M': while ( ++j < M.length ) val += 0.5 * ctgc[M[j]][index2]; return val; 
		    case 'R': while ( ++j < R.length ) val += 0.5 * ctgc[R[j]][index2]; return val; 
		    case 'W': while ( ++j < W.length ) val += 0.5 * ctgc[W[j]][index2]; return val; 
		    case 'S': while ( ++j < S.length ) val += 0.5 * ctgc[S[j]][index2]; return val; 
		    case 'Y': while ( ++j < Y.length ) val += 0.5 * ctgc[Y[j]][index2]; return val; 
		    case 'K': while ( ++j < K.length ) val += 0.5 * ctgc[K[j]][index2]; return val; 
		    case 'B': while ( ++j < B.length ) val += THD * ctgc[B[j]][index2]; return val; 
		    case 'D': while ( ++j < D.length ) val += THD * ctgc[D[j]][index2]; return val; 
		    case 'H': while ( ++j < H.length ) val += THD * ctgc[H[j]][index2]; return val; 
		    case 'V': while ( ++j < V.length ) val += THD * ctgc[V[j]][index2]; return val; 
		    case 'X': while ( ++j < SIZE ) val += 0.25 * ctgc[j][index2]; return val; 
		    }
		}
		if ( (index1 == -1) && (index2 == -1) ) {
		    i = -1; 
		    switch (charState1) {
		    case 'M':
			while( ++i < M.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[M[i]][M[j]]; return val; 
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[M[i]][R[j]]; return val; 
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[M[i]][W[j]]; return val; 
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[M[i]][S[j]]; return val; 
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[M[i]][Y[j]]; return val; 
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[M[i]][K[j]]; return val; 
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[M[i]][B[j]]; return val; 
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[M[i]][D[j]]; return val; 
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[M[i]][H[j]]; return val; 
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[M[i]][V[j]]; return val; 
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[M[i]][j]; return val; 
			    }
			}
			break;
		    case 'R':
			while( ++i < R.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[R[i]][M[j]]; return val; 
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[R[i]][R[j]]; return val; 
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[R[i]][W[j]]; return val; 
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[R[i]][S[j]]; return val; 
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[R[i]][Y[j]]; return val; 
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[R[i]][K[j]]; return val; 
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[R[i]][B[j]]; return val; 
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[R[i]][D[j]]; return val; 
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[R[i]][H[j]]; return val; 
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[R[i]][V[j]]; return val; 
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[R[i]][j]; return val; 
			    }
			}
			break;
		    case 'W':
			while( ++i < W.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[W[i]][M[j]]; return val; 
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[W[i]][R[j]]; return val; 
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[W[i]][W[j]]; return val; 
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[W[i]][S[j]]; return val; 
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[W[i]][Y[j]]; return val; 
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[W[i]][K[j]]; return val; 
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[W[i]][B[j]]; return val; 
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[W[i]][D[j]]; return val; 
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[W[i]][H[j]]; return val; 
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[W[i]][V[j]]; return val; 
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[W[i]][j]; return val; 
			    }
			}
			break;
		    case 'S':
			while( ++i < S.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[S[i]][M[j]]; return val; 
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[S[i]][R[j]]; return val; 
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[S[i]][W[j]]; return val; 
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[S[i]][S[j]]; return val; 
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[S[i]][Y[j]]; return val; 
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[S[i]][K[j]]; return val; 
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[S[i]][B[j]]; return val; 
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[S[i]][D[j]]; return val; 
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[S[i]][H[j]]; return val; 
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[S[i]][V[j]]; return val; 
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[S[i]][j]; return val; 
			    }
			}
			break;
		    case 'Y':
			while( ++i < Y.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[Y[i]][M[j]]; return val; 
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[Y[i]][R[j]]; return val; 
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[Y[i]][W[j]]; return val; 
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[Y[i]][S[j]]; return val; 
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[Y[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[Y[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[Y[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[Y[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[Y[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[Y[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[Y[i]][j]; return val;
			    }
			}
			break;
		    case 'K':
			while( ++i < K.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.25 * ctgc[K[i]][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += 0.25 * ctgc[K[i]][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += 0.25 * ctgc[K[i]][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += 0.25 * ctgc[K[i]][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += 0.25 * ctgc[K[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += 0.25 * ctgc[K[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += STH * ctgc[K[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += STH * ctgc[K[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += STH * ctgc[K[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += STH * ctgc[K[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += 0.125 * ctgc[K[i]][j]; return val;
			    }
			}
			break;
		    case 'B':
			while( ++i < B.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += STH * ctgc[B[i]][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += STH * ctgc[B[i]][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += STH * ctgc[B[i]][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += STH * ctgc[B[i]][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += STH * ctgc[B[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += STH * ctgc[B[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += NTH * ctgc[B[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += NTH * ctgc[B[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += NTH * ctgc[B[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += NTH * ctgc[B[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += TTH * ctgc[B[i]][j]; return val;
			    }
			}
			break;
		    case 'D':
			while( ++i < D.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += STH * ctgc[D[i]][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += STH * ctgc[D[i]][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += STH * ctgc[D[i]][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += STH * ctgc[D[i]][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += STH * ctgc[D[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += STH * ctgc[D[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += NTH * ctgc[D[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += NTH * ctgc[D[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += NTH * ctgc[D[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += NTH * ctgc[D[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += TTH * ctgc[D[i]][j]; return val;
			    }
			}
			break;
		    case 'H':
			while( ++i < H.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += STH * ctgc[H[i]][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += STH * ctgc[H[i]][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += STH * ctgc[H[i]][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += STH * ctgc[H[i]][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += STH * ctgc[H[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += STH * ctgc[H[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += NTH * ctgc[H[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += NTH * ctgc[H[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += NTH * ctgc[H[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += NTH * ctgc[H[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += TTH * ctgc[H[i]][j]; return val;
			    }
			}
			break;
		    case 'V':
			while( ++i < V.length ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += STH * ctgc[V[i]][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += STH * ctgc[V[i]][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += STH * ctgc[V[i]][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += STH * ctgc[V[i]][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += STH * ctgc[V[i]][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += STH * ctgc[V[i]][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += NTH * ctgc[V[i]][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += NTH * ctgc[V[i]][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += NTH * ctgc[V[i]][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += NTH * ctgc[V[i]][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += TTH * ctgc[V[i]][j]; return val;
			    }
			}
			break;
		    case 'X' :
			while( ++i < SIZE ) {
			    j = -1; 
			    switch (charState2) {
			    case 'M': while ( ++j < M.length ) val += 0.125 * ctgc[i][M[j]]; return val;
			    case 'R': while ( ++j < R.length ) val += 0.125 * ctgc[i][R[j]]; return val;
			    case 'W': while ( ++j < W.length ) val += 0.125 * ctgc[i][W[j]]; return val;
			    case 'S': while ( ++j < S.length ) val += 0.125 * ctgc[i][S[j]]; return val;
			    case 'Y': while ( ++j < Y.length ) val += 0.125 * ctgc[i][Y[j]]; return val;
			    case 'K': while ( ++j < K.length ) val += 0.125 * ctgc[i][K[j]]; return val;
			    case 'B': while ( ++j < K.length ) val += TTH * ctgc[i][B[j]]; return val;
			    case 'D': while ( ++j < K.length ) val += TTH * ctgc[i][D[j]]; return val;
			    case 'H': while ( ++j < K.length ) val += TTH * ctgc[i][H[j]]; return val;
			    case 'V': while ( ++j < K.length ) val += TTH * ctgc[i][V[j]]; return val;
			    case 'X': while ( ++j < SIZE ) val += 0.0625 * ctgc[i][j]; return val;
			    }
			}
			break;
		    }
		}
	    }
	}
	return 0;
    }

    public double get(int row , int col) {
	if ( (row >= 0) && (row < SIZE) && (col >= 0) && (col < SIZE) ) return ctgc[row][col];
	return -1.0;
    }

    public void set(int row , int col , double x) {
	if ( (row >= 0) && (row < SIZE) && (col >= 0) && (col < SIZE) && (x >= 0) ) ctgc[row][col] = x;
    }

    public int size() { 
	return SIZE; 
    }

    /*public char getChar(int charIndex) {
      switch (charIndex) {
      case A: return 'A';
      case C: return 'C';
      case G: return 'G';
      case T: return 'T';
      default: return '-';
      }
      }*/

    public NucleotideMeter copy() {
	num_copy = new NucleotideMeter();
	i = -1; while ( ++i < SIZE ) { j = -1; while ( ++j < SIZE ) num_copy.set(i , j , ctgc[i][j]); }
	return num_copy;
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
