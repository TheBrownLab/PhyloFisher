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

public class StatisticalTest {

    // Wilcoxon
    private static int size;
    private static ArrayList<WilcoxonTuple> wt_array;
    private static ArrayList<Double> rank, distribution, array_temp;
    private static double n, wp, wm, w, pValue;
    private static double r, norm;
    private static WilcoxonTuple wt;

    // standard normal distribution
    private static double snd_epsilon, snd_z, snd_p, snd_f, snd_ff;

    // sign test
    private static int st_plus, st_minus, st_n, st_k;
    private static double st_f, st_ff, st_p;

    // lnGamma
    private final static double[] LANCZOS = {76.18009172947146 , -86.50532032941677 , 24.01409824083091 , 
					     -1.231739572450155 , 0.001208650973866179 , -0.000005395239384953 };
    private static double lg_z, lg_f, lg_ff;

    // incomplete gamma
    private static double ig_k, ig_f, ig_s, ig_b, ig_c, ig_d, ig_an;



    // khi2 distribution
    private static double k2d_p, k2d_cte1, k2d_cte2, k2d_cte3, k2d_f, k2d_ff, k2d_s, k2d_z, k2d_epsilon;



    private static int i, j, b;


    private StatisticalTest() {}


    public static double getWilcoxonMatchedPairsSignedRanksPvalue(ArrayList<Double> sample1 , ArrayList<Double> sample2) {
	size = Math.min(sample1.size() , sample2.size());
	if ( size == 0 ) return -1;
	//### creating the sorted list of WilcoxonTuple ###
	wt_array = new ArrayList<WilcoxonTuple>(0); n = 0;
	i = -1;
	while ( ++i < size ) {
	    wt = new WilcoxonTuple(sample1.get(i).doubleValue() , sample2.get(i).doubleValue());
	    if ( wt.sign != 0 ) { wt_array.add( wt ); n = 1; break; }
	}
	while ( ++i < size ) {
	    wt = new WilcoxonTuple(sample1.get(i).doubleValue() , sample2.get(i).doubleValue());
	    if ( wt.sign != 0 ) {
		if ( wt.absoluteDifference >= wt_array.get(((int) (n - 1))).absoluteDifference ) wt_array.add( wt ); // wt is the largest
		else {
		    j = -1;  // inserting wt inside wt_array
		    while ( ++j < n ) {
			if ( wt.absoluteDifference < wt_array.get(j).absoluteDifference ) {
			    wt_array.add(j , wt); 
			    break;
			}
		    }
		}
		n++; // increasing the size of wt_array
	    }
	}
	//### computing the ranks ###
	rank = new ArrayList<Double>(((int) n));
	i = 0; while ( ++i <= n ) rank.add(new Double(i));
	i = -1;
	while ( ++i < n ) {
	    wt = wt_array.get(i);
	    j = i;
	    while ( ++j < n ) {
		if ( wt_array.get(j).absoluteDifference != wt.absoluteDifference ) break;
	    }
	    if ( --j > i ) {
		r = 0; size = 0;
		b = i - 1; while ( ++b <= j ) { r += rank.get(b).doubleValue(); size++; }
		r /= (double) size;
		b = i - 1; while ( ++b <= j ) rank.set(b , new Double(r));
		i = j;
	    }
	}
	//### computing Wilcoxon parameters and p-value ###
	wp = 0; wm = 0;
	i = -1;
	while ( ++i < n ) {
	    if ( wt_array.get(i).sign < 0 ) wm += rank.get(i).doubleValue();
	    else wp += rank.get(i).doubleValue();
	    /*System.out.println( wt_array.get(i).sign + " " + wt_array.get(i).absoluteDifference + " " + rank.get(i).doubleValue() );*/
	}
	//System.out.println(n + "  " + wp + " " + wm);
	pValue = 0;
	if ( n <= 50 ) {        // computing the exact distribution
	    distribution = new ArrayList<Double>(0);
	    distribution.add(new Double(1.0));
	    distribution.add(new Double(1.0));
	    norm = 2.0;
	    b = 1;
	    while ( ++b <= n ) {
		size = b*(b+1)/2 + 1;
		norm *= 2.0;
		array_temp = new ArrayList<Double>(0);
		i = -1; 
		while ( ++i < distribution.size() ) array_temp.add( distribution.get(i) );
		while ( ++i <= size ) { array_temp.add( new Double(0) ); distribution.add( new Double(0) ); }
		i = -1;
		while ( ++i < size ) 
		    distribution.set(i , new Double(array_temp.get(i).doubleValue() 
						    + array_temp.get(size - i - 1).doubleValue()));
	    }
	    array_temp = null;
	    w = Math.min(wm , wp);
	    i = -1; while ( ++i <= w ) pValue += distribution.get(i).doubleValue();
	    pValue = 2.0 * pValue / norm;
	}
	else {                  // approximating with standard normal distribution
	    w = Math.max(wm , wp); 
	    pValue = getStandardNormalDistributionPvalue( (w - 0.5 - n*(n+1)/4) / Math.sqrt( n*(n+1)*(2*n+1)/24 ) ); 
	}
	return pValue;
    }

    private static class WilcoxonTuple {
	public double absoluteDifference, sign;
	public WilcoxonTuple(double x1 , double x2) {
	    this.absoluteDifference = Math.abs(x1 - x2); 
	    if ( this.absoluteDifference == 0 ) this.sign = 0;
	    else this.sign = (x1 - x2) / this.absoluteDifference;
	}
    }

    public static double getSignTestPvalue(ArrayList<Double> sample1 , ArrayList<Double> sample2) {
	size = Math.min(sample1.size() , sample2.size());
	if ( size == 0 ) return -1;
	st_plus = 0; st_minus = 0; 
	i = -1; 
	while ( ++i < size ) {
	    if ( sample1.get(i).doubleValue() < sample2.get(i).doubleValue() ) st_minus++;
	    else {
		if ( sample1.get(i).doubleValue() > sample2.get(i).doubleValue() ) st_plus++;
	    }
	}
	//System.out.println(" " + st_plus + " " + st_minus + " ");
	if ( (st_plus == 0) && (st_minus == 0) ) return 1;
	//### binomial distribution ###
	st_n = st_plus + st_minus; st_k = Math.min(st_plus , st_minus); st_f = 1; st_p = st_f;
	i = 1; while ( i <= st_k ) { st_f = st_f * (st_n--) / (i++); st_p += st_f; }
	st_p /= Math.pow(2 , st_plus + st_minus - 1 );
	if ( ! (Double.isNaN(st_p) || Double.isInfinite(st_p)) ) return st_p;
	//### approximation with normal distribution ###
	return getStandardNormalDistributionPvalue( ((double) (Math.abs(st_minus - st_plus) - 1)) / Math.sqrt(st_plus + st_minus) );
    }


    public static double getStandardNormalDistributionPvalue(double z) {
	snd_epsilon = 0.00001; snd_z = z; snd_p = 0; snd_f = Math.exp(-0.5 * snd_z * snd_z);
	while ( snd_z < 5 ) {
	    snd_z += snd_epsilon; snd_ff = Math.exp(-0.5 * snd_z * snd_z);
	    if ( Double.isNaN(snd_ff) || Double.isInfinite(snd_ff) ) break;
	    snd_p +=  snd_f + snd_ff; snd_f = snd_ff;
	}
	snd_p *= snd_epsilon / Math.sqrt(2 * Math.PI);
	return snd_p;
    }


    public static double getLnGammaFunction(double z) {
	if ( z == 0 ) return -1;
	if ( z == 1 ) return 0;
	lg_z = z;
	lg_ff = (--lg_z) + 5.5;
	lg_ff -= (lg_z + 0.5) * Math.log(lg_ff);
	i = -1;
	lg_f = 1; while ( ++i < 6 ) lg_f += LANCZOS[i] / (++lg_z);
	//return Math.log( Math.sqrt(2 * Math.PI) * lg_f ) - lg_ff;
	return Math.log( 2.5066282746310005024157652848110 * lg_f ) - lg_ff;
    }

    public static double getIncompleteGammaFunction(double k , double z) {
	if ( (k <= 0) || (z < 0) ) return -1;
	if ( z < k + 1.0 ) {
	    ig_k = k; ig_f = 1.0 / ig_k; ig_s = ig_f;
	    while ( true ) {
		ig_f *= z / (++ig_k); 
		if ( Double.isNaN(ig_f) || Double.isInfinite(ig_f) || (ig_f < 0.000000001) || (ig_f < 0) ) break;
		ig_s += ig_f;
	    }
	    return ig_s * Math.exp(k * Math.log(z) - z - getLnGammaFunction(k));
	}
	else {
	    ig_b = z + 1 - k; ig_c = Double.MAX_VALUE; ig_d = 1.0 / ig_b; ig_s = ig_d;
	    i = 0;
	    while ( true ) {
		ig_an = ((double) ++i) * (k - ((double)i));
		ig_b++; ig_b++;
		ig_d = ig_an * ig_d + ig_b; if ( Math.abs(ig_d) < Float.MIN_VALUE ) ig_d = Double.MIN_VALUE;
		ig_c = ig_an / ig_c + ig_b; if ( Math.abs(ig_c) < Float.MIN_VALUE ) ig_c = Float.MIN_VALUE;
		ig_d = 1.0 / ig_d;
		ig_f = ig_d * ig_c;
		if ( Double.isNaN(ig_f) || Double.isInfinite(ig_f) || (ig_f == 1.0) ) break;
		ig_s *= ig_f;
	    }
	    return 1.0 - ig_s * Math.exp(k * Math.log(z) - z - getLnGammaFunction(k));
	}
    }
    
    public static double getKhi2Pvalue(double z , int dof) {
	if ( z * dof == 0 ) return 1;
	k2d_p = 1.0 - getIncompleteGammaFunction(((double) dof)/2.0 , z/2.0);
	if ( k2d_p > 1.0 ) System.out.println(z + "(" + dof + ")=" + k2d_p);
	return k2d_p;
    }














    
    public static double getKhi2Pvalue_(double z , int dof) {
	//### the solution is approximatively 1.0 ###
	if ( (z == 0) || ((z < 1.0) && (dof >= 10)) ) return 0.9999999999;
	//### exact computing of the khi2 p-value ###
	if ( dof <= 35 ) {
	    //### analytical formula when dof/2 is an integer ###
	    if ( dof % 2 == 0 ) {
		k2d_cte1 = z / 2; k2d_cte2 = Math.exp(-k2d_cte1);
		k2d_p = k2d_cte2; b = dof / 2; k2d_f = 1; 
		i = 0;
		while ( ++i < b ) {
		    k2d_f *= k2d_cte1;
		    k2d_ff = k2d_f * k2d_cte2 / getGammaFunctionValue(i+1);
		    if ( Double.isNaN(k2d_ff) || Double.isInfinite(k2d_ff) ) break;
		    k2d_p += k2d_ff;
		}
		return k2d_p;
	    }
	    /*k2d_cte1 = z / 2; k2d_cte2 = ((double) dof) / 2.0; k2d_cte3 = getGammaFunctionValue(k2d_cte2);
	    k2d_f = Math.pow(k2d_cte1 , k2d_cte2); k2d_ff = 1; k2d_p = k2d_f / (k2d_cte2 * k2d_cte3); i = 0;
	    k2d_f = k2d_cte2 * Math.log(k2d_cte1); 
	    while ( true ) {
		//k2d_f *= k2d_cte1; k2d_ff *= (++i);
		//k2d_s = k2d_f / ((++k2d_cte2) * k2d_ff * k2d_cte3);
		k2d_f += StrictMath.log(k2d_cte1) - StrictMath.log(++i); 
		k2d_s = StrictMath.exp(k2d_f) / ((++k2d_cte2) * k2d_cte3);
		//System.out.println("  " + k2d_f + " " + k2d_ff + " " + k2d_s);
		if ( Double.isNaN(k2d_s) || Double.isInfinite(k2d_s) || (k2d_s < 0.0000001) ) { System.out.print(k2d_s + "  "); break;}
		if ( i % 2 == 0 ) k2d_p += k2d_s;
		else k2d_p -= k2d_s;
		if ( Math.abs(k2d_p) <= 1.0 ) return 1.0 - k2d_p;
		System.out.print("#" + k2d_f + " " + k2d_p + "#");
		}*/

	    k2d_cte1 = z / 2; k2d_cte2 = ((double) dof) / 2.0;
	    i = 51; k2d_p = 0;
	    while ( --i > 0 ) {
		if ( i % 2 == 0 ) k2d_p = (i/2) * k2d_cte1 / (k2d_cte2 + i + k2d_p);
		else k2d_p = - (k2d_cte2 + (i-1)/2) * k2d_cte1 / (k2d_cte2 + i + k2d_p);
	    }
	    k2d_p = 1.0 - Math.pow(k2d_cte1 , k2d_cte2) * Math.exp(-k2d_cte1) / ( (k2d_cte2 + k2d_p) * getGammaFunctionValue(k2d_cte2) );
	    if ( (0.0001 < k2d_p) && (k2d_p < 1.0) ) return k2d_p;

	    k2d_cte2 = ((double) dof) / 2.0; k2d_cte1 = Math.pow(2 , 1.0 + k2d_cte2) * getGammaFunctionValue(k2d_cte2); k2d_cte2 -= 1.0;
	    k2d_p = 0; k2d_epsilon = 0.0005; k2d_z = z;
	    k2d_f = Math.pow(k2d_z , k2d_cte2) / Math.exp(k2d_z / 2.0);
	    while ( true ) {
		k2d_z += k2d_epsilon;
		k2d_ff = Math.pow(k2d_z , k2d_cte2) / Math.exp(k2d_z / 2.0);
		k2d_s = (3 * Math.max(k2d_f , k2d_ff) - Math.min(k2d_f , k2d_ff)) / k2d_cte1;
		if ( Double.isNaN(k2d_s) || Double.isInfinite(k2d_s) || (k2d_s < 0.000001) ) break;
		k2d_p += k2d_epsilon * k2d_s;
		k2d_f = k2d_ff;
	    }
	    return k2d_p;
	}
	return 0.5 * getStandardNormalDistributionPvalue( ( Math.pow(z/((double)dof),1.0/3.0) - 1.0 + 2.0/(9.0*((double)dof)) ) 
							  / Math.sqrt(2.0/(9.0*((double)dof))) );
    }		  
		    
		    


    
    private static double getGammaFunctionValue(double k) {  // only for k such that 1 <= 2*k <= 40 and 2*k is an integer
	switch ((int)(2*k)) {                                // if k is an integer, then this returns (k-1)!
	case 1: return Math.sqrt(Math.PI);
	case 2: return 1.0;
	case 3: return Math.sqrt(Math.PI) / 2.0;
	case 4: return 1.0;
	case 5: return (3.0 / 4.0) * Math.sqrt(Math.PI);
	case 6: return 2.0;
	case 7: return (15.0 / 8.0) * Math.sqrt(Math.PI);
	case 8: return 6.0;
	case 9: return (105.0 / 16.0) * Math.sqrt(Math.PI);
	case 10: return 24.0;
	case 11: return (945.0 / 32.0) * Math.sqrt(Math.PI);
	case 12: return 120.0;
	case 13: return (10395.0 / 64.0) * Math.sqrt(Math.PI);
	case 14: return 720.0;
	case 15: return (135135.0 / 128.0) * Math.sqrt(Math.PI);
	case 16: return 5040.0;
	case 17: return (2027025.0 / 256.0) * Math.sqrt(Math.PI);
	case 18: return 40320.0;
	case 19: return (34459425.0 / 512.0) * Math.sqrt(Math.PI);
	case 20: return 362880.0;
	case 21: return (654729075.0 / 1024.0) * Math.sqrt(Math.PI);
	case 22: return 3628800.0;
	case 23: return (13749310575.0 / 2048.0) * Math.sqrt(Math.PI);
	case 24: return 39916800.0;
	case 25: return (316234143225.0 / 4096.0) * Math.sqrt(Math.PI);
	case 26: return 479001600.0;
	case 27: return (7905853580625.0 / 8192.0) * Math.sqrt(Math.PI);
	case 28: return 6227020800.0;
	case 29: return (213458046676875.0 / 16384.0) * Math.sqrt(Math.PI);
	case 30: return 87178291200.0;
	case 31: return (6190283353629374.0 / 32768.0) * Math.sqrt(Math.PI);
	case 32: return 1307674368000.0;
	case 33: return (191898783962510624.0 / 65536.0) * Math.sqrt(Math.PI);
	case 34: return 20922789888000.0;
	case 35: return (6332659870762850304.0 / 131072.0) * Math.sqrt(Math.PI);
	case 36: return 355687428096000.0;
	case 37: return (221643095476699758592.0 / 262144.0) * Math.sqrt(Math.PI);
	case 38: return 6402373705728000.0;
	case 39: return (8200794532637890838528.0 / 524288.0) * Math.sqrt(Math.PI);
	case 40: return 121645100408832000.0;
	}
	return -1.0;
    }


}

