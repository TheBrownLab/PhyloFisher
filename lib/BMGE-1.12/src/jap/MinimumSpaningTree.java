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
import java.util.*;

/**
 *  The <code>MinimumSpaningTree</code> class implements a minimum spaning tree from a symetrical matrix representation of a valued graph. No branch must be represented by an Integer.MAX_VALUE edge length
 **/

public class MinimumSpaningTree {

    private TreeSet<DissimilarityEntry> dissimilarityTuple, minimumSpaningTree;
    private Comparator<Object> comparator;
    private DissimilarityEntry de;
    private ArrayList<Integer> partition, partitionName, bipartition1, bipartition2;
    
    private int i, j, p, k, pi, pj;
    
    public MinimumSpaningTree(Matrix mat) {
	comparator = new Comparator<Object>() {
	    public int compare (Object o1 , Object o2) {
		if ( (! (o1 instanceof DissimilarityEntry)) || (! (o2 instanceof DissimilarityEntry)) ) throw new ClassCastException();
		if ( ((DissimilarityEntry) o1).distance < ((DissimilarityEntry) o2).distance ) return -1;
		if ( ((DissimilarityEntry) o1).distance > ((DissimilarityEntry) o2).distance ) return 1;
		if ( ((DissimilarityEntry) o1).label1 < ((DissimilarityEntry) o2).label1 ) return -1;
		if ( ((DissimilarityEntry) o1).label1 > ((DissimilarityEntry) o2).label1 ) return 1;
		if ( ((DissimilarityEntry) o1).label2 < ((DissimilarityEntry) o2).label2 ) return -1;
		if ( ((DissimilarityEntry) o1).label2 > ((DissimilarityEntry) o2).label2 ) return 1;
		return 0;
	    }
	    public boolean equals (Object o1 , Object o2) {
		if ( (! (o1 instanceof DissimilarityEntry)) || (! (o2 instanceof DissimilarityEntry)) ) throw new ClassCastException();
		return (((DissimilarityEntry) o1).distance == ((DissimilarityEntry) o2).distance)
		&& (((DissimilarityEntry) o1).label1 == ((DissimilarityEntry) o2).label1)
		&& (((DissimilarityEntry) o1).label2 == ((DissimilarityEntry) o2).label2);
	    }
	};
	dissimilarityTuple = new TreeSet<DissimilarityEntry>(comparator);
	k = Math.min(mat.getRowDimension() , mat.getColumnDimension());
	i = -1;
	while ( ++i < k ) {
	    j = -1;
	    while ( ++j < k ) 
		dissimilarityTuple.add(new DissimilarityEntry(i , j , (mat.get(i , j) + mat.get(j , i))/2.0));
	}
	
	minimumSpaningTree = new TreeSet<DissimilarityEntry>(comparator);
	bipartition1 = new ArrayList<Integer>(0);
	bipartition2 = new ArrayList<Integer>(0);
	partition = new ArrayList<Integer>(0);
	partitionName = new ArrayList<Integer>(0);
	p = -1;
	while ( minimumSpaningTree.size() < k-1 ) { // a MST has k-1 edges, where k is the number of nodes
	    de = dissimilarityTuple.first();
	    if ( (! partition.contains(de.label1)) && (! partition.contains(de.label2)) ) {
		partition.add(new Integer(de.label1)); partitionName.add(new Integer(++p));
		partition.add(new Integer(de.label2)); partitionName.add(new Integer(p));
		minimumSpaningTree.add(de);
	    }
	    else {
		if ( partition.contains(de.label1) && (! partition.contains(de.label2)) ) {
		    partition.add(new Integer(de.label2)); partitionName.add( partitionName.get(partition.indexOf(de.label1)) );
		    minimumSpaningTree.add(de);
		}
		else {
		    if ( (! partition.contains(de.label1)) && partition.contains(de.label2) ) {
			partition.add(new Integer(de.label1)); partitionName.add( partitionName.get(partition.indexOf(de.label2)) );
			minimumSpaningTree.add(de);
		    }
		    else {
			pi = partitionName.get(partition.indexOf(de.label1)).intValue();
			pj = partitionName.get(partition.indexOf(de.label2)).intValue();
			if ( pi != pj ) {
			    if ( minimumSpaningTree.size() == k-2 ) { // 'de' is the last edge of the MST
				i = -1;
				while ( ++i < partitionName.size() ) {
				    if ( partitionName.get(i).intValue() == pi ) bipartition1.add(partition.get(i));
				    else bipartition2.add(partition.get(i));
				}
			    }
			    else {
				i = -1;
				while ( ++i < partitionName.size() ) 
				    if ( (partitionName.get(i).intValue() == pi) || (partitionName.get(i).intValue() == pj) ) 
					partitionName.set(i , new Integer(Math.min(pi , pj)));
			    }
			    minimumSpaningTree.add(de);
			}
		    }
		}
	    }
	    dissimilarityTuple.remove(de);
	}
	partition = null; partitionName = null; dissimilarityTuple = null;
    }
    
    private class DissimilarityEntry {
	public int label1, label2; public double distance;
	public DissimilarityEntry(int label1, int label2, double distance) {
	    this.label1 = label1; this.label2 = label2; this.distance = distance;
	}
    }
    
    public ArrayList<Integer> getFirstBipartition() {
	return bipartition1;
    }
    
    public ArrayList<Integer> getSecondBipartition() {
	return bipartition2;
    }
    
    public TreeSet<DissimilarityEntry> getMinimumSpaningTree() {
	return minimumSpaningTree;
    }
}
