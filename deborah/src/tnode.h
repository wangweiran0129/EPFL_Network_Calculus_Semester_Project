/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   *
 *   luca.bisti@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef DEBORAHTNODE_H
#define DEBORAHTNODE_H

// enable the following to dump the PI and PI_EQ curves at each tnode during computation of LUDB
//#define TNODE_VERBOSE_LUDB

// enable this to print debug messages in the experimental LUDB computation function
//#define TNODE_DEBUG_LUDB_EXPERIMENTAL

// enable this to dump the constraints set for each simplex performed
//#define TNODE_DEBUG_LUDB_CONSTRAINTS


#define TNODE_LUDB_EPSILON	0.0000001
#define TNODE_LUDB_INFINITY	9999999999.0

#define LUDB_MODE_EXACT		-100
#define LUDB_MODE_HEURISTIC 0x7FFFFFFF


#include "types.h"
#include "flow.h"
#include "node.h"
#include "curve.h"
#include "linearsegment.h"
#include "pseudoaffine.h"
#include <vector>
#include "cpoly.h"

namespace deborah {

	typedef std::vector<uint> IndexSet;
	typedef std::vector<CPoly> ConstraintSet;
	typedef std::vector<IndexSet> IndexSetVector;


/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class TNode{
public:
    TNode();
	TNode(uint num_flows);
	TNode(uint flow_id, uint num_tnodes, uint num_leaves);
	bool has_leaf();
	uint num_children();
    ~TNode();

	// computes a specific instance of the PI
	PseudoAffine& getPI(const IndexSet &idx_set, ConstraintSet &c_set);
	// computes the equivalent service curve for the given PI
	PseudoAffine& computeEquivalent(PseudoAffine &p, const IndexSet &set, ConstraintSet &c_set);
	// compute a specific instance of the equivalent service curve at this node
	PseudoAffine &getE(const IndexSet &idx_set);
	// compute the Delay Bound at this node
	Poly& getDelayBound(const IndexSet &idx_set, ConstraintSet &c_set);

	// compute LUDB (old-fashioned algorithm, serial enumeration)
	double computeLUDB(Poly &delay, Poly &solution, int nvars, IndexSetVector *good_combos);
	// compute LUDB (new implementation, recursive enumeration)
	uint64 computeLUDB_experimental(Poly &delay, Poly &solution, int nvars, IndexSetVector &good_combos, int max_good_combos, bool eta);

	void Print(bool subtree);
	void PrintIndexSet(const IndexSet &set, bool compact);

	// index and pointer to a Flow object, representing flow (h,k)
	uint flow;
	Flow *pflow;

	// array of TNode children, representing the set S(h,k)
	TNode **tnodes;
	uint num_tnodes;

	// array of node indexes, representing set C(h,k)
	uint *leaves;
	uint num_leaves;

	// The PI(h,k) pseudo-affine service curve
	PseudoAffine *pi;
	PseudoAffine *pi_eq;
	// The PI_C(h,k) rate-latency service curve, obtained as (X)Bj
	LinearSegment pi_c;

	// Structured set of sets of indexes for the current node:
	IndexSet index_set;
	// Each element is the number of indexes belonging to the corresponding node
	// the first node is the current one, then children follow
	std::vector<uint> index_schema;
	// Each element is the TNode number associated to the corresponding index
	IndexSet index_nodes;

	// Methods devoted to handle the recursive set data structure
	bool getCombo(IndexSet &set, uint64 n);
	int incCombo(IndexSet &set, int index);
	uint64 getNumCombos();
	// Return the Ith subset of the node's index set
	bool getSubset(IndexSet &set, int n_child);
	bool getSubset(const IndexSet &set, int n_child, IndexSet &subset);

};

}

#endif
