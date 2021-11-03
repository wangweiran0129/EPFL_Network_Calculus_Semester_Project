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
#include "tnode.h"
#include <iostream>
#include "simplex.h"
#include <math.h>
#include <time.h>
#include "rng.h"
#include <stdio.h>

// Random Number Generator used in LUDB routines
static RNG rng_ludb(9981);


namespace deborah {

TNode::TNode()
{
	num_tnodes = 0;
	num_leaves = 0;
	tnodes = NULL;
	leaves = NULL;
	flow = 0;
	pflow = NULL;
	pi = NULL;
	pi_eq = NULL;
}


TNode::TNode(uint flow_id, uint nt, uint nl)
{
	flow = flow_id;
	num_tnodes = nt;
	num_leaves = nl;
	if(num_tnodes>0)
		tnodes = new TNode* [num_tnodes];
	else tnodes = NULL;
	if(num_leaves>0)
		leaves = new uint [num_leaves];
	else leaves = NULL;
	pflow = NULL;
	pi = NULL;
	index_set.clear();
	index_schema.clear();
	index_nodes.clear();
}


TNode::~TNode()
{
	//Warning: don't use pflow, since the pointed Flow object could have been already destroyed
	//delete child tnodes
	if(tnodes != NULL)
	{
		//delete children
		if(num_tnodes>0)
			for(int i=0;i<num_tnodes;i++)
				delete tnodes[i];
		//delete array of pointers to children
		delete [] tnodes;
		tnodes = NULL;
	}
	//delete leaves array
	if(leaves != NULL)
	{
		delete [] leaves;
		leaves = NULL;
	}

	if(pi != NULL)
	{
		delete pi;
		pi = NULL;
	}

	if(pi_eq != NULL)
	{
		delete pi_eq;
		pi_eq = NULL;
	}

	flow = 0;
	pflow = NULL;
}


void TNode::Print(bool subtree)
{
	if(pflow == NULL)
	{
		printf("TNode: flow id #%02i has NULL pointer to Flow, internal error!\n",flow);
		return;
	}
	printf("TNode: flow (%i,%i)   flow_id=#%02i\n", pflow->src+1, pflow->exit+1, flow);

	//Print leaves
	printf("       leaf nodes:");
	if(num_leaves == 0) printf(" none");
	else {
		for(int i=0; i<num_leaves; i++)
			printf(" %02i",leaves[i]+1);
	}
	printf("\n");

	//Print directly nested flows: set S(h,k)
	printf("       child flows:");
	if(num_tnodes == 0) printf(" none");
	else {
		for(int i=0; i<num_tnodes; i++)
			printf(" (%i,%i)", tnodes[i]->pflow->src+1, tnodes[i]->pflow->exit+1);
	}
	printf("\n");

	//Print PI_C service curve associated to this node
	printf("       PI_C: ");
	if(num_leaves == 0) printf("none\n");
	else printf("latency=%.2f rate=%.2f\n",pi_c.x,pi_c.rate);

	//printf("       PI: "); pi->Print();

	//Print set of indexes available at this node
	printf("       Indexes: [");
	uint idx_pos = 0;
	uint idx_combos=1;
	for(int i=0; i<index_schema.size(); i++)
	{
		uint setsize = index_schema[i];
		printf(" [");
		for(int j=idx_pos; j<idx_pos+setsize; j++)
		{
			printf(" %i",index_set[j]);
			idx_combos *= index_set[j];
		}
		printf(" ]");
		idx_pos += setsize;
	}
	printf(" ] --> combinations = %i\n",idx_combos);

	//Recurse printing the subtree, if requested
	if(subtree)
	{
		for(int i=0; i<num_tnodes; i++)
			tnodes[i]->Print(true);
	}
}


void TNode::PrintIndexSet(const IndexSet &set, bool compact)
{
	printf("[");
	uint idx_pos = 0;
	if(compact)
	{
		for(int i=0; i<index_schema.size(); i++)
		{
			uint setsize = index_schema[i];
			printf(" [");
			for(int j=idx_pos; j<idx_pos+setsize; j++)
			{
				printf(" %i",set[j]);
			}
			printf(" ]");
			idx_pos += setsize;
		}
	}
	else
	{
		for(int i=0; i<index_schema.size(); i++)
		{
			uint setsize = index_schema[i];
			printf(" [");
			for(int j=idx_pos; j<idx_pos+setsize; j++)
			{
				printf(" %i/%i",set[j],index_set[j]);
			}
			printf(" ]");
			idx_pos += setsize;
		}
	}
	printf(" ]");
}

bool TNode::has_leaf()
{
	if(num_leaves > 0) return true;
	return false;
}


uint TNode::num_children()
{
	return num_tnodes;
}


/*	This method computes the PI() pseudo-affine service curve
	for the current t-node.
	1. Retrieve the PI curves of all the child flows
	2. Compute their equivalent service curves (--> at the child node!)
	3. Convolve them all, convolve also with the local PI and return
*/
PseudoAffine& TNode::getPI(const IndexSet &idx_set, ConstraintSet &c_set)
{
#ifdef TNODE_VERBOSE_LUDB
	printf(">>>>>> TNode #%02i: getPI() --> IndexSet = ",flow); PrintIndexSet(idx_set,false); printf("\n");
#endif
	//if this is a leaf flow, return the node's PI directly

	if(num_tnodes == 0)
	{
#ifdef TNODE_VERBOSE_LUDB
		printf("<<<<<< TNode #%02i: getPI() finished (no children)\n",flow);
#endif
		return *pi;
	}

	// initialize the PI with the local pi_c curve
	if(has_leaf())
		pi->MakeLatencyRate(pi_c);
	else pi->Zero();

	// now convolve it with the equivalent service curve of each child flow
	for(int i=0; i<num_tnodes; i++)
	{
		// get index set for this child
		IndexSet child_set;
		if(getSubset(idx_set, i, child_set) == false) continue;
		// get PI curve of the child flow
		PseudoAffine &child_pi = tnodes[i]->getPI(child_set, c_set);
		// compute equivalent service curve at the child node
		PseudoAffine &child_eq = tnodes[i]->computeEquivalent(child_pi, child_set, c_set);
		// convolve the equivalent service curve into the local PI curve
		pi->Convolve(child_eq);
	}

#ifdef TNODE_VERBOSE_LUDB
	printf("<<<<<< TNode #%02i: getPI() finished\n",flow);
#endif
	return *pi;
}



/*	This method computes an equivalent service queue according to corollary 2.7
	Input parameters:
		- the pseudo-affine curve that will be combined with the
	  	flow corresponding to the current node in order to determine the equivalent
		service curve experienced by the flow
		- the index set for the current node, specifying the index value selected for the max() expression
		- the constraint set, which will be filled in accordingly
*/
PseudoAffine& TNode::computeEquivalent(PseudoAffine &p, const IndexSet &set, ConstraintSet &c_set)
{
	if(pi_eq == NULL)
	{
		printf("TNode::computeEquivalent(tnode flow #%02i): pi_eq == NULL!!!\n",flow);
		return p;
	}

#ifdef TNODE_VERBOSE_LUDB
	printf(">>> TNode #%02i: computeEquivalent(",flow); PrintIndexSet(set, false); printf(")\n");
#endif

	pi_eq->Zero();

	//determine which is the ith element of the max() expression that we have been asked for
	uint max_index = set[0]; // max_index is the index of the chosen 'ith' element
	uint max_index_limit = index_set[0];
	//range-check (for bugs)
	if(max_index >= max_index_limit)
	{
		printf("!!! BUG !!! Requested Index=%i but highest available is %i\n",max_index,max_index_limit);
		return *pi_eq;
	}

	/*****  now compute the value of the chosen element of the max() expression  *****/

	uint nvars = p.delay.getNumVariables();
	Poly max_term(nvars);  // max_term = Zero
	CPoly diseq(nvars,flow); // will store the constraints in the form "ax + b >= 0"
	Poly ptemp(nvars);

	// an index value exceeding the number of stages of the PI corresponds to the case max()={0}
	if(max_index < p.num_stages)
	{
		// this is the case where max(0, X1,..Xn) == Xi
		max_term.SetB(pflow->burst); //max_term = sigma (of the flow represented by the current t-node)
		// compute the term [ (sigma - sigma_i) / rho_i ]+  (where i == max_index)
		max_term = (max_term - p.stages2[max_index].sigma) / p.stages2[max_index].rho;
#ifdef TNODE_VERBOSE_LUDB
		printf("MAX_TERM = "); max_term.Print(); printf("\n");
#endif
		//add the constraint (sigma - sigma_i) > 0 to implement the + apice, i.e. max(0, max_term)
		diseq = max_term;
		// store the constraint only if it's not a tautology
		if(!diseq.isTautology()) c_set.push_back(diseq);

		// now add the constraints (sigma - sigma_i)/rho_i >= all the others stages of the local PI
		for(int i=0; i<p.num_stages; i++)
		{
			if(i == max_index) continue;
			ptemp.Zero();
			ptemp.SetB(pflow->burst);
			ptemp = (ptemp - p.stages2[i].sigma) / p.stages2[i].rho;
			// new constraint: max_term >= ptemp
			diseq = max_term - ptemp;
			// store the constraint only if it's not a tautology
			if(!diseq.isTautology()) c_set.push_back(diseq);
		}
	}
	else {
		// this is the case where max(0, X1,..Xn) == 0
#ifdef TNODE_VERBOSE_LUDB
		printf("               case MAX()=0 hit\n");
#endif
		max_term.Zero();
		// add the constraints sigma - sigma_i < 0 for all the stages of the local PI
		for(int i=0; i<p.num_stages; i++)
		{
			ptemp.Zero();
			ptemp.SetB(pflow->burst);
			ptemp = (ptemp - p.stages2[i].sigma) / p.stages2[i].rho;
			diseq = max_term - ptemp;
			// store the constraint only if it's not a tautology
			if(!diseq.isTautology()) c_set.push_back(diseq);
		}
	}

	/*****  now max_term contains the chosen element of the max(0, X1,..,Xn) term  *****/

	// now make max_term = S(i,j) + (sigma - sigma_i) / rho_i
	max_term.SumCoeff(flow, 1.0);

	// compute the delay term of the equivalent service curve: D = PI.D + max() + S(i,j)
	pi_eq->delay = p.delay + max_term; // D = PI.D + max() +S(i,j)

	// compute the leaky-bucket stages of the equivalent service curve
	LeakyBucket2 newstage(nvars);
	for(int i=0; i<p.num_stages; i++) //the number of stages is the same as the PI
	{
		// new_sigma = rho_x * [max() + S(i,j)] - (sigma - sigma_x)
		newstage.sigma = (max_term * p.stages2[i].rho) - pflow->burst + p.stages2[i].sigma;
		// new_rho = rho_x - rho
		newstage.rho = p.stages2[i].rho - pflow->rate;
		// store computed stage
		pi_eq->AddStage(newstage);
	}

#ifdef TNODE_VERBOSE_LUDB
	printf("<<< TNode #%02i: Equivalent service curve: ",flow); pi_eq->Print();
#endif
	return *pi_eq;
}


/*	This method computes the equivalent service curve
	for the current t-node.
	1. Retrieve all the equivalent service curves of the child nodes
	2. Convolve them all and convolve also with the node's PI_C, obtaining the local PI
	3. Compute the equivalent service curve of the local PI and return it
*/
PseudoAffine & TNode::getE(const IndexSet &idx_set)
{
	return *pi_eq;
}


/*	This method returns a reference to a Poly which contains the delay bound
	computed at the current t-node for the given set of indexes
	Input:
		idx_set is an instance of the index set corresponding to the desired combination
		c_set will be filled in with the constraints resulting from the desired combination
	Output:
		a Poly reference to the delay bound expression: h(PI,a)
*/
Poly& TNode::getDelayBound(const IndexSet &idx_set, ConstraintSet &c_set)
{
	PseudoAffine &p = getPI(idx_set, c_set);
	PseudoAffine &pe = computeEquivalent(p, idx_set, c_set);
	//subtract the S(h,k) term at the root node
	pe.delay.SumCoeff(flow, -1.0);
	return pe.delay;
}


// Convert a number into the corresponding IndexSet combination
bool TNode::getCombo(IndexSet &set, uint64 n)
{
	uint q, r = n;
	set.clear();
	int size = index_set.size();
	for(int i=size-1; i>=0; i--)
	{
		q = r % index_set[i];
		set.insert(set.begin(),q);
		r = r / index_set[i];
	}
	return true;
}


/*	Increments the IndexSet and return the position of the highest digit that has been incremented
	'set' is the IndexSet to be incremented
	'index' is the first index position to be incremented (higher positions will be unchanged)
	If 'index' is negative, the increment will be applied to the whole set of indexes normally
	The method returns the position of the highest index incremented, or -1 if a wrap-around has occurred
*/
int TNode::incCombo(IndexSet &set, int index)
{
	int size = index_set.size();
	// if we are not starting from the bottom, zero the last index positions
	if(index < 0) index = size - 1;
	else for(int i=index+1; i<size; i++) set[i] = 0;
	// perform the increment
	while(index >= 0)
	{
		if((++set[index]) < index_set[index]) return index;
		// carry set: increment also the next index
		set[index] = 0;
		index--;
	}
	// signal wrap-around
	return -1;
}


// Returns the total number of combinations corresponding to the given set of indexes
uint64 TNode::getNumCombos()
{
	uint64 num = 1;
	for(int i=0; i<index_set.size(); i++)
		num *= index_set[i];
	return num;
}


/* Returns the index subset corresponding to the Ith child of this node
 * Parameters: n_child = number of set desired (0..num_tnodes-1)
 */
bool TNode::getSubset(IndexSet &subset, int n_child)
{
	return getSubset(index_set, n_child, subset);
}


bool TNode::getSubset(const IndexSet &set, int n_child, IndexSet &subset)
{
	subset.clear();
	uint count,size;

	if((n_child<0) || (n_child >= num_tnodes))
		return false;

	// entry #0 in index_schema is the index for the current node
	count = 0;
	// now sum the number of indexes of the preceding children
	for(int i=0; i<=n_child; i++)
		count += index_schema[i];

	// fetch the number of indexes for the child requested
	size = index_schema[n_child+1];
	// copy the subset of indexes of the child into the new set
	for(int i=0; i<size; i++)
		subset.push_back(set[count+i]);

	return true;
}


/*	Computes the LUDB at the current TNode level, rather than at the Tandem level
	For the parameters, see the comments for computeLUDB_experimental.
	'good_combos' can be NULL, in this case information about good combination will be discarded
*/
double TNode::computeLUDB(Poly &delay, Poly &solution, int nvars, IndexSetVector *good_combos)
{
	uint num_combos = getNumCombos();
	IndexSet iset;
	ConstraintSet constraints;
	int iresult;
	double MinDelayBound = TNODE_LUDB_INFINITY;
	Poly optimum_bak(nvars), delay_bak(nvars);

	for(uint64 combo=0; combo<num_combos; combo++)
	{
		// get the combination of indexes corresponding to this iteration
		getCombo(iset, combo);
		// clear the set of constraints
		constraints.clear();
		// compute the delay bound expression at the root t-node for this iteration
		Poly &delay = getDelayBound(iset, constraints);
		// compute the simplex
		Poly optimum(delay.getNumVariables());

		if(constraints.size()==0 && delay.getNumVariables()==1)
		{
			optimum.SetB(delay.GetB());
			iresult = SIMPLEX_RESULT_OK;
			//printf("TNode::computeLUDB(): trivial simplex skipped.\n");
		}
		else iresult = Simplex::minimum(delay, constraints, flow, optimum);

		if(iresult != SIMPLEX_RESULT_OK) continue;

		// update the least delay bound, if necessary
		if(optimum.GetB() < MinDelayBound)
		{
			MinDelayBound = optimum.GetB();
			optimum_bak = optimum;
			delay_bak = delay;
		}
		if(good_combos != NULL) good_combos->push_back(iset);
	}

	solution = optimum_bak;
	delay = delay_bak;

	return MinDelayBound;
}


/*	This method computes the LUDB at the current t-node by recursively computing it at the children t-nodes and trying
	only the combinations reported by them.
	'delay' will contain the expression of the objective function of the "optimum simplex"
	'solution' will contain the solution which gives the optimum
	'nvars' is the number of variables at the top level (# of flows)
	'good_combos' will contain the IndexSets that have produces good simplexes at the lower levels
	Returns the number of simplexes computed in this sub-tree
*/
uint64 TNode::computeLUDB_experimental(Poly &delay, Poly &solution, int nvars, IndexSetVector &good_combos, int max_good_combos, bool eta)
{
	Poly d(nvars);
	Poly s(nvars);
	Poly sol(nvars);
	Poly optimum(nvars), delay_opt(nvars);
	uint64 num_combos = getNumCombos();
	IndexSet set;
	ConstraintSet constraints;
	double MinDelayBound = TNODE_LUDB_INFINITY;
	uint64 num_simplexes = 0;
	uint num_good_combos = 0;

	IndexSetVector children_combos[num_tnodes];
	int children_num_combos[num_tnodes];
	int children_combo[num_tnodes];

	clock_t t_start = clock();

	// terminate recursion when a t-node has no children
	if(num_tnodes == 0)
	{
		//printf("TNode::computeLUDB_fast(t-node=%i): no child t-nodes, compute local simplex\n",flow);
		computeLUDB(delay, solution, nvars, &good_combos);
		return 1;
	}

	// compute LUDB at all child nodes
	// each child tells the father the set of good combinations to be tried on him
	// the father will try all the indexes obtained by all the possible juxtapositions of the children's good combos,
	// plus its own index. In turn, the father will tell its own grandfather only the combinations that were feasible
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
	printf("TNode::computeLUDB_experimental(t-node=%i): computing simplexes at %i child t-nodes, max_good_combos = %i\n",flow,num_tnodes,max_good_combos);
#endif
	num_combos = index_set[0];
	for(int i=0; i<num_tnodes; i++)
	{
		d.Zero();
		s.Zero();
		num_simplexes += tnodes[i]->computeLUDB_experimental(d,s,nvars, children_combos[i], max_good_combos, eta);
		children_num_combos[i] = children_combos[i].size();
		num_combos *= children_num_combos[i];
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
		printf("TNode%02i.child%02i returned %u good combo(s):\n", flow, i, children_num_combos[i]);
		for(int j=0; j<children_num_combos[i]; j++)
		{
			printf("    %02i:  [ ");
			int size = children_combos[i][j].size();
			for(int k=0; k<size; k++) printf("%i ",children_combos[i][j][k]);
			printf("]\n");
		}
#endif
	}

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
	printf("TNode::computeLUDB_experimental(t-node=%i): computing %Lu simplexes:\n",flow,num_combos);
#endif
	if(num_combos > 1000000) eta = true;
	for(uint64 combo=0; combo<num_combos; combo++)
	{
		// check if it's time to update ETA
		if(eta && ((combo & 0x1FFF) == 0x1000))
		{
			clock_t t_now = clock();
			double t_elapsed = ((double)t_now - (double)t_start) / (double)CLOCKS_PER_SEC;
			printf("%06Lu/%06Lu: ETA = %.0lf s   LUDB = %lf\n", combo, num_combos, ((num_combos-combo)*t_elapsed)/combo, MinDelayBound);
		}

		// compute indexes for children
		uint64 q, r;
		r = combo;
		for(int i=num_tnodes-1; i>=0; i--)
		{
			q = r % children_num_combos[i];
			children_combo[i] = q;
			r = r / children_num_combos[i];
		}
		q = r % index_set[0];

		// assemble top-level combination
		set.clear();
		set.push_back(q);
		for(int i=0; i<num_tnodes; i++)
		{
			int v = children_combo[i];
			IndexSet cs = (children_combos[i])[v];
			int size = cs.size();
			for(int j=0; j<size; j++)
			{
				set.push_back(cs[j]);
			}
		}

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
		printf("TNode %02i:  ",flow); PrintIndexSet(set,true);
#endif

		// compute top-level simplex
		constraints.clear();
		Poly &del = getDelayBound(set, constraints);
		int iresult = Simplex::minimum(del, constraints, flow, sol);
#ifdef TNODE_DEBUG_LUDB_CONSTRAINTS
		printf(" --> %i constraints:\n", constraints.size());
		for(int con=0; con<constraints.size(); con++)
		{
			printf("C%03i: ",con); constraints[con].Print(); printf("\n");
		}
		printf("Obj.Fun: "); del.Print(); printf("\n");
		printf("Result");
#endif
		if(iresult != SIMPLEX_RESULT_OK)
		{
			// simplex was unfeasible
#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
			printf(" --> N/F\n");
#endif
			continue;
		}
		if(max_good_combos != LUDB_MODE_EXACT)
		{
			// heuristic mode
			if(sol.GetB() < (MinDelayBound - TNODE_LUDB_EPSILON))
			{
				MinDelayBound = sol.GetB();
				optimum = sol;
				delay_opt = del;
				good_combos.clear();
				good_combos.push_back(set);
				num_good_combos = 1;
#ifndef TNODE_DEBUG_LUDB_EXPERIMENTAL
				//if(eta) printf("TNode %02i:  ",flow); PrintIndexSet(set,true); printf(" --> new LUDB = %lf\n",MinDelayBound);
#endif
			}
			else if(fabs(sol.GetB() - MinDelayBound) <= TNODE_LUDB_EPSILON)
			{
				good_combos.push_back(set);
				num_good_combos++;
			}
		}
		else {
			// exact mode
			if(sol.GetB() < MinDelayBound)
			{
				MinDelayBound = sol.GetB();
				optimum = sol;
				delay_opt = del;
#ifndef TNODE_DEBUG_LUDB_EXPERIMENTAL
				//if(eta) printf("TNode %02i:  ",flow); PrintIndexSet(set,true); printf(" --> new LUDB = %lf\n",MinDelayBound);
#endif
			}
			good_combos.push_back(set);
		}

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
		printf(" --> %lf\n",sol.GetB());
#endif

	}

	solution = optimum;
	delay = delay_opt;
	num_simplexes += num_combos;

	// in heuristic mode, trim the set of good combinations to the desired size
	if(max_good_combos != LUDB_MODE_EXACT)
	{
		// delete exceeding combinations randomly
		while(num_good_combos > max_good_combos)
		{
			int rnd_combo = floor(rng_ludb.uniform(0.0, num_good_combos));
			if(rnd_combo >= num_good_combos) {printf("TNode::computeLUDB_experimental(): RNG right limit trespassed!\n"); rnd_combo = num_good_combos-1;}
			good_combos.erase(good_combos.begin() + rnd_combo);
			num_good_combos--;
		}
	}

#ifdef TNODE_DEBUG_LUDB_EXPERIMENTAL
	printf("TNode::computeLUDB_experimental(): finished, %u simplexes computed in this sub-tree\n",num_simplexes);
#endif
	return num_simplexes;
}


}
