/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   				   *
 *   luca.bisti@iet.unipi.it                                               *
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
#include "tandem.h"
#include "simplex.h"
#include "timing.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <set>
#include <algorithm>


// The following directive enables sanity checks in the LowerBound functions (recommended for testing)
#define TANDEM_DEBUG_CHECKS

// The following directive enables very verbose debug messages during the LowerBound computation
//#define TANDEM_DEBUG_LOWERBOUND

// The following directive enables debugging information in the LowerBound compute_CDF_I() method
//#define TANDEM_DEBUG_CDF_I

// The following directive enables debugging in LowerBound combinations cutting code
//#define TANDEM_DEBUG_LBOPT

// The directive below enables the printing of some info about each iteration during the LUDB computation
//#define TANDEM_VERBOSE_LUDB

// The directive below enables very verbose debug messages during the LUDB computation
//#define TANDEM_DEBUG_LUDB

// The directive below enables debug prints about the tandem cutting functions
//#define TANDEM_DEBUG_CUTS

// The directive below enables verbose debug prints about the output arrival curves computing functions
//#define TANDEM_DEBUG_OUTPUT_ARR_CURVES

//#define TANDEM_DEBUG_TIMINGS


#ifdef TANDEM_CATCH_EXCEPTIONS
#include <exception>
#endif



namespace deborah {

#define _max(x,y) ((x)>(y)?(x):(y))
#define _min(x,y) ((x)<(y)?(x):(y))

/* ************************************************************************************************************* */
/* ********************                     General purpose methods                         ******************** */
/* ************************************************************************************************************* */


Tandem::Tandem()
{
	num_nodes = 0;
	num_flows = 0;
	nesting_level = 0;
	nodes_arr = NULL;
	flows_arr = NULL;
	tagged_flow=0;
	nesting_tree = NULL;
	cached_ludb = -1.0;
}


Tandem::Tandem(uint nodes, uint flows)
{
	num_nodes = nodes;
	num_flows = flows;
	nesting_level = 0;
	nodes_arr = new Node [num_nodes];
	flows_arr = new Flow [num_flows];
	tagged_flow=0;
	nesting_tree = NULL;
	cached_ludb = -1.0;
}

Tandem::Tandem(const Tandem &t2)
{
	nesting_level = 0;
	nesting_tree = NULL;
	num_nodes = t2.num_nodes;
	num_flows = t2.num_flows;
	tagged_flow = t2.tagged_flow;
	cached_ludb = t2.cached_ludb;

	if(num_nodes > 0) {
		nodes_arr = new Node [num_nodes];
		for(uint i=0; i<num_nodes; i++)
		{
			nodes_arr[i] = t2.nodes_arr[i];
		}
	} else nodes_arr = NULL;

	if(num_flows > 0) {
		flows_arr = new Flow [num_flows];
		for(uint i=0; i<num_flows; i++)
		{
			flows_arr[i] = t2.flows_arr[i];
		}
	} else flows_arr = NULL;
}


Tandem::~Tandem()
{
	Clear();
}


uint Tandem::NumNodes()
{
	return num_nodes;
}


uint Tandem::NumFlows()
{
	return num_flows;
}


void Tandem::Clear()
{
	if(nesting_tree != NULL) delete nesting_tree;
	nesting_tree = NULL;
	if(nodes_arr != NULL) delete [] nodes_arr;
	if(flows_arr != NULL) delete [] flows_arr;
	nodes_arr = NULL;
	flows_arr = NULL;
	num_nodes = 0;
	num_flows = 0;
	nesting_level = 0;
	cached_ludb = -1.0;
}


Node* Tandem::getNode(uint index)
{
	if(index < num_nodes) return &nodes_arr[index];
	else return NULL;
}


// Return flow at #index
Flow* Tandem::getFlow(uint index)
{
	if(index < num_flows) return &flows_arr[index];
	else return NULL;
}


// Find a flow by UID
Flow* Tandem::getFlowId(uint id)
{
	for(uint i=0; i<num_flows; i++)
	{
		if(flows_arr[i].uid == id) return &flows_arr[i];
	}
	// not found
	return NULL;
}


// Find a flow by UID
uint Tandem::getFlowIdx(uint id)
{
	for(uint i=0; i<num_flows; i++)
	{
		if(flows_arr[i].uid == id) return i;
	}
	// not found
	return FLOW_INVALID_ID;
}

/*	This method creates a balanced-tree topology.
	'bt.order' is the number of children for each t-node
	'bt.depth' is the number of nesting levels
	If bt.order is 1, then a sink-tree is generated.
*/
bool Tandem::CreateTree(const balanced_tree_descriptor &bt, RNG &rng)
{
	uint curr_flows;
	uint curr_node, nodes_span;
	bool make_sink_tree = false;

	if((bt.depth <= 0) || (bt.order <= 0))
		return false;

	Clear();

	if(bt.order < 2) make_sink_tree = true;

	if(!make_sink_tree)
	{
		// balanced tree
		num_nodes = (uint) pow(bt.order, bt.depth);;
		num_flows = 1;
		for(int i=1; i<=bt.depth; i++) num_flows += (uint) pow(bt.order,i);
	}
	else {
		// sink tree
		num_nodes = bt.depth;
		num_flows = num_nodes;
	}

	nodes_arr = new Node [num_nodes];
	if(nodes_arr == NULL) return false;
	flows_arr = new Flow [num_flows];
	if(flows_arr == NULL)
	{
		Clear();
		return false;
	}

	// create flows
	if(make_sink_tree)
	{
		// create sink tree
		for(unsigned int i=0; i<num_flows; i++)
		{
			flows_arr[i].src = i;
			flows_arr[i].exit = num_flows-1;
			flows_arr[i].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
			flows_arr[i].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
			flows_arr[i].uid = i;
		}
	}
	else
	{
		// create balanced tree
		num_flows = 0;
		for(int i=0; i<=bt.depth; i++)
		{
			curr_flows = (uint) pow(bt.order, i); // number of tnodes at this level
			nodes_span = num_nodes / curr_flows;
			curr_node = 0;
			for(uint j=0; j<curr_flows; j++)
			{
				flows_arr[num_flows].src = curr_node;
				curr_node += nodes_span;
				flows_arr[num_flows].exit = curr_node-1;
				flows_arr[num_flows].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
				flows_arr[num_flows].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
				flows_arr[num_flows].uid = num_flows;
				num_flows++;
			}
		}
	}

	// create nodes
	for(unsigned int i=0; i<num_nodes; i++)
	{
		double node_rate_min = 0.0;
		for(unsigned int j=0; j<num_flows; j++)
		{
			if(isNodeCrossedByFlow(i,j)) node_rate_min += flows_arr[j].rate;
		}
		nodes_arr[i].latency = rng.uniform(bt.node_latency_min, bt.node_latency_max);
		nodes_arr[i].rate = node_rate_min * rng.uniform(bt.node_rate_overprov_min, bt.node_rate_overprov_max);
	}

	tagged_flow = 0;

	ComputeNestingLevels();

	return true;
}

/*	This method creates a tandem comprising 'n' nodes and 'k' flows
 *  for all possible combinations of (src, dst) nodes.
 *  Note: the resulting tandem may be non-nested.
 */
bool Tandem::CreateFullTandem(const balanced_tree_descriptor &bt, RNG &rng)
{
	bool random = false;
	Clear();
	if(bt.order < 3) return false;

	// allocate nodes
	num_nodes = bt.order;
	nodes_arr = new Node [num_nodes];

	// allocate flows
	num_flows = num_nodes*(num_nodes-1)/2;
	if(bt.depth > 0 && bt.depth < 100)
	{
		num_flows *= (float)bt.depth / 100.0;
		random = true;
	}
	else if(bt.depth < 0 && bt.depth > -num_flows)
	{
		num_flows = -bt.depth;
		random = false;
	}
	flows_arr = new Flow [num_flows];

	// configure flows
	uint fid = 0;

	if(random) while (fid<num_flows)
	{	// only a subset of flows
		if(fid==0)
		{
			// set tagged flow manually
			flows_arr[fid].src=0;
			flows_arr[fid].exit=num_nodes-1;
		}
		else {
			flows_arr[fid].src = rng.uniform(0, num_nodes-1);
			flows_arr[fid].exit = rng.uniform(flows_arr[fid].src+1, num_nodes);
			if(flows_arr[fid].src == flows_arr[fid].exit) continue;
			bool found = false;
			for(uint i=0; i<fid; i++)
				if(flows_arr[i].src==flows_arr[fid].src && flows_arr[i].exit==flows_arr[fid].exit)
				{
					found = true;
					break;
				}
			if(found) continue;
		}
		flows_arr[fid].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
		flows_arr[fid].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
		flows_arr[fid].uid = fid;
		fid++;
	}
	else for(int flen=num_nodes; flen>=2; flen--)
	{	// all possible flows
		int nf = num_nodes-flen+1;
		for(int i=0; i<nf; i++)
		{
			if(fid<num_flows)
			{
				flows_arr[fid].src = i;
				flows_arr[fid].exit = i+flen-1;
				flows_arr[fid].burst = rng.uniform(bt.flow_burst_min, bt.flow_burst_max);
				flows_arr[fid].rate = rng.uniform(bt.flow_rate_min, bt.flow_rate_max);
				flows_arr[fid].uid = fid;
			}
			fid++;
		}
	}

	tagged_flow = 0;
	// sanity check
	if(fid != num_flows)
	{
		printf("Tandem::CreateFullTandem(): BUG: %i flows created, expected count = %i\n",fid,num_flows);
	}

	// configure nodes
	for(unsigned int i=0; i<num_nodes; i++)
	{
		double node_rate_min = 0.0;
		for(unsigned int j=0; j<num_flows; j++)
		{
			if(isNodeCrossedByFlow(i,j)) node_rate_min += flows_arr[j].rate;
		}
		nodes_arr[i].latency = rng.uniform(bt.node_latency_min, bt.node_latency_max);
		nodes_arr[i].rate = node_rate_min * rng.uniform(bt.node_rate_overprov_min, bt.node_rate_overprov_max);
	}

	return true;
}


/*	This method sets the current object according to the configuration file
 *  In case of error, it returns an empy object
 */
bool Tandem::Load(char *config_filename)
{
	FILE *fconf;
	char buff[1024];
	char directive[1024];
	bool result;
	int tagged=-1;

	int nodes=0, flows=0;

	fconf = fopen(config_filename,"rt");
	if(fconf == NULL)
	{
		return false;
	}

	//Delete previous data, if present
	Clear();

	while(true)
	{
		uint src_node, exit_node;
		double theta,R,sigma,rho;

		buff[0]=0;
		if(fgets(buff,1024,fconf) == 0) break;
		if(sscanf(buff,"%s",directive) <= 0) continue;
		if(directive[0] == '#') continue;
		if(!strcasecmp(directive,"TANDEM"))
		{
			sscanf(&buff[6],"%i %i",&num_nodes, &num_flows);
			if((num_nodes<=0) || (num_flows<=0)) break;
			//Allocate structures
			nodes_arr = new Node [num_nodes];
			flows_arr = new Flow [num_flows];
			continue;
		}
		if(!strcasecmp(directive,"NODE"))
		{
			uint node_id;
			if(nodes >= num_nodes) continue;
			theta=0.0; R=0.0;
			sscanf(&buff[5],"%i %lf %lf",&node_id,&theta,&R);
			//printf("NODE %i: %lf, %lf\n",node_id,theta,R);
			if(node_id-1 > num_nodes)
			{
				printf("Warning: invalid node ID %02i in configuration file, skipped",node_id);
				continue;
			}
			//Add node to tandem
			nodes_arr[node_id-1].latency = theta;
			nodes_arr[node_id-1].rate = R;
			nodes++;
			continue;
		}
		if(!strcasecmp(directive,"FLOW") || !strcasecmp(directive,"TFLOW"))
		{
			if(flows >= num_flows) continue;
			// check for tagged-flow directive
			if((directive[0]=='T') || (directive[0]=='t'))
				tagged=flows;
			src_node=-1; exit_node=-1;
			sigma=0.0; rho=0.0;
			sscanf(&buff[5],"%i %i %lf %lf", &src_node, &exit_node, &sigma, &rho);
			//printf("FLOW (%i,%i): %lf, %lf\n",src_node,exit_node,sigma,rho);
			if((src_node == 1) && (exit_node == num_flows)) tagged_flow = flows;
			//Add flow to tandem
			flows_arr[flows].src = src_node-1;
			flows_arr[flows].exit = exit_node-1;
			flows_arr[flows].burst = sigma;
			flows_arr[flows].rate = rho;
			flows_arr[flows].uid = flows;
			flows++;
			continue;
		}
	}

	//Close configuration file
	fclose(fconf);

	//Check results
	result = true;
	if((nodes<=0) || (flows<=0)) result = false;
	if((nodes<num_nodes) || (flows<num_flows))
	{
		printf("Error: config file is not consistent: missing declarations for nodes and/or flows\n");
		result = false;
	}

	if(tagged >= 0) tagged_flow = tagged;

	//Now update the auxiliary data structures (nesting tree, etc)
	ComputeNestingLevels();

	//In case of error, empty the current Tandem
	if(!result) Clear();

	return result;
}


/*	Writes the current tandem topology to a configuration file
*/
bool Tandem::Save(char *config_filename)
{
	FILE *config;

	// some preliminary checks
	if(config_filename == NULL) return false;
	if(!num_nodes || !num_flows) return false;

	config = fopen(config_filename, "wt");
	if(config == NULL) return false;

	// write header
	fprintf(config, "# DEBORAH CONFIGURATION FILE v1.2\n\n");
	fprintf(config, "# TANDEM <nodes #> <flows #>\n");
	fprintf(config, "TANDEM %i %i\n\n", num_nodes, num_flows);

	// write nodes
	fprintf(config,"# NODE <node #> <latency> <rate>\n");
	for(int i=0; i<num_nodes; i++)
	{
		fprintf(config, "NODE	%i	%lf	%lf\n", i+1, nodes_arr[i].latency, nodes_arr[i].rate);
	}
	fprintf(config,"\n");

	// write flows
	fprintf(config,"# FLOW <source node> <sink node> <burst> <rate>\n");
	for(int i=0; i<num_flows; i++)
	{
		fprintf(config, "FLOW	%u	%u	%lf	%lf\n", flows_arr[i].src+1, flows_arr[i].exit+1, flows_arr[i].burst, flows_arr[i].rate);
	}
	fprintf(config,"\n");

	fclose(config);
	return true;
}


void Tandem::Print()
{
	printf("------------------------------------------------------------\n");
	printf("Tandem: %i nodes, %i flows, nesting_level=%i\n",num_nodes,num_flows, nesting_level);
	if(num_nodes>0 && num_flows>0)
	{
		printf("        Tagged flow = (%i,%i)  flow_id=#%i\n",flows_arr[tagged_flow].src+1,flows_arr[tagged_flow].exit+1,tagged_flow);
		printf("        Tandem is %snested.\n", (isNested()?"":"not "));
		for(int i=0;i<num_nodes;i++) printf("NODE %02i: latency=%.2lf, rate=%.2lf\n",i+1,nodes_arr[i].latency,nodes_arr[i].rate);
		for(int i=0;i<num_flows;i++) printf("FLOW %02i (%i,%i): burst=%.2lf, rate=%.2lf; nesting_level=%i; uid=%x\n",i,flows_arr[i].src+1,flows_arr[i].exit+1,flows_arr[i].burst,flows_arr[i].rate,flows_arr[i].nesting_level,flows_arr[i].uid);
	}
	printf("------------------------------------------------------------\n");

}


bool Tandem::isSinkTree()
{
	if(num_flows != num_nodes) return false;
	for(int i=0; i<num_flows; i++)
	{
		if((flows_arr[i].src != i) || (flows_arr[i].exit != num_nodes-1))
			return false;
	}
	return true;
}


bool Tandem::isNested()
{
	// check that no two flows intersect each other
	for(unsigned int i=0; i<num_flows; i++)
	{
		for(unsigned int j=0; j<num_flows; j++)
		{
			if(i == j) continue;
			if(flows_arr[i].intersects(flows_arr[j])) return false;
		}
	}
	return true;
}


bool Tandem::checkProvisioning()
{
	double rate;
	for(unsigned int i=0; i<num_nodes; i++)
	{
		rate = 0.0;
		for(unsigned int j=0; j<num_flows;j++)
			if(isNodeCrossedByFlow(i,j))
				rate += flows_arr[j].rate;
		if(rate > nodes_arr[i].rate) return false;
	}
	return true;
}


bool Tandem::scaleProvisioning(double flow_mult, double node_mult)
{
	if((flow_mult <= 0.0) || (node_mult <= 0.0))
		return false;
	for(unsigned int i=0; i<num_nodes; i++)
	{
		nodes_arr[i].rate *= node_mult;
		for(unsigned int j=0; j<num_flows; j++)
			if(isNodeCrossedByFlow(i,j))
				flows_arr[j].rate *= flow_mult;
	}
	return checkProvisioning();
}


bool Tandem::analyzeProvisioning(double &min_prov, int &nodemin, double &max_prov, int &nodemax)
{
	if(!num_nodes || !num_flows) return false;
	double rate, overprov;
	min_prov = 10000000.0; max_prov = 0.0;
	nodemin = -1; nodemax = -1;
	for(unsigned int i=0; i<num_nodes; i++)
	{
		rate = 0.0;
		for(unsigned int j=0; j<num_flows;j++)
			if(isNodeCrossedByFlow(i,j))
				rate += flows_arr[j].rate;
		overprov = nodes_arr[i].rate / rate;
		if(overprov < min_prov) {
			min_prov = overprov;
			nodemin = i;
		}
		if(overprov > max_prov) {
			max_prov = overprov;
			nodemax = i;
		}
	}
	return true;
}



/* ************************************************************************************************************* */
/* ********************                       LUDB related methods                          ******************** */
/* ************************************************************************************************************* */




/* This method computes the nesting level of each flow in the tandem
 * For simplicity, values are stored in each flow object
 */
void Tandem::ComputeNestingLevels()
{
	if((num_flows <=0) || (flows_arr==NULL)) return;

	nesting_level = 0;
	for(int i=0; i<num_flows; i++)
	{
		flows_arr[i].nesting_level = 0;
		for(int j=0 ;j<num_flows; j++)
		{
			if(flows_arr[i].is_nested(flows_arr[j]))
			{
				flows_arr[i].nesting_level++;
				if(flows_arr[i].nesting_level > nesting_level) nesting_level = flows_arr[i].nesting_level;
			}
		}
	}
}


/*	This method builds the complete nesting tree for this tandem
 */
bool Tandem::BuildNestingTree()
{
	if(flows_arr==NULL) return false;
	if(nodes_arr==NULL) return false;

	if(nesting_tree != NULL) delete nesting_tree;

	nesting_tree = build_nesting_subtree(tagged_flow);
	if(nesting_tree == NULL) return false;

	return true;
}


/*	Auxiliary method which builds a sub-tree
*/
TNode* Tandem::build_nesting_subtree(uint flow_id)
{
	TNode *tnode =  NULL;
	uint n_leaves;
	uint n_tnodes;
	int nest_level;
	bool contained;

	//pointer to considered flow (h,k)
	Flow *f;
	//pointer to flow(s) directly nested into (h,k)
	Flow *nf;
	//temp pointer to a node
	Node *pn;

	if(flow_id >= num_flows) return NULL;
	f = &flows_arr[flow_id];
	nest_level = f->nesting_level;

	tnode = new TNode();
	tnode->flow = flow_id;
	tnode->pflow = f;

	//Find number of directly nested flows
	n_tnodes = 0;
	for(int i=0; i<num_flows; i++)
	{
		if(flows_arr[i].nesting_level == (nest_level+1))
		{
			if(flows_arr[i].is_nested(*f))
			{
				n_tnodes++;
			}
		}
	}
	tnode->num_tnodes = n_tnodes;

	//Create child nodes for directly nested flows
	if(tnode->num_tnodes > 0)
	{
		tnode->tnodes = new TNode* [tnode->num_tnodes];
		n_tnodes = 0;
		for(int i=0; i<num_flows; i++)
		{
			if(flows_arr[i].nesting_level == (nest_level+1))
			{
				if(flows_arr[i].is_nested(*f))
				{
					tnode->tnodes[n_tnodes] = build_nesting_subtree(i);
					//check errors in the costruction of the subtree
					if(tnode->tnodes[n_tnodes] == NULL) return NULL;
					n_tnodes++;
				}
			}
		}
	}
	else tnode->tnodes = NULL;

	//Find number of leaves
	n_leaves=0;
	for(int i=f->src; i<=f->exit; i++)
	{
		contained = false;
		for(int j=0;j<tnode->num_tnodes;j++)
		{
			nf = &flows_arr[tnode->tnodes[j]->flow];
			if(nf->contains(i)) contained = true;
		}
		if(!contained) n_leaves++;
	}
	tnode->num_leaves = n_leaves;

	//Create the actual leaves
	if(tnode->num_leaves > 0)
	{
		tnode->leaves = new uint[tnode->num_leaves];
		n_leaves = 0;
		for(int i=f->src; i<=f->exit; i++)
		{
			contained = false;
			for(int j=0;j<tnode->num_tnodes;j++)
			{
				nf = &flows_arr[tnode->tnodes[j]->flow];
				if(nf->contains(i)) contained = true;
			}
			if(!contained)
			{
				tnode->leaves[n_leaves] = i;
				n_leaves++;
			}
		}

	}

	// Compute service curve PI_C(h,k) for leaves
	if(tnode->has_leaf())
	{
		pn = getNode(tnode->leaves[0]);
		tnode->pi_c.setLatencyRate(pn->latency, pn->rate);
		for(int i=1; i<tnode->num_leaves; i++)
		{
			pn = getNode(tnode->leaves[i]);
			tnode->pi_c.convolve(pn->latency, pn->rate);
		}
	}

	// Initialize service curves PI(h,k) and PI_EQ(h,k) for the current t-node
	if(tnode->has_leaf())
	{
		// Initialize the service curve at this "child flow" with the PI_C
		tnode->pi = new PseudoAffine(num_flows);
		tnode->pi_eq = new PseudoAffine(num_flows);
		tnode->pi->MakeLatencyRate(tnode->pi_c);
		//tnode->pi->Print();
	}
	else {
		// Initialize an empty PI and PI_EQ
		tnode->pi = new PseudoAffine(num_flows);
		tnode->pi_eq = new PseudoAffine(num_flows);
	}
	// When you compute PI_EQ you introduce new indexes for the current node!

	// Build the set of sets of indexes of the current t-node
	tnode->index_set.clear();
	tnode->index_schema.clear();
	tnode->index_nodes.clear(); tnode->index_nodes.push_back(tnode->flow);
	// first add the (single) index for our own t-node
	tnode->index_schema.push_back(1);
	// initialize the index limit as the number of stages of the equivalent service curve of this t-node
	uint index_limit = tnode->has_leaf() ? 1 : 0;
	// then add the sets of indexes from every child t-node
	for(int i=0; i<tnode->num_tnodes; i++)
	{
		std::vector<uint> child_schema = tnode->tnodes[i]->index_schema;
		IndexSet child_set = tnode->tnodes[i]->index_set;
		for( int j=0; j<child_set.size(); j++)
		{
			tnode->index_set.push_back(child_set[j]);
			tnode->index_nodes.push_back(tnode->tnodes[i]->index_nodes[j]);
		}
		tnode->index_schema.push_back(child_set.size());
		// add the number of stages of the child's equivalent service curve
		// note: if the child has no leaf, then the actual number of stages of his PI is
		// child_set[0]-1 because its index limit takes into account the extra case "max()=0"
		if(tnode->tnodes[i]->has_leaf())
			index_limit += child_set[0];
		else index_limit += child_set[0] - 1;
	}
	// if we have no leaf, count ONE MORE choice representing "case max()=0"
	if(!tnode->has_leaf()) index_limit++;
	// head-insert the determined index limit for our t-node
	tnode->index_set.insert(tnode->index_set.begin(), index_limit);

	return tnode;
}


void Tandem::PrintNestingTree()
{
	if(nesting_tree != NULL) nesting_tree->Print(true);
	else printf("Nesting Tree not computed for this tandem!\n");
}


uint64 Tandem::LUDB_TotalCombos()
{
	if(nesting_tree == NULL) return 0;
	return nesting_tree->getNumCombos();
}


void Tandem::PrintAllCombinations()
{
	if(nesting_tree == NULL) return;

	IndexSet set;
	uint64 num = nesting_tree->getNumCombos();
	for(uint i=0; i<num; i++)
	{
		nesting_tree->getCombo(set,i);
		printf("%03i: [ ",i);
		for(int j=0;j<set.size();j++) printf("%i ",set[j]);
		printf("]\n");
	}

	// alternate enumeration method
	printf("--- alternate method\n");
	nesting_tree->getCombo(set, 0);
	for(uint i=0; i<num; i++)
	{
		int h;
		printf("%03i: [ ",i);
		for(int j=0;j<set.size();j++) printf("%i ",set[j]);
		printf("]\n");
		h = nesting_tree->incCombo(set, -1);
	}
}


void Tandem::Test_LUDB_code()
{
	uint num_combos = nesting_tree->getNumCombos();
	uint combo = 0;
	IndexSet set;
	bool result;

	ConstraintSet constraints;

	printf("Tandem::Test()\n");
	printf("----------------------------------------------------\n");

	nesting_tree->getCombo(set, combo);
	constraints.clear();
	Poly &delay = nesting_tree->getDelayBound(set, constraints);
	int num_constraints = constraints.size();

	printf("\n");
	printf("Combination: %i/%i\n",combo,num_combos);
	printf("*** PI = "); nesting_tree->pi->Print(); printf("\n");
	printf("*** EQ = "); nesting_tree->pi_eq->Print();printf("\n");
	printf("*** Delay = "); delay.Print();printf("\n");
	printf("    Number of associated constraints: %i\n",num_constraints);

	// check if the constraints are feasible
	result = true;
	for(int i=0; i<num_constraints; i++)
	{
		if(constraints[i].isOxymoron())
		{
			result = false;
			printf("Impossible constraint found, skip this simplex!!!\n");
			constraints[i].Print();
			break;
		}
	}
	if(result) printf("Constraints check passed.\n");

	Simplex::test_es4();
}


/*	OBSOLETED! Not used anymore in the code.
	This is the dumb version of the LUDB computation function, which tries ALL simplexes.
	It is expected that the nesting tree has been already computed.
*/
double Tandem::LUDB()
{
	uint num_combos;
	IndexSet set;
	bool result;
	int iresult;

	ConstraintSet constraints;
	int num_constraints;

	double MinDelayBound = LINEARSEGMENT_NEGATIVE_INFINITE;
	Poly optimum(num_flows);

	Poly optimum_bak(num_flows);
	Poly delay_bak(num_flows);
	uint optimum_combo=0;

	clock_t t_start, t_end;
	double t_elapsed;
	uint64 clocks_start, clocks_end, clocks_elapsed;
	uint64 clocks_per_feasible, clocks_per_unfeasible;

	if(nesting_tree == NULL) return TANDEM_LUDB_ERROR_NESTINGTREE;
	num_combos = nesting_tree->getNumCombos();

	printf("Tandem::LUDB() algorithm started:\n");
	printf("\n----------------------------------------------------\n");

	uint num_unfeasible = 0; // unbounded + empty
	uint num_unbounded = 0;
	uint num_avoided = 0;
	uint num_constraints_max = 0;
	uint64 clocks_unfeasible = 0, clocks_feasible = 0;

	t_start = clock();

	for(int combo=0; combo<num_combos; combo++)
	{
		// get the combination of indexes corresponding to this iteration
		nesting_tree->getCombo(set, combo);

		// clear the set of constraints
		constraints.clear();

#ifdef TANDEM_DEBUG_LUDB
		printf("Iteration %03i/%03i: ",combo,num_combos);
		nesting_tree->PrintIndexSet(set, false);
		printf("\n");
#endif

		// compute the delay bound expression at the root t-node for this iteration
		Poly &delay = nesting_tree->getDelayBound(set, constraints);
		// get the number of constraints
		num_constraints = constraints.size();
		if(num_constraints > num_constraints_max) num_constraints_max = num_constraints;

#ifdef TANDEM_DEBUG_LUDB
		printf("Delay bound expression: "); delay.Print(); printf("    %i constraint(s), dead variable=X%i\n",num_constraints,nesting_tree->flow);
#endif

		// sanity check for the constraints
		result = true;
		for(int i=0; i<num_constraints; i++)
		{
			if(constraints[i].isOxymoron())
			{
				result = false;
				printf("[!] constraints sanity check failed, skip this simplex (%i)!!!\n",combo);
				printf("    offending constraint: "); constraints[i].Print(); printf("\n");
				break;
			}
#ifdef TANDEM_DEBUG_LUDB
			printf("C #%03i: ",i); constraints[i].Print(); printf("\n");
#endif
		}
		// if the sanity check has failed, go to the next iteration
		if(!result)
		{
			printf("Constraints sanity check failed, skipping simplex\n");
			printf("----------------------------------------------------\n");
			num_avoided++;
			continue;
		}

		// perform the simplex
		READTSC(&clocks_start);
		iresult = Simplex::minimum(delay, constraints, nesting_tree->flow, optimum);
		READTSC(&clocks_end);
		clocks_elapsed = clocks_end - clocks_start;

#ifdef TANDEM_VERBOSE_LUDB
		printf("%03i/%03i: ",combo,num_combos);
		nesting_tree->PrintIndexSet(set,true);
		if(iresult == SIMPLEX_RESULT_OK) printf(" --> %lf\n", optimum.GetB());
		else printf(" --> N/F\n");
#endif

		// did we find a solution?
		if(iresult != SIMPLEX_RESULT_OK)
		{
#ifdef TANDEM_DEBUG_LUDB
			printf("Simplex was unfeasible\n");
			printf("----------------------------------------------------\n");
#endif
			num_unfeasible++;
			if(iresult == SIMPLEX_RESULT_UNBOUNDED) num_unbounded++;
			clocks_unfeasible += clocks_elapsed;
			continue;
		}
		else clocks_feasible += clocks_elapsed;
#ifdef TANDEM_DEBUG_LUDB
		printf("Simplex output: "); optimum.Print(); printf("   (optimum delay = %lf)\n",optimum.GetB());
		printf("----------------------------------------------------\n");
#endif

		// update the least delay bound, if necessary
		if(optimum.GetB() < MinDelayBound)
		{
			MinDelayBound = optimum.GetB();
			optimum_bak = optimum;
			optimum_combo = combo;
			delay_bak = delay;
			nesting_tree->PrintIndexSet(set,true); printf(" --> New LUDB = %lf\n", optimum.GetB());
		}
	}

	t_end = clock();
	t_elapsed = ((double)t_end - (double)t_start) / (double)CLOCKS_PER_SEC;
	if(num_combos != num_unfeasible) clocks_per_feasible = clocks_feasible/(num_combos-num_unfeasible);
	else clocks_per_feasible = 0;
	if(num_unfeasible > 0) clocks_per_unfeasible = clocks_unfeasible/num_unfeasible;
	else clocks_per_unfeasible = 0;

	nesting_tree->getCombo(set, optimum_combo);
	printf("LUDB = %lf *** found at iteration %03i/%03i: ",MinDelayBound,optimum_combo,num_combos);
	nesting_tree->PrintIndexSet(set, false);
	printf("\n");
	printf("Delay bound expression: "); delay_bak.Print(); printf("\n");
	printf("Simplexes: computed %i  /  avoided %i  /  unfeasible %i (unbounded %i)\n",num_combos-num_avoided,num_avoided,num_unfeasible, num_unbounded);
	printf("Number of constraints per simplex: %i (max)\n",num_constraints_max);
	printf("Computation time: %.3lf seconds\n", t_elapsed);
	printf("Average CPU clocks per simplex: feasible = ");
	if(clocks_per_feasible > 0)
		printf("%qd",clocks_per_feasible);
	else printf("N/A");
	printf(", unfeasible = ");
	if(clocks_per_unfeasible > 0)
		printf("%qd",clocks_per_unfeasible);
	else printf("N/A");
	printf("\n");
	printf("Solution:\n"); PrintSolution(optimum_bak,tagged_flow); printf("\n");

	cached_ludb = MinDelayBound;
	return MinDelayBound;
}



double Tandem::LUDB_experimental(int max_good_combos)
{
	ConstraintSet constraints;
	uint64 n_simplexes;
	Poly solution(num_flows);
	Poly delay(num_flows);
	IndexSetVector good_combos;
	clock_t t_start, t_end;

	if(nesting_tree == NULL) return -1.0;
	if(num_flows > NMAX) return -2.0;

	printf("%s LUDB algorithm started:\n", ((max_good_combos != LUDB_MODE_EXACT) ? "Heuristic":"Exact"));
	if(max_good_combos != LUDB_MODE_EXACT) printf("Maximum number of combinations returned by each t-node: %i\n",max_good_combos);
	printf("\n----------------------------------------------------\n");

	double t_elapsed;
	t_start = clock();
	n_simplexes = nesting_tree->computeLUDB_experimental(delay, solution, num_flows, good_combos, max_good_combos, true);
	t_end = clock();
	t_elapsed = ((double)t_end - (double)t_start) / (double)CLOCKS_PER_SEC;

	printf("LUDB = %lf\n",solution.GetB());
	printf("Computation time: %.3lf seconds\n", t_elapsed);
	printf("Simplexes computed: %Lu\n", n_simplexes);
	printf("Delay bound expression: "); delay.Print(); printf("\n");
	printf("Solution:\n"); PrintSolution(solution,tagged_flow); printf("\n");
	cached_ludb = solution.GetB();
	return solution.GetB();
}


double Tandem::LUDB_experimental_quiet(int max_good_combos, uint64 &n_simplexes)
{

	return LUDB_experimental_quiet(nesting_tree, max_good_combos, n_simplexes);
	/*
	ConstraintSet constraints;
	Poly solution(num_flows);
	Poly delay(num_flows);
	IndexSetVector good_combos;

	if(nesting_tree == NULL) return TANDEM_LUDB_ERROR_NESTINGTREE;

	n_simplexes = nesting_tree->computeLUDB_experimental(delay, solution, num_flows, good_combos, max_good_combos, false);

	return solution.GetB();
	*/
}


double Tandem::LUDB_experimental_quiet(TNode *root, int max_good_combos, uint64 &n_simplexes)
{
	ConstraintSet constraints;
	Poly solution(num_flows);
	Poly delay(num_flows);
	IndexSetVector good_combos;

	if(root == NULL) return TANDEM_LUDB_ERROR_NESTINGTREE;

	//printf("root=%x   max good combos = %i\n",root,max_good_combos);
	n_simplexes = root->computeLUDB_experimental(delay, solution, num_flows, good_combos, max_good_combos, false);

	cached_ludb = solution.GetB();
	return solution.GetB();
}


/*	This method prints a solution vector for the LUDB algorithm
*/
void Tandem::PrintSolution(Poly &p, uint dead_var)
{
	char stmp[32];
	int nvars = p.getNumVariables();

	for(int i=0; i<nvars; i++)
	{
		if(i==dead_var) continue;
		flows_arr[i].to_string(stmp);
		printf("S%s = %lf\n",stmp,p.GetCoeff(i));
	}
}


bool Tandem::SetTaggedFlow(uint tf)
{
	if(tf >= num_flows) return false;
	tagged_flow = tf;
	return true;
}


uint Tandem::GetTaggedFlow()
{
	return tagged_flow;
}


/*	Joins a vector of sub-tandems into the original tandem
 * 	nearly reverting the operation of cut()
 * 	Return: zero or one if the resulting tandem is nested or not
 * 	(negative integer for error)
 */
int Tandem::join(std::vector<Tandem> &cuts)
{
	Clear();

	uint size = cuts.size();
	if(!size) return DEBORAH_ERROR_SIZE;
	for(uint i=0; i<size; i++)
	{
		num_nodes += cuts[i].num_nodes;
		num_flows += cuts[i].num_flows;
	}
	if(!num_nodes || (num_flows <= size-1))
		return DEBORAH_ERROR_SAFECHECK;
	num_flows -= size - 1; // count tagged flow only once
	flows_arr = new Flow [num_flows];
	nodes_arr = new Node [num_nodes];
	if(!flows_arr || !nodes_arr) return DEBORAH_ERROR_ALLOC;
	// reserve flow #0 for tagged flow
	uint curr_flow = 1;
	uint curr_node = 0;
	for(uint i=0; i<size; i++)
	{
		Tandem &ti = cuts[i];
		uint tf = ti.tagged_flow;
		// compile flows
		for(uint j=0; j<ti.num_flows; j++)
		{
			if(j == tf) continue; // skip tagged flow
			flows_arr[curr_flow].burst = ti.flows_arr[j].burst;
			flows_arr[curr_flow].rate = ti.flows_arr[j].rate;
			//flows_arr[curr_flow].uid = ti.flows_arr[j].uid;
			flows_arr[curr_flow].uid = curr_flow;
			flows_arr[curr_flow].src = curr_node + ti.flows_arr[j].src;
			flows_arr[curr_flow].exit = curr_node + ti.flows_arr[j].exit;
			curr_flow++;
		}
		// compile nodes
		for(uint j=0; j<cuts[i].num_nodes; j++)
		{
			nodes_arr[curr_node].latency = cuts[i].nodes_arr[j].latency;
			nodes_arr[curr_node].rate = cuts[i].nodes_arr[j].rate;
			curr_node++;
		}
	}
	// finally set the tagged flow
	tagged_flow = 0;
	Flow *tf = cuts[0].getFlow(cuts[0].tagged_flow);
	flows_arr[tagged_flow].burst = tf->burst;
	flows_arr[tagged_flow].rate = tf->rate;
	flows_arr[tagged_flow].src = tf->src;
	flows_arr[tagged_flow].exit = num_nodes-1;
	//flows_arr[tagged_flow].uid = tf->uid;
	flows_arr[tagged_flow].uid = 0;

	if(isNested())
	{
		ComputeNestingLevels();
		return 1;
	}

	return 0;
}


/* Compute the exhaustive set of cuts for this (non-nested) tandem
 * Cuts longer than 'max_length' nodes will be pruned
 */
std::vector<CutNodelist> Tandem::compute_cutsets()
{
	std::vector<CutNodelist> result;
	// for non-nested tandems only
	if(isNested()) return result;
	std::set<std::vector<uint> > cutpos;

	BasicTimer tcut, tcomp;

	tcut.mark(true);
	// fill-in the conflicts matrix
	for(uint i=0; i<num_flows; i++)
	{
		for(uint j=0; j<num_flows; j++)
		{
			uint s1=flows_arr[i].src;
			uint e1=flows_arr[i].exit;
			uint s2=flows_arr[j].src;
			uint e2=flows_arr[j].exit;
			if(s1<s2 && s2<=e1 && e1<e2)
			{
				//printf("intersection %02i with %02i: ",i,j);
				// flows intersect between nodes [s2,e1] included
				std::vector<uint> alt;
				for(uint k=s2; k<=e1+1; k++)
				{
					alt.push_back(k);
					//printf("%i ",k);
					//node_tag[k]++;
				}
				//printf("\n");
				cutpos.insert(alt);
			}
		}
	}

	/*
	for(int i=0; i<num_deps; i++)
	{
		printf("d%i:  ",i);
		for(int j=0; j<cutpos[i].size(); j++) printf("%i ",cutpos[i][j]);
		printf("\n");
	}
	*/

	// setup the node vs dependency matrix
	uint num_deps = cutpos.size(); // number of dependencies
	std::set<uint> nodedeps[num_nodes];

	// fill-in
	for(uint n=0; n<num_nodes; n++)
	{
		//printf("node %i:  ",n);
		uint d=0;
		for(std::set<std::vector<uint> >::iterator st=cutpos.begin(); st!=cutpos.end(); st++)
		{
			if(std::find(st->begin(), st->end(), n) != st->end())
			{
				nodedeps[n].insert(d);
				//printf("d%i ",d);
			}
			d++;
		}
		//printf("\n");
	}

#ifdef TANDEM_DEBUG_TIMINGS
	tcut.mark(false);
	printf("compute_cutsets(): t_dep=%.3lf  n_deps=%i\n",tcut.time(), num_deps);
#endif

	// compute the sets of cuts
#ifdef TANDEM_CUTCOMP_OLDVERSION
	bool res[num_nodes]; for(uint i=0; i<num_nodes; i++) res[i]=false;
	std::set<uint> cdeps;
#endif

	std::set<std::set<uint> > allres;

	compstatus stat;
	stat.cnode = 0;
	stat.cres = 0x0000;
	if(num_deps > MAXDEPS*32) {
		printf("Tandem::compute_cutsets(): ERROR, max %i dependencies supported (%i found).\n",MAXDEPS*32,num_deps);
		return result;
	}
	depsse alldeps[num_nodes];

	int a = num_deps / MAXDEPS;
	int b = num_deps % MAXDEPS;
	//printf("num_deps=%i   a=%i b=%i\n",num_deps,a,b);
	for(uint i=0; i<num_nodes; i++)
	{
		for(int j=0; j<=a; j++) alldeps[i].dep[j] = 0;
		for(int k=0; k<32-b; k++) {alldeps[i].dep[a] >>= 1; alldeps[i].dep[a] |= 0x80000000;}
		for(int j=a+1; j<MAXDEPS; j++) alldeps[i].dep[j] = 0xFFFFFFFF;
		for(std::set<uint>::iterator jt=nodedeps[i].begin(); jt!=nodedeps[i].end(); jt++)
		{
			uint nd = *jt;
			int a = nd >> 5;
			int b = nd & 0x1F;
			//printf("node %02i: dep %03i --> a=%i b=%i\n",i,nd,a,b);
			alldeps[i].dep[a] |= (1L << b);
		}
		//printf("alldeps[%02i] = ",i); for(int j=0;j<MAXDEPS;j++) printf("%08x ",alldeps[i].dep[j]); printf("\n");
	}
	for(int i=0; i<MAXDEPS; i++) stat.cdep.dep[i] = 0;
	tcomp.mark(true);
	cutcomp(&stat, alldeps, allres);

#ifdef TANDEM_DEBUG_TIMINGS
	tcomp.mark(false);
	printf("compute_cutsets(): t_comp_new=%.3lf  n_sets=%i\n", tcomp.time(), allres.size());
#endif

#ifdef TANDEM_CUTCOMP_OLDVERSION
	tcomp.mark(true);
	cutcomp_slow(res, num_deps, nodedeps, allres);
#ifdef TANDEM_DEBUG_TIMINGS
	tcomp.mark(false);
	printf("compute_cutsets(): t_comp_old=%.3lf  n_sets=%i\n", tcomp.time(), allres.size());
#endif
#endif

	// eliminate supersets in 'allres'
	for(std::set<std::set<uint> >::iterator it=allres.begin(); it!=allres.end(); it++)
	{
		CutNodelist cn;
		bool found = false;
		// find any subsets of this
		for(std::set<std::set<uint> >::iterator kt=allres.begin(); kt!=allres.end(); kt++)
		{
			if(it == kt) continue;
			if(includes(it->begin(), it->end(), kt->begin(), kt->end()))
			{
				found = true;
				break;
			}
		}
		if(found) continue;

		for(std::set<uint>::iterator jt=it->begin(); jt!=it->end(); jt++)
			cn.push_back(*jt);
		cn.push_back(num_nodes);
		result.push_back(cn);
	}

#ifdef TANDEM_DEBUG_CUTS
	printf("Tandem::compute_cutsets():  %i sets found:\n", result.size());
	for(uint i=0; i<result.size(); i++)
	{
		printf("    set %02i = { ", i);
		for(uint j=0; j<result[i].size(); j++)
		{
			printf("%i ", result[i][j]);
		}
		printf("}\n");
	}
#endif

	return result;
}


void Tandem::filter_cutsets(std::vector<CutNodelist> &cs, int deltalen)
{
	uint len_max = 0;
	uint len_min = 0xFFFFFFFF;
	for(int i=0; i<cs.size(); i++)
	{
		uint len = cs[i].size();
		if(len > len_max) len_max = len;
		if(len < len_min) len_min = len;
	}

	len_max = len_min + deltalen;

	for(int i=0; i<cs.size(); i++)
	{
		uint len = cs[i].size();
		if(len > len_max) { cs.erase(cs.begin()+i); i--; }
	}
}

void Tandem::filter_cutsets_v2(std::vector<CutNodelist> &cs, int deltalen)
{
	std::vector<CutNodelist> result;
	uint len_min = 0xFFFFFFFF;
	uint n = cs.size();
	for(uint i=0; i<n; i++)
	{
		uint len = cs[i].size();
		if(len < len_min) len_min = len;
	}

	uint len_max = len_min + deltalen;

	for(uint i=0; i<n; i++)
	{
		if(cs[i].size() <= len_max)
			result.push_back(cs[i]);
	}
	cs = result;
}


bool Tandem::cutcomp_slow(bool *res2, uint ndeps, std::set<uint> nodedeps[], std::set<std::set<uint> > &allres)
{
	if(res2==NULL) return false;
	byte cdeps[ndeps]; // current set of satisfied dependencies
	memset(cdeps,0,ndeps);
	uint currnode=0; // highest order node which has been checked so far

	// compute cdeps
	for(uint i=1; i<num_nodes; i++)
	{
		if(!res2[i]) continue;
		currnode = i;
		std::set<uint> &d = nodedeps[i]; // dependencies satisfied by node *it
		for(std::set<uint>::iterator jt=d.begin(); jt!=d.end(); jt++)
			cdeps[*jt] = 0x01;
	}

	//printf("res={ ");for(std::set<uint>::iterator it=res.begin(); it!=res.end(); it++)printf("%i ",*it);printf("}\n");
	//printf("cdeps={ ");for(int i=0; i<ndeps; i++)if(cdeps[i]) printf("%i ",i);printf("}\n");
	for(uint n=currnode+1; n<num_nodes; n++)
	{
		//printf("check node %i\n",n);
		// don't consider nodes already into result
		if(res2[n]) continue;
		// avoid the case of three cuts at consecutive nodes
		if((n>=2) && res2[n-1] && res2[n-2]) continue;
		// check if this node satisfies new dependencies
		std::set<uint> &dep = nodedeps[n];
		int found = -1;
		for(std::set<uint>::iterator it=dep.begin(); it!=dep.end(); it++)
		{
			if(!cdeps[*it])
			{
				found = *it;
				break;
			}
		}
		if(found<0) continue; // no new dependencies, prune
		// add node n
		//printf("add node %i (new dep=%i) --> ",n,found);
		res2[n]=true;
		byte ddeps[ndeps]; memcpy(ddeps, cdeps, ndeps);
		// update list of satisfied dependencies
		for(std::set<uint>::iterator jt=dep.begin(); jt!=dep.end(); jt++)
					ddeps[*jt] = 0x01;
		// check if all dependencies are now satisfied (termination condition)
		bool finished = true;
		for(uint i=0; i<ndeps; i++) if(!ddeps[i]) { finished=false; break; }
		//printf("ddeps={ ");for(int i=0; i<ndeps; i++)if(ddeps[i]) printf("%i ",i);printf("} - finished=%s\n",(finished?"true":"false"));

		if(!finished) finished = cutcomp_slow(res2, ndeps, nodedeps, allres);
		else {
			//printf("found set={ ");for(int h=0;h<num_nodes;h++)if(res2[h])printf("%i ",h);printf("}\n");
			std::set<uint> tmpset;
			for(uint i=0; i<num_nodes; i++) if(res2[i]) tmpset.insert(i);
			allres.insert(tmpset);
			//printf("set max_depth = %i\n", maxdepth);
		}

		res2[n]=false;
	}
	return false;
}


void Tandem::cutcomp(compstatus *stat, depsse deps[], std::set<std::set<uint> > &allres)
{
	compstatus stat2;
	bool found;
	uint64 tmp;
	if(stat == NULL) return;
	for(uint n=stat->cnode+1; n<num_nodes; n++)
	{
		//printf("checking node %i\n",n);
		if(stat->cres & (1L << n)) continue;
		if(n>=2) tmp = 3L << (n-2);
		if((n>=2) && ((stat->cres & tmp)==tmp)) continue;
		depsse &newdep = deps[n];
		uint depsok = 0;
		found = false;
		for(int i=0; i<MAXDEPS; i++)
		{
			stat2.cdep.dep[i] = stat->cdep.dep[i] | newdep.dep[i];
			//printf("dep%02i: old=%08x curr=%08x  res=%08x\n",i,stat->cdep.dep[i],newdep.dep[i],stat2.cdep.dep[i]);
			if(stat2.cdep.dep[i] == 0xFFFFFFFF) depsok++;
			if(stat2.cdep.dep[i] != stat->cdep.dep[i]) { found=true; }
		}
		if(!found) continue;
		//printf("adding node %i; depsok=%i/%i\n",n,depsok,MAXDEPS);
		stat2.cres = stat->cres | (1L << n);
		stat2.cnode = n;
		if(depsok<MAXDEPS)
		{
			cutcomp(&stat2, deps, allres);
		}
		else {
			std::set<uint> tmpset;
			for(uint i=0; i<num_nodes; i++)
			{
				if(stat2.cres & (1L << i)) tmpset.insert(i);
			}
			//printf("found set={ ");for(std::set<uint>::iterator jt=tmpset.begin();jt!=tmpset.end();jt++)printf("%i ",*jt);printf("}\n");
			allres.insert(tmpset);
		}
	}
}


/* Cut the tandem into a vector of sub-tandems according to the specified list of cuts
 * Nodes in 'cutset' must appear in ascending order (sorted vector)
 * Returns a vector of sub-tandems, flows are initialized with original parameters
 */
std::vector<Tandem> Tandem::cut(CutNodelist &cutset)
{
	std::vector<Tandem> result;
	int ncuts = cutset.size();
	if(ncuts == 0)
	{
		result.push_back(*this);
		return result;
	}
#ifdef TANDEM_DEBUG_CUTS
	printf("Tandem::cut(): set contains the following cuts: { ");
	for(int i=0; i<ncuts; i++)
		printf("%i ", cutset[i]);
	printf("}\n");
#endif
	// make the cuts
	uint cnode = 0;
	for(int i=0; i<ncuts; i++)
	{
		Tandem t;
		t.Clear(); // this call is redundant...
		// fill-in the nodes array: determine the subset of nodes
		uint cutnode = cutset[i]; // zero-based numbering
		if(cutnode > num_nodes) continue; // little safety check
		t.num_nodes = cutnode-cnode;
		t.nodes_arr = new Node [t.num_nodes];
#ifdef TANDEM_DEBUG_CUTS
		printf("Sub-tandem %i: nodes %i-%i (len=%i)\n",i,cnode,cutnode-1,t.num_nodes);
#endif
		// copy the subset of nodes
		for(uint n=0; n<t.num_nodes; n++)
		{
			t.nodes_arr[n] = nodes_arr[cnode+n];
		}
		// fill-in the flows array: split the flows affected
		t.flows_arr = new Flow [num_flows]; // we'll have at most 'num_flows' flows in the sub-tandem
		// the sub-tandem spans from node 'cnode' (included) to 'cutset[i]' (excluded)
		t.num_flows = 0;
		for(uint f=0; f<num_flows; f++)
		{
			// handle the tagged flow assignement
			if(f==tagged_flow) t.tagged_flow = t.num_flows;
			// if the flow does not cross the current sub-tandem, skip it
			if((flows_arr[f].exit < cnode) || (flows_arr[f].src >= cutnode))
			{
#ifdef TANDEM_DEBUG_CUTS
				printf("  skip (%i,%i)\n",flows_arr[f].src,flows_arr[f].exit);
#endif
				continue;
			}
			// this flow intersects the sub-tandem: trim it (i,j) --> (i,h) or (k,j)
			t.flows_arr[t.num_flows] = flows_arr[f];
			t.flows_arr[t.num_flows].src = _max(cnode, flows_arr[f].src) - cnode;
			t.flows_arr[t.num_flows].exit = _min(cutnode-1, flows_arr[f].exit) - cnode;
			bool cut_out = flows_arr[f].exit > cutnode-1;
			bool cut_in = flows_arr[f].src < cnode;
#ifdef TANDEM_DEBUG_CUTS
			printf("  cut  (%i,%i) as (%i,%i)  cut_in=%i cut_out=%i\n",flows_arr[f].src,flows_arr[f].exit,t.flows_arr[t.num_flows].src,t.flows_arr[t.num_flows].exit, cut_in, cut_out);
#endif
			t.num_flows++;
		}
		t.ComputeNestingLevels();
		// store the new sub-tandem
		result.push_back(t);
		// move the current node pointer after the cut
		cnode = cutnode;
	}
	return result;
}


uint Tandem::aggregate_flows(bool merge_tf)
{
	//printf("Tandem::aggregate_flows(%s): %i flows, tagged_flow=%i\n",(merge_tf?"true":"false"),num_flows,tagged_flow);
	for(uint i=0; i<num_flows-1; i++)
	{
		// check if this flow can be aggregated with others
		for(uint j=i+1; j<num_flows; j++)
		{
			if((flows_arr[i].src == flows_arr[j].src) && (flows_arr[i].exit == flows_arr[j].exit))
			{
				if(!merge_tf && ((j==tagged_flow) || (i==tagged_flow))) continue;
				// aggregate j into i
				flows_arr[i].burst += flows_arr[j].burst;
				flows_arr[i].rate += flows_arr[j].rate;
				if(j==tagged_flow) tagged_flow = i;
				else if(tagged_flow>j) tagged_flow--;
				// shift down the flows array to eliminate j
				for(uint k=j; k<num_flows-1; k++)
				{
					flows_arr[k] = flows_arr[k+1];
				}
				// write junk into the left space, for safety...
				flows_arr[i].uid = FLOW_INVALID_ID;
				flows_arr[num_flows-1].src = flows_arr[num_flows-1].exit = FLOW_INVALID_ID;
				num_flows--;
				j--; // compensate for shifted indexes
			}
		}
	}
	ComputeNestingLevels();
	// correct the nesting level of the tagged flow (when ar != 0)
	flows_arr[tagged_flow].nesting_level = 1;
	//printf("Tandem::aggregate_flows(): updated tagged_flow = %i, new num flows = %i\n",tagged_flow,num_flows);
	return tagged_flow;
}

/*	Compute the output arrival curve for flow 'fid' at the exit of the current tandem
 * 	and store the parameters (sigma, rho) into 'fdest'
 */
bool Tandem::compute_output_arrival_curve(uint fid, Flow *fdest, int mode)
{
	if((fid>=num_flows) || (fdest==NULL))
	{
		printf("Tandem::compute_output_arrival_curve(fid=%i, fdest=%x): invalid parameter(s), bug?\n", fid, fdest);
		return false;
	}
	Flow &f = flows_arr[fid];
	if(f.exit != num_nodes-1)
	{
		printf("Tandem::compute_output_arrival_curve(fid=%i): flow does not exit the tandem, bug?\n", fid);
		return false;
	}
	double Dmin = 0.0;
	// check for special case: single-flow tandem
	if(num_flows == 1)
	{
		// do the easy math
		Dmin = LUDB_trivial();
		fdest->rate = f.rate;
		fdest->burst = f.burst + f.rate * Dmin;
		return true;
	}
	// compute Dmin according to equation (13), section 4.1
	uint nchildren = nesting_tree->num_children();
	if(!nchildren)
	{
		printf("Tandem::compute_output_arrival_curve(fid=%i): error, root node has no direct children\n",fid);
		printf("Affected tandem is:\n"); Print();
		return false;
	}
	// compute the sum of LUDBs of the first-order child tnodes
	for(uint i=0; i<nchildren; i++)
	{
		uint64 nsimpl;
		Dmin += LUDB_experimental_quiet(nesting_tree->tnodes[i], mode, nsimpl);
	}
	// add the latency of C(1,N) of the possible leaves of the root node
	for(uint i=0; i<nesting_tree->num_leaves; i++)
	{
		uint n=nesting_tree->leaves[i];
		Dmin += nodes_arr[n].latency;
	}
	fdest->rate = f.rate;
	fdest->burst = f.burst + f.rate*Dmin;
	return true;
}


/* Compute all the necessary output arrival curves at the current tandem
 * Input: 't_next' is a reference to the subsequent tandem in the chain
 */
void Tandem::output_arrival_curves(Tandem &t_next, int mode)
{
	for(uint f=0; f<num_flows; f++)
	{
		Tandem tt = *this; // a copy of the original sub-tandem
		Flow *tf = tt.getFlow(f);
		if(tf == NULL) continue;
		// determine if the arrival curve of flow f is needed at the next tandem
		Flow *fnext = t_next.getFlowId(tf->uid);
		if(fnext == NULL) continue;
		// more checks
		if(tf->exit != tt.NumNodes()-1) continue; // flow exits midway, skip
		// this output curve is actually needed
		if(tf->src > 0)
		{
			// flow enters halfway, must split this sub-tandem recursively
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			printf("    !!!  flow=#%i uid=%i enters at node %i, must cut this sub-tandem there\n",f,tf->uid,tf->src);
#endif
			// perform the cut
			CutNodelist cs;
			cs.push_back(tf->src);
			cs.push_back(tt.NumNodes());
			std::vector<Tandem> subt = tt.cut(cs);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			printf("    splitted the subtandem into %i sub-subtandems (should be 2)\n",subt.size());
			for(uint i=0; i<subt.size(); i++) subt[i].Print();
#endif
			if(subt.size() != 2)
			{
				printf("    wrong size of sub-splitted tandem\n");
				continue;
			}
			// recursively compute the OAC's required for the 2nd sub-tandem
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			printf("    recursion enter\n");
#endif
			subt[0].output_arrival_curves(subt[1], mode);
			// now proceed with the computation of the original OAC (2nd sub-tandem only is involved)
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			printf("    recursion exit\n");
#endif
			uint tagged = subt[1].getFlowIdx(tf->uid); // remap the tagged flow in original tandem into 2nd sub-subtandem
			if(tagged == FLOW_INVALID_ID) {printf("BUG: unable to map tagged flow in sub-subtandem!\n"); continue;}
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			printf("    mapped tagged flow (uid=%i) in sub-subtandem found at fid=%i\n",tf->uid,tagged);
#endif
			subt[1].SetTaggedFlow(tagged);
			uint tf2=subt[1].aggregate_flows(false);
			subt[1].BuildNestingTree();
			bool ok = subt[1].compute_output_arrival_curve(tf2, fnext, mode);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
			if(ok) printf("    [OK] output arrival curve for tmp_tf=%i (uid=%i): burst=%.3lf, rate=%.3lf\n", tf2, tf->uid, fnext->burst, fnext->rate);
#endif
			continue;
		}
		tt.SetTaggedFlow(f);
		uint tmptf = tt.aggregate_flows(false);
		tt.BuildNestingTree();
		// compute the output arrival curve
		bool res;
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
		printf("    output_arrival_curve(): compute_output_arrival_curve(tmp_tf=%i uid=%i)\n", tmptf, tf->uid);
#endif
		res = tt.compute_output_arrival_curve(tmptf, fnext, mode);
#ifdef TANDEM_DEBUG_OUTPUT_ARR_CURVES
		if(!res) printf("   compute_output_arrival_curve(tmp_tf=%i) returned error\n",tmptf);
		else printf("    [OK] output arrival curve for tmp_tf=%i (uid=%i): burst=%.3lf, rate=%.3lf\n", tmptf, tf->uid, fnext->burst, fnext->rate);
#endif
	}
}


// Handle LUDB computation in trivial cases
double Tandem::LUDB_trivial()
{
	// case 1): one flow crossing the whole tandem
	if(num_flows == 1)
	{
		double latency = 0.0;
		double rate = TNODE_LUDB_INFINITY;
		for(uint i=0; i<num_nodes; i++)
		{
			latency += nodes_arr[i].latency;
			rate = _min(rate, nodes_arr[i].rate);
		}
#ifdef TANDEM_DEBUG_LUDB
		printf("Tandem::LUDB_trivial(): case 1) applied.\n");
#endif
		cached_ludb = latency + flows_arr[0].burst / rate;
		return cached_ludb;
	}
	// no supported cases detected, return error
	return -1000.0;
}

/* ************************************************************************************************************* */
/* ********************                   Lower-Bound related methods                       ******************** */
/* ************************************************************************************************************* */



// Auxiliary function which returns the n-th bit of the value x
bool inline Tandem::bit(uint64 x, uint n)
{
	uint64 mask = (uint64) pow(2,n);
	return (x & mask) ? true : false;
}

// Returns true is Flow f traverses Node n
bool inline Tandem::isNodeCrossedByFlow(uint n, uint f)
{
	if(f >= num_flows) return false;
	return ((flows_arr[f].src <= n) && (flows_arr[f].exit >= n));
}

// Returns true is Flow f enters the tandem at Node n
bool inline Tandem::isNodeEnteredByFlow(uint n, uint f)
{
	if(f >= num_flows) return false;
	return (flows_arr[f].src == n);
}


// Returns true is Flow f leaves the tandem at Node n
bool inline Tandem::isNodeLeftByFlow(uint n, uint f)
{
	if(f >= num_flows) return false;
	return (flows_arr[f].exit == n);
}


//	Check if the given flow interferes with the tagged flow at the specified node
bool inline Tandem::isInterferingFlow(uint n, uint f, uint tagged)
{
	uint n0 = flows_arr[tagged].src;
	if(!isNodeCrossedByFlow(n, f) || isNodeCrossedByFlow(n0, f)) return false;
	// more conditions
	if(flows_arr[f].exit < n0) return false;
	if(flows_arr[f].src > flows_arr[tagged].exit) return false;
	return true;
}


/*	This function performs FIFO de-multiplexing at the node exit point.
	The CDF_i contribution of the specified flow is devised from the overall CDF
	considering the CAF_i and CAF.
*/
bool Tandem::compute_CDF_i(uint node_id, uint flow_id)
{
	Curve caf, cdf, caf_i, cdf_i;
	double caf_latency=0.0, cdf_latency=0.0, caf_i_latency=0.0;
	int s_caf, s_caf_i, s_cdf, s_cdf_i;
	double x_caf, y_caf, x_caf_i, x_cdf, x_cdf_i;
	LinearSegment segcdf_i, segcaf, segcdf, segcaf_i;
	LinearSegment segcdf_next, segcaf_next;
	bool res;

	if(node_id >= num_nodes)
	{
#ifdef TANDEM_DEBUG_CDF_I
		printf("Tandem::compute_CDF_i() invoked with invalid node ID.\n");
#endif
		return false;
	}

	Node &node = nodes_arr[node_id];

	caf_i = node.CAF_i[flow_id];
	caf = node.CAF;
	cdf = node.CDF;

#ifdef TANDEM_DEBUG_CDF_I
	printf("Computing Node%i.CDF_%02i:\n",node_id,flow_id);
	printf("CAF-original = ");caf.Print();
	printf("CAF_%02i-original = ",flow_id);caf_i.Print();
	printf("CDF-original = ");cdf.Print();
#endif

	int nsegs_caf = caf.num_segments();
	int nsegs_cdf = cdf.num_segments();
	int nsegs_caf_i = caf_i.num_segments();
	int nsegs_cdf_i = 0;

	if((nsegs_caf == 0) || (nsegs_cdf == 0))
	{
		printf("Tandem::compute_CDF_i(): empty CAF or CDF at node %i\n",node_id);
		return false;
	}

	// CAF and CAF_i should have the same latency, so as CDF and CDF_i
	// Compared to them, CDF is retarded by the delay of the node's service curve (beta)
	cdf_latency = node.CDF.getLatency();
	caf_latency = node.CAF.getLatency();
	caf_i_latency = node.CAF_i[flow_id].getLatency();

	caf.removeLatency();
	if(caf_latency > 0.0)
	{
#ifdef TANDEM_DEBUG_CDF_I
		printf("shifting CAF_i by %lf horizontally\n",-caf_latency);
#endif
		caf_i.shiftRight(-caf_latency);
	}
	cdf.removeLatency();

	// smooth CAF (remove colinear segments and recompute rates)
	caf.remove_colinear_segments();

	// update segments counts
	nsegs_caf = caf.num_segments();
	nsegs_cdf = cdf.num_segments();
	nsegs_caf_i = caf_i.num_segments();

#ifdef TANDEM_DEBUG_CDF_I
	printf("Original latencies: CAF=%lf, CAF_i=%lf, CDF=%lf\n",caf_latency,caf_i_latency,cdf_latency);
	printf("Latencies removed, all curves now centered in origin\n");
	printf("CAF = ");caf.Print();
	printf("CAF_%02i = ",flow_id);caf_i.Print();
	printf("CDF = ");cdf.Print();
#endif


	// pre-process CAF to add any missing y-breakpoints which are present in CDF instead
	res = false;
	for(int i=0; i<nsegs_cdf; i++)
	{
		segcdf = cdf.getSeg(i);
		if(caf.has_BPY_at(segcdf.y) < 0)
		{
#ifdef TANDEM_DEBUG_CDF_I
			printf("CAF pre-processing: CDF has a BP at Y=%lf (X=%lf) which is not present in CAF, added\n",segcdf.y,segcdf.x);
#endif
			caf.add_BPY(segcdf.y);
			res = true;
		}
	}
	if(res)
	{
		nsegs_caf = caf.num_segments();
#ifdef TANDEM_DEBUG_CDF_I
		printf("CAF-updated = ");caf.Print();
#endif
	}

	// pre-process CAF to add any missing y-breakpoints which are present in CAF_i instead
	res = false;
	for(int i=0; i<nsegs_caf_i; i++)
	{
		segcaf_i = caf_i.getSeg(i);
		if(caf.has_BPX_at(segcaf_i.x) < 0)
		{
#ifdef TANDEM_DEBUG_CDF_I
			printf("CAF pre-processing: CAF_i has a BP at X=%lf (Y=%lf) which is not present in CAF, added\n",segcaf_i.x,segcaf_i.y);
#endif
			caf.add_BPX(segcaf_i.x);
			res = true;
		}
	}
	if(res)
	{
		nsegs_caf = caf.num_segments();
#ifdef TANDEM_DEBUG_CDF_I
		printf("CAF-updated = ");caf.Print();
#endif
	}

	cdf_i.Zero();

	s_caf = 0;
	s_cdf = 0;
	x_caf = 0.0; y_caf = 0.0; x_cdf = 0.0; x_cdf_i = 0.0;
	double txbits, mytxbits;

	segcdf_i.Zero();
	while(s_caf < nsegs_caf)
	{
		// fetch current and next CAF segments
		segcaf = caf.getSeg(s_caf);
		if(s_caf < nsegs_caf-1) segcaf_next = caf.getSeg(s_caf+1);
		else segcaf_next = segcaf;
		// determine the number of bits transmitted in the current segment
		txbits = segcaf_next.y - segcaf.y;
		y_caf = segcaf.y + txbits / 2.0; // mid-point y in CAF segment of interest
		x_cdf = cdf.f_inv(y_caf, false);
		s_cdf = cdf.getSegmentDefining(x_cdf);
		segcdf = cdf.getSeg(s_cdf);

#ifdef TANDEM_DEBUG_CDF_I
		printf("CAF segment: (%lf, %lf) <--> (%lf, %lf)   total bits transmitted = %lf\n",segcaf.x,segcaf.y,segcaf_next.x,segcaf_next.y,txbits);
#endif

		if((IS_ZERO(segcaf.rate)) || (txbits < CURVE_EPSILON))
		{
			s_caf++;
			continue;
		}

		if(segcaf.id >= 0)
		{
			// someone's burst is being transmitted
			if(segcaf.id == flow_id)
			{
				// this flow is bursting, assign all CDF's rate to it
				segcdf_i.rate = segcdf.rate;
				mytxbits = txbits;
#ifdef TANDEM_DEBUG_CDF_I
				printf("my burst: %.3lf bits transmitted by me\n",mytxbits);
#endif
			}
			else {
				// a different flow is bursting, this flow's rate is zero
				segcdf_i.rate = 0.0;
				mytxbits = 0.0;
#ifdef TANDEM_DEBUG_CDF_I
				printf("alien burst (flow %02i), %lf bits transmitted by me\n",segcaf.id,mytxbits);
#endif
			}
			x_caf = segcaf.x;
			x_cdf = cdf.f_inv(segcaf_next.y, false);
		}
		else {
			// no burst is being transmitted
			x_caf = (segcaf_next.x + segcaf.x) / 2.0; // mid-point x in CAF interval of interest
			double caf_i_slope = caf_i.getSlope(x_caf, false);
			//x_cdf = cdf.f_inv(y_caf, false); // left or right limit???
			double cdf_slope = cdf.getSlope(x_cdf, false);
			if(segcaf.rate == 0.0)
			{
				// nobody was transmitting in the considered interval
				segcdf_i.rate = 0.0;
				mytxbits = 0.0;
				printf("Tandem::compute_CDF_i(): should never get here, inconsistency detected!\n");
			}
			else {
				// compute the rate quota of this flow
				x_cdf = cdf.f_inv2(segcaf_next.y,false);
				segcdf_i.rate = cdf_slope * (caf_i_slope / segcaf.rate);
				mytxbits = segcdf_i.rate * (x_cdf - segcdf_i.x);
			}
#ifdef TANDEM_DEBUG_CDF_I
			printf("Midpoint for CDF: Y=%lf   Midpoint for CAF/CAF_i: X=%lf\n",y_caf,x_caf);
			printf("affected CDF interval: %lf <--> %lf   slope = %lf ",cdf.f_inv(segcaf.y,false),x_cdf,cdf_slope);
			printf("  (segment: "); segcdf.Print(); printf(")\n");
			printf("affected CAF_i interval: %lf <--> %lf   slope = %lf   my bits transmitted = %lf\n",segcaf.x,segcaf_next.x,caf_i_slope,mytxbits);
#endif
		}
		//printf("cdf_i.add_segment(): ");segcdf_i.Print();printf("\n");
		cdf_i.add_segment(segcdf_i);
		segcdf_i.y += mytxbits;
		segcdf_i.x = x_cdf;
		s_caf++;
	}
	segcdf_i.rate = 0.0;
	//printf("add-last-segment: ");segcdf_i.Print();printf("\n");
	cdf_i.add_segment(segcdf_i);

#ifdef TANDEM_DEBUG_CDF_I
	printf("CDF_%02i = ",flow_id); cdf_i.Print();
#endif

	// sanity checks
	res = true;
#ifdef TANDEM_DEBUG_CHECKS
	if(fabs(cdf_i.getTotalTraffic() - caf_i.getTotalTraffic()) > LINEARSEGMENT_EPSILON)
	{
		printf("Tandem::compute_CDF_%02i(): sanity check error: bits count mismatch at Node%i between CAF_%02i (%.3lf) and CDF_%02i (%.3lf), bug somewhere!\n",flow_id,node_id,flow_id,caf_i.getTotalTraffic(),flow_id,cdf_i.getTotalTraffic());
		printf("difference: %lf\n",caf_i.getTotalTraffic() - cdf_i.getTotalTraffic());
		res = false;
	}
#endif

	// finally re-add the CDF's latency to the CDF_i
	cdf_i.shiftRight(cdf_latency);
	// cleanup the resulting curve by removing colinear segments and recomputing rates (smoothing)
	cdf_i.remove_colinear_segments();

	// post-process the curve segments in order to reduce nasty rounding errors, etc.
#ifdef TANDEM_DEBUG_CHECKS
	cdf_i.SanityCheck();
#endif

#ifdef TANDEM_DEBUG_CDF_I
	printf("Tandem::compute_CDF_%02i() finished\n",flow_id);
#endif
	node.CDF_i[flow_id] = cdf_i;
	return res;

}


/*	This method processes all nodes in the tandem (starting from the node where the tagged flow enters)
	and computes the delay experienced by the last bit of the tagged flow at the exit of the tandem
	Each interfering flow which enters the path along the way can transmit a burst at the beginning or at the end of the interval;
	the particular choice made for this run is encoded by the i-th bit of the input parameter
	Returns a negative number in case of errors.
*/
double Tandem::compute_delay(bool *combo)
{
	double delay = 0.0;
	double x1, x2;
	double tcaf,tcdf;
	bool res;

#ifdef TANDEM_DEBUG_LOWERBOUND
	printf("Tandem::compute_delay(): combo ");
	for(int i=0; i<num_flows; i++) printf("%i",(int)combo[i]);
	printf("\n");
#endif
	if(combo == NULL) return TANDEM_LB_ERROR_INVALID_PARAM;

	// let's start from the node where the tagged flow enters
	uint current_node = flows_arr[tagged_flow].src;
	uint last_node = flows_arr[tagged_flow].exit;
	// the active interval for incoming flows is initially null
	x1=0.0;
	x2=0.0;

	while(current_node <= last_node)
	{
#ifdef TANDEM_DEBUG_LOWERBOUND
		printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Tandem::compute_delay(): processing Node %i <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n",current_node);
#endif
		Node &node = nodes_arr[current_node];

		// first, we must determine the interval (x1,x2) for the CAF_i of new flows entering the tandem at this node
		if(current_node > flows_arr[tagged_flow].src)
		{
			double xx1 = CURVE_POS_INFINITY;
			double xx2 = 0.0;
			for(int i=0; i<num_flows; i++)
			{
				if(isNodeCrossedByFlow(current_node,i) && !isNodeEnteredByFlow(current_node,i))
				{
					double x = nodes_arr[current_node-1].CDF_i[i].getLastBit(tagged_flow);
					if(x > xx2) xx2 = x;
					x = nodes_arr[current_node-1].CDF_i[i].getFirstBit(tagged_flow);
					if(x < xx1) xx1 = x;
				}
			}
			x1 = nodes_arr[current_node-1].CDF.getFirstBit(tagged_flow);
			x2 = nodes_arr[current_node-1].CDF.getLastBit(tagged_flow);
			x1 = xx1; x2 = xx2;
#ifdef TANDEM_DEBUG_LOWERBOUND
			printf("Node%i.CAF is active between x1=%lf and x2=%lf\n",current_node,xx1,xx2);
			printf("New flows will transmit between x1=%lf and x2=%lf\n",x1,x2);
#endif
		}

		// let's setup CAF and the various CAF_i in input at this node
		node.CAF.Zero();
		//Curve caftmp;
		for(int i=0; i<num_flows; i++)
		{
			bool entering = isNodeEnteredByFlow(current_node, i);
			bool crossing = isNodeCrossedByFlow(current_node, i);

			// check if this flow is entering the tandem at this node
			if(entering)
			{
				// create its brand-new CAF_i
				flows_arr[i].create_caf(x1, x2, combo[i], i);
				//TODO: maybe we should check if this is exactly the tagged flow and do something special? (probably not)
				node.CAF_i[i] = flows_arr[i].caf;
#ifdef TANDEM_DEBUG_LOWERBOUND
				printf("Node%i.CAF%i (new flow) = ",current_node,i);node.CAF_i[i].Print();
#endif
			}
			// if not, check whether this flow has entered the tandem at an upstream node
			else if(crossing)
			{
				// copy its CAF_i from the CDF_i at the previous node
				if(current_node < 1) {printf("compute_delay(): SANITY CHECK 1 FAILED - CRASH!!!\n"); break;}
				node.CAF_i[i] = nodes_arr[current_node-1].CDF_i[i];
#ifdef TANDEM_DEBUG_LOWERBOUND
				printf("Node%i.CAF%i = Node%i.CDF%i = ",current_node,i,current_node-1,i);node.CAF_i[i].Print();
#endif
			}
			// this flow is not traversing the node at all, zero() the corresponding CAF_i entry
			else node.CAF_i[i].Zero();

			// sum the CAF_i just computed to the node's global CAF
			res = node.CAF.add(node.CAF_i[i]);
		}
#ifdef TANDEM_DEBUG_LOWERBOUND
		printf("Node%i.CAF = ",current_node);node.CAF.Print();
		printf("Node%i.Beta = ",current_node);node.beta.Print();
#endif
		// now we can compute CDF = CAF (x) beta
		node.CDF = node.CAF;
		res = node.CDF.ConvolveWithLatencyRate(node.beta);

		// perform some sanity checks about the computed CDF
#ifdef TANDEM_DEBUG_CHECKS
		tcaf = node.CAF.getTotalTraffic();
		tcdf = node.CDF.getTotalTraffic();
		if(!IS_EQUAL(tcaf, tcdf))
		{
			printf("Tandem::compute_delay(): bitcount mismatch at Node%i between CAF (%.3lf) and CDF (%.3lf) - bug in convolution?\n",current_node,tcaf,tcdf);
			printf("Computed CDF is:\nNode%i.CDF = ",current_node); node.CDF.Print();
			return TANDEM_LB_ERROR_INTERNAL;
		}
#endif

#ifdef TANDEM_DEBUG_LOWERBOUND
		printf("Node%i.CDF = ",current_node);node.CDF.Print();
#endif

		// and finally we compute the various CDF_i... note that this step is not necessary at the last node!
		for(int i=0; i<num_flows; i++)
		{
			bool crossing = isNodeCrossedByFlow(current_node, i);
			bool leaving = isNodeLeftByFlow(current_node, i);

			// check if the flow has just crossed the node and it's not leaving the tandem here
			if(crossing && !leaving)
			{
				// we must compute its CDF_i
				res = compute_CDF_i(current_node, i);
				//res = compute_CDF_i_new(current_node, i);
#ifdef TANDEM_DEBUG_LOWERBOUND
				printf("Node%i.CDF_%i = ",current_node,i); node.CDF_i[i].Print();
#endif
			}
			// this flow is leaving or it's not traversing the node at all: we don't need to compute the CDF_i
			else node.CDF_i[i].Zero();
		}

		// and that's it, we have finished to process this node
		current_node++;
	}

	delay = nodes_arr[last_node].CDF.getLastBit(tagged_flow);

	return delay;
}


/*	This is the main function which performs the LowerBound algorithm on the tandem
*/
double Tandem::LowerBound(bool quiet, uint64 start_combo, uint64 combo_range_len, float percentage)
{
	uint64 num_combos, eff_num_combos;
	double lowerbound = 0.0;
	clock_t t_start, t_end;
	double t_elapsed;
	RNG rng(1400);
	int max_flows = 8 * sizeof(uint64);
	bool b_combo[num_flows];
	std::vector<uint64> lb_combos;

	if(!quiet)
	{
		printf("Tandem::LowerBound() algorithm started.\n");
		printf("\n----------------------------------------------------\n");
	}

	//compute_output_rates(tagged_flow);

	if(num_flows > max_flows)
	{
		// replace uint with an IndexSet to overcome the limitation
		printf("ERROR: the current implementation does not support more than %i flows.\n", max_flows);
		return -1.0;
	}

	// clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
	for(uint i=0; i<num_nodes; i++) nodes_arr[i].clear_curves(num_flows);

	num_combos = (uint64) (pow(2,num_flows)-1); // we actually count <= 2^N - 1
	if(!quiet) printf("Number of combinations: %Lu  (2 ^ %i)\n",num_combos+1,num_flows);

	// to save time, we will only try "greedy" for flows which enter the tandem at node 0
	// in fact, "greedy" or "late" is the same at node 0 because flows only transmit the burst
	uint64 combo_mask = 0ULL;
	for(int i=0; i<num_flows; i++) if(isNodeEnteredByFlow(0, i)) combo_mask |= (1ULL << i);
	if(!quiet) printf("These flows transmit only the burst and will be considered greedy only:  ");
	eff_num_combos = 1ULL;
	for(int i=0; i<num_flows; i++)
	{
		if(isNodeEnteredByFlow(0,i))
		{
			if(!quiet) printf("%02i ",i);
			eff_num_combos <<= 1;
		}
	}
	eff_num_combos = (num_combos+1) / eff_num_combos;

	//start_combo = 1536;  start_combo = 98344;
	//excel2 #32000 problemi cdf_i al nodo 6
	//nest8-20 #156 problemi convoluzione nodo 6
	if(combo_range_len != 0) num_combos = start_combo + combo_range_len - 1;

	uint64 combo_base = start_combo;
	double combo_fstep = 100.0 / percentage;
	double combo_accum = 0.0, combo_istep = 0.0;
	uint64 combo_step = (uint64) (combo_fstep < 1.01) ? 1+floor(combo_fstep) : 1;
	uint64 combo = combo_base;

	if(!quiet)
	{
		printf("\nRandom combinations percentage: ");
		if(percentage >= 0.00001) printf("%.3f %%\n",percentage);
		else printf("%f %%\n",percentage);
		if(combo_fstep < 1.01) printf("Effective number of combinations that will be computed: %Lu\n",eff_num_combos);
		else {
			eff_num_combos = (uint64) floor(num_combos * (percentage / 100.0));
			printf("Approximate number of combinations that will be computed: %Lu\n", eff_num_combos);
		}
	}
	eff_num_combos = 0;

	// start the benchmark timer
	t_start = clock();
	// start the lower-bound computations!
	while(combo_base <= num_combos)
	{
#ifdef TANDEM_LB_RANDOM_COMBOS
		if(combo_fstep > 1.01)
		{
			combo_accum += combo_fstep;
			double r = modf(combo_accum, &combo_istep);
			if(combo_istep < 1.0) combo_step = 1;
			else {combo_step = (uint64)combo_istep; combo_accum = r;}
			combo = combo_base + (uint64) floor(combo_step * rng.uniform(0.0, 1.0));
			if(combo > num_combos) continue;
		}
		else {
			// optimization: compute one combination only involving the flows which enter the tandem at node 0
			if(combo_base & combo_mask)
			{
				combo_base++;
				continue;
			}
			combo = combo_base;
		}
#else
		if(combo_base & combo_mask)
		{
			combo_base++;
			continue;
		}
		combo = combo_base;
#endif
		// call the function which computes the delay for this combination and update the maximum
		double delay;
		for(int i=0; i<num_flows; i++) b_combo[i] = bit(combo,i);
#ifdef TANDEM_CATCH_EXCEPTIONS
		try {
			delay = compute_delay(b_combo);
		} catch (std::exception &ex) {
			printf("Tandem::LowerBound(combo=%Lu): exception caught inside compute_delay(): %s\n",combo,ex.what());
			delay = -1.0;
		}
#else
		delay = compute_delay(b_combo);
#endif
		if(delay < 0.0)
		{
			printf("Tandem::LowerBound(combo=%Lu): computation failed.\n",combo);
			combo_base += combo_step;
			continue;
		}
		if(delay > lowerbound + LINEARSEGMENT_EPSILON)
		{
			lb_combos.clear();
			lb_combos.push_back(combo);
			lowerbound = delay;
			if(!quiet) printf("#%08Lu: new LowerBound = %lf\n",combo,lowerbound);
		}
		else if(fabs(delay - lowerbound) < LINEARSEGMENT_EPSILON)
			lb_combos.push_back(combo);
		//lowerbound = _MAX(lowerbound, delay);
		//printf("Combo %06i: delay computed = %lf   curr_lowerbound = %.3lf\n",combo,delay,lowerbound);
		combo_base += combo_step;
		eff_num_combos++;
	}
	t_end = clock();
	t_elapsed = ((double)t_end - (double)t_start) / (double)CLOCKS_PER_SEC;

	if(!quiet)
	{
		printf("Lower Delay Bound = %lf obtained with the following combinations:\n", lowerbound);
		for(int i=0; i<lb_combos.size(); i++)
		{
			printf("  %03i: #%08Lu = ",i,lb_combos[i]);
			for(int j=0; j<num_flows; j++) printf("%i",(int)bit(lb_combos[i],j));
			printf("\n");
		}
		printf("Effective combinations computed: %Lu\n",eff_num_combos);
		printf("Computation time: %.3lf seconds", t_elapsed);
		if(t_elapsed > 0.1) printf("  (%.2lf combos/s)",(double)(eff_num_combos)/t_elapsed);
		printf("\n");
	}

	return lowerbound;
}


/*	This function compute the max output rates of flows at each node
*/
bool Tandem::compute_output_rates(uint tagged)
{
	int node_start = flows_arr[tagged].src;
	int node_end = flows_arr[tagged].exit;
	double rates[num_flows];
	double R;
	bool rate_reduction = false;

	if(tagged >= num_flows) return false;

	// handle the first node(s) separately
	for(int i=0; i<=node_start; i++)
		nodes_arr[i].set_output_rates(num_flows, NULL, false);
	nodes_arr[node_start].set_output_rate(tagged, nodes_arr[node_start].rate);

	// handle remaining nodes in the path of the tagged flow
	for(int n=node_start+1; n<=node_end; n++)
	{
#ifdef TANDEM_DEBUG_LBOPT
		printf("Node %02i:",n+1);
#endif
		R = 0.0;
		for(int f=0; f<num_flows; f++)
		{
#ifdef TANDEM_DEBUG_LBOPT
			printf(" "); flows_arr[f].Print();
#endif
			if(f == tagged)
			{
				rates[f] = nodes_arr[n-1].get_output_rate(f);
				R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
				printf("/tf=%.3lf ",rates[f]);
#endif
				continue;
			}
			if(!isInterferingFlow(n, f, tagged))
			{
				rates[f] = -1.0;
#ifdef TANDEM_DEBUG_LBOPT
				printf("/NI ");
#endif
				continue;
			}
			if(isInterferingFlow(n-1, f, tagged))
			{
				rates[f] = nodes_arr[n-1].get_output_rate(f);
				R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
				printf("/r%02i=%.3lf ",n,rates[f]);
#endif
			}
			else if(isNodeEnteredByFlow(n, f))
			{
				rates[f] = flows_arr[f].rate;
				R += rates[f];
#ifdef TANDEM_DEBUG_LBOPT
				printf("/rho=%.3lf ",rates[f]);
#endif
			}
		}
#ifdef TANDEM_DEBUG_LBOPT
		printf(" - R=%.3lf",R);
#endif
		if(R > nodes_arr[n].rate + LINEARSEGMENT_EPSILON)
		{
			// downscale rates
			double scale = nodes_arr[n].rate / R;
			for(int f=0; f<num_flows; f++)
			{
				if(rates[f] >= 0.0) rates[f] *= scale;
			}
			rate_reduction = true;
#ifdef TANDEM_DEBUG_LBOPT
			printf(", scale by %lf",scale);
#endif
		}
		else rate_reduction = false;
#ifdef TANDEM_DEBUG_LBOPT
		printf("\n");
#endif
		nodes_arr[n].set_output_rates(num_flows, rates, rate_reduction);
	}

	return true;
}


#define RULE_BURST_EARLY	0x00
#define RULE_BURST_LATE		0x01
#define RULE_BURST_BOTH		0x02
#define RULE_ORDER_NONE		0x00
#define RULE_ORDER_BEFORE	0x01
#define RULE_ORDER_AFTER	0x02
#define RULE_ORDER_CLASH	0x03

uint64 Tandem::compute_lowerbound_rules(uint tagged, byte *burst_time)
{
	int node;
	byte order_rules[num_flows][num_flows];
	memset(burst_time, RULE_BURST_BOTH, num_flows);
	memset(order_rules,RULE_ORDER_NONE, num_flows*num_flows);

	for(int f=0; f<num_flows; f++)
	{
#ifdef TANDEM_DEBUG_LBOPT
		printf("Tagged flow = "); flows_arr[f].Print(); printf("\n");
#endif
		compute_output_rates(f);
#ifdef TANDEM_DEBUG_LBOPT
		printf("\n");
#endif
	}

	for(int f=0; f<num_flows; f++)
	{
		compute_output_rates(f);
		for(int n = flows_arr[f].src; n<=flows_arr[f].exit; n++)
		{
			if(nodes_arr[n].is_rate_reducing())
			{
				// all interfering flows which leave the tandem after node <n> must trasmit their burst before
				// the first bit of flow <f> enters the node
				for(int l=0; l<num_flows; l++)
				{
					if((l == f) || (flows_arr[l].exit >= n)) continue;
					for(int m=flows_arr[l].src; m<=flows_arr[l].exit; m++)
					{
						if(isInterferingFlow(m,l,f))
						{
							// <l> transmits before <f>
							order_rules[l][f] |= RULE_ORDER_BEFORE;
						}
					}
				}
			}
			else {
				// all interfering flows which enter the tandem at node <n> must transmit their burst at the end,
				// just before the last bit of flow <f> enters the node
				for(int l=0; l<num_flows; l++)
				{
					if(l == f) continue;
					if(isNodeEnteredByFlow(n,l) && isInterferingFlow(n,l,f))
					{
						// <l> transmits at the end of <f>
						order_rules[l][f] |= RULE_ORDER_AFTER;
					}
				}
			}
		}
	}
	// now we have filled in the order_rules matrix
#ifdef TANDEM_DEBUG_LBOPT
	printf("Dump of rules matrix:\n");
	printf("   ");
	for(int i=0; i<num_flows; i++) printf("%02i|",i);
	printf("\n");
	for(int i=0; i<num_flows; i++)
	{
		printf("%02i|",i);
		for(int j=0; j<num_flows; j++)
		{
			switch(order_rules[i][j])
			{
				case	RULE_ORDER_NONE:
					printf("-");
					break;
				case	RULE_ORDER_BEFORE:
					printf("<");
					break;
				case	RULE_ORDER_AFTER:
					printf(">");
					break;
				case	RULE_ORDER_CLASH:
					printf("X");
					break;
				default: break;
			}
			printf(" |");
		}
		printf("\n");
	}
#endif
	// now determine the actual combinations
	uint64 num_combos = 1;
	for(int i=0; i<num_flows; i++)
	{
		if(order_rules[i][tagged] == RULE_ORDER_BEFORE)
			burst_time[i] = RULE_BURST_EARLY;
		else if(order_rules[i][tagged] == RULE_ORDER_CLASH)
		{
			burst_time[i] = RULE_BURST_BOTH;
			num_combos *= 2;
		}
		else burst_time[i] = RULE_BURST_LATE;
	}
	// extra reduction of combinations for flows which enter the tandem at node zero
	for(int i=0; i<num_flows; i++)
	{
		if(burst_time[i] != RULE_BURST_BOTH) continue;
		if(isNodeEnteredByFlow(flows_arr[tagged].src, i))
		{
			burst_time[i] = RULE_BURST_EARLY;
			num_combos /= 2;
		}
	}

	return num_combos;
}


/*	This method applies the exclusion criteria to determine
*/
double Tandem::LowerBound_experimental(float percentage)
{
	byte burst_time[num_flows];
	bool b_combo[num_flows];
	double lower_bound = 0.0;
	double delay = 0.0;
	uint64 num_combos;
	uint bitpos;
	RNG rng(1400);

	printf("Tandem::LowerBound_experimental() algorithm started.\n");
	printf("\n----------------------------------------------------\n");

	// clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
	for(int i=0; i<num_nodes; i++) nodes_arr[i].clear_curves(num_flows);

	num_combos = compute_lowerbound_rules(tagged_flow, burst_time);
	for(int i=0; i<num_flows; i++) b_combo[i] = (burst_time[i] == RULE_BURST_LATE);

	printf("\nEffective number of combinations: %Lu, total %Lu (2 ^ %i)\n",num_combos+1,(uint64)(pow(2.0,num_flows)),num_flows);

	for(uint64 combo=0; combo<num_combos; combo++)
	{
#ifdef TANDEM_LB_RANDOM_COMBOS
		if(percentage < 99.999)
		{
			double p = rng.uniform(0.0, 100.0);
			if(p > percentage) continue;
		}
#endif

		// set b_combo
		bitpos = 0;
#ifdef TANDEM_DEBUG_LBOPT
		printf("combo %06Lu: bitfield ",combo);
#endif
		for(int i=0; i<num_flows; i++)
		{
			if(burst_time[i] == RULE_BURST_BOTH)
			{
				b_combo[i] = bit(combo,bitpos);
				bitpos++;
#ifdef TANDEM_DEBUG_LBOPT
				printf(" %i ",(int)b_combo[i]);
#endif
			}
#ifdef TANDEM_DEBUG_LBOPT
			else printf("(%i)",(int)b_combo[i]);
#endif
		}

		// compute lower bound for this combination
		delay = compute_delay(b_combo);
		if(delay < 0.0)
		{
			printf("Combination %Lu: LowerBound computation failed.\n",combo);
			continue;
		}
		if(delay > lower_bound + LINEARSEGMENT_EPSILON)
		{
			lower_bound = delay;
		}
#ifdef TANDEM_DEBUG_LBOPT
		printf("  -->  delay = %.3lf,  LowerBound = %lf\n",delay,lower_bound);
#endif
	}

	return lower_bound;
}



double Tandem::LowerBound_experimental_quiet(uint64 &num_combos, float percentage)
{
	byte burst_time[num_flows];
	bool b_combo[num_flows];
	double lower_bound = 0.0;
	double delay = 0.0;
	uint bitpos;
	RNG rng(1400);

	// clear the curves at each node: beta, CAF, CAF_i, CDF and CDF_i
	for(int i=0; i<num_nodes; i++) nodes_arr[i].clear_curves(num_flows);

	num_combos = compute_lowerbound_rules(tagged_flow, burst_time);
	for(int i=0; i<num_flows; i++) b_combo[i] = (burst_time[i] == RULE_BURST_LATE);

	for(uint64 combo=0; combo<num_combos; combo++)
	{
#ifdef TANDEM_LB_RANDOM_COMBOS
		if(percentage < 99.999)
		{
			double p = rng.uniform(0.0, 100.0);
			if(p > percentage) continue;
		}
#endif

		// set b_combo
		bitpos = 0;
		for(int i=0; i<num_flows; i++)
		{
			if(burst_time[i] == RULE_BURST_BOTH)
			{
				b_combo[i] = bit(combo,bitpos);
				bitpos++;
			}
		}

		// compute lower bound for this combination
		delay = compute_delay(b_combo);
		if(delay < 0.0)
		{
			printf("Combination %Lu: LowerBound computation failed.\n",combo);
			continue;
		}
		if(delay > lower_bound + LINEARSEGMENT_EPSILON)
		{
			lower_bound = delay;
		}
	}

	return lower_bound;
}



/*	*** TEST code for LowerBound-related classes and operations ***
*/
void Tandem::Test_LowerBound_code()
{
	Curve beta, caf1, caf2, cdf;
	LinearSegment ls;
	bool res;

	Clear();
	flows_arr = new Flow[2];
	flows_arr[0].rate = 40.0;
	flows_arr[0].burst = 4.0;
	flows_arr[1].rate = 60.0;
	flows_arr[1].burst = 6.0;

	printf("\n--------- TWO FINAL BURSTS -------\n");

	flows_arr[0].create_caf(10.0, 11.0, true,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(10.0, 11.0, true,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();


	printf("\n--------- TWO INITIAL BURSTS -------\n");

	flows_arr[0].create_caf(10.0, 11.0, false,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(10.0, 11.0, false,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();

	printf("\n--------- TWO AT BEGINNING, TWO AT END -------\n");

	flows_arr[0].create_caf(10.0, 11.0, true,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(10.0, 11.0, false,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	caf1 = caf2;
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();

	printf("\n--------- ONE AT BEGINNING, TAGGED FLOW -------\n");

	flows_arr[0].create_caf(0.0, 0.0, true,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(0.0, 0.0, true,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();

	printf("\n--------- ONE CONTINUOS, ONE GAP -------\n");

	flows_arr[0].burst = 0.0;
	flows_arr[0].create_caf(10.0, 12.0, true,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(10.0, 11.0, true,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();

	printf("\n--------- TWO CONTINUOS -------\n");

	flows_arr[0].burst = 0.0; flows_arr[1].burst = 0.0;
	flows_arr[0].create_caf(10.0, 12.0, true,0); caf1 = flows_arr[0].caf;
	printf("CAF1 "); caf1.Print();
	caf2.create_token_bucket(0.0, 10.0); // sigma, rho
	flows_arr[1].create_caf(11.0, 13.0, true,1); caf2 = flows_arr[1].caf;
	printf("CAF2 "); caf2.Print();
	caf2.add(caf1);
	printf("CAF total "); caf2.Print();
	beta.create_latency_rate(10.0, 10.0); // theta, R
	printf("beta "); beta.Print();
	cdf = caf2;
	res = cdf.ConvolveWithLatencyRate(beta);
	if(!res) printf("Convolution failed\n");
	printf("CDF "); cdf.Print();

	Clear();
}

}
