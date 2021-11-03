/***************************************************************************
 *   Copyright (C) 2008 by Luca Bisti   								   *
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
#ifndef DEBORAHTANDEM_H
#define DEBORAHTANDEM_H

#include "types.h"
#include "node.h"
#include "flow.h"
#include "tnode.h"
#include "rng.h"
#include "deborah.h"
#include <set>


namespace deborah {

// This is the return value used to signal an error in LUDB-related functions
// (which deal with non-negative results only)
#define TANDEM_LUDB_ERROR_NESTINGTREE	-10000.0
#define TANDEM_LB_ERROR_INVALID_PARAM	-10000.0
#define TANDEM_LB_ERROR_INTERNAL	-20000.0

// When the following is defined, C++ exceptions handling is enabled in certain parts of the code (small performance penalty)
#define TANDEM_CATCH_EXCEPTIONS

// When the following is defined, the "--lb-random-combos" switch will be obeyed
#define TANDEM_LB_RANDOM_COMBOS

#define TANDEM_LUDB_MAXCUTS	0x7FFFFFFF

typedef struct {
	int order;
	int depth;
	double node_latency_min;
	double node_latency_max;
	double node_rate_overprov_min;
	double node_rate_overprov_max;
	double flow_rate_min;
	double flow_rate_max;
	double flow_burst_min;
	double flow_burst_max;
} balanced_tree_descriptor;


typedef std::vector<uint> CutNodelist;



// ----- Class declaration -----

class Tandem{
public:
    Tandem();
    Tandem(const Tandem &t2);
    Tandem(uint nodes, uint flows);

    ~Tandem();

	bool CreateTree(const balanced_tree_descriptor &bt_desc, RNG &rng);
	bool CreateFullTandem(const balanced_tree_descriptor &bt_desc, RNG &rng);
	bool Load(char *config_filename);
	bool Save(char *config_filename);
	void Clear();
	void Print();

	bool BuildNestingTree();
	void PrintNestingTree();
	bool isSinkTree();
	bool isNested();
	bool checkProvisioning();
	bool analyzeProvisioning(double &min, int &nodemin, double &max, int &nodemax);
	bool scaleProvisioning(double flow_mult, double node_mult);
	void PrintAllCombinations();
	void PrintSolution(Poly &p, uint dead_var);
	uint64 LUDB_TotalCombos();

	Node *getNode(uint index);
	Flow *getFlow(uint index);
	Flow *getFlowId(uint uid);
	uint getFlowIdx(uint uid);

	uint NumNodes();
	uint NumFlows();

	bool SetTaggedFlow(uint tf);
	uint GetTaggedFlow();

	// Non-nested tandem splitting methods
	std::vector<CutNodelist> compute_cutsets();
	void filter_cutsets(std::vector<CutNodelist> &cs, int deltalen);
	void filter_cutsets_v2(std::vector<CutNodelist> &cs, int deltalen);
	std::vector<Tandem> cut(CutNodelist &cutset);
	int join(std::vector<Tandem> &cuts);
	uint aggregate_flows(bool merge_tf);
	bool compute_output_arrival_curve(uint fid, Flow *fdest, int mode);
	void output_arrival_curves(Tandem &t_next, int mode);

	void Test_LUDB_code();
	void Test_LowerBound_code();

	// Old & dumb LUDB computation, which tries ALL simplexes
	double LUDB();

	// Optimized LUDB computation (verbose & quiet versions provided)
	double LUDB_experimental(int max_good_combos);
	double LUDB_experimental_quiet(int max_good_combos, uint64 &simplexes);
	double LUDB_experimental_quiet(TNode *root, int max_good_combos, uint64 &simplexes);
	double LUDB_trivial();
	double LUDB_getLast() { return cached_ludb; };
	void LUDB_setLast(double ludb) { cached_ludb = ludb; };

	// LowerBound computation
	double LowerBound(bool quiet = false, uint64 combo_range_start = 0, uint64 combo_range_len = 0, float percentage = 100.0);

	// Optimized LowerBound computation with reduced number of combinations to try (UNSTABLE)
	double LowerBound_experimental(float percentage = 100.0);
	double LowerBound_experimental_quiet(uint64 &num_combos, float percentage = 100.0);
	bool compute_output_rates(uint tagged); // temporarily public method
	uint64 compute_lowerbound_rules(uint tagged, byte *burst_time); // temporarily public method

private:
	uint num_nodes;
	uint num_flows;
	uint tagged_flow;

	Node *nodes_arr;
	Flow *flows_arr;

	void ComputeNestingLevels();
	int nesting_level;

	TNode *nesting_tree;
	TNode* build_nesting_subtree(uint flow_id);

	double compute_delay(bool *combo);
	bool bit(uint64 x, uint n);
	bool compute_CDF_i(uint node_id, uint flow_id);

	bool isNodeCrossedByFlow(uint node, uint flow);
	bool isNodeEnteredByFlow(uint node, uint flow);
	bool isNodeLeftByFlow(uint node, uint flow);
	bool isInterferingFlow(uint node, uint flow, uint tagged);

	#define MAXDEPS 32
	typedef struct {
		uint dep[MAXDEPS];
	} depsse;
	typedef struct {
		uint64 cres;
		depsse cdep;
		uint cnode;
	} compstatus;
	bool cutcomp_slow(bool *res, uint ndeps, std::set<uint> nodedeps[], std::set<std::set<uint> > &allres);
	void cutcomp(compstatus *stat, depsse deps[], std::set<std::set<uint> > &allres);

	// Cached data
	double cached_ludb;
};

}
#endif
