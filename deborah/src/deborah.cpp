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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "deborah.h"
#include "tandem.h"
#include "curve.h"
#include "poly.h"
#include "rng.h"
#include "simplex.h"
#include "timing.h"


using namespace std;
using namespace deborah;


// Enable the following to split benchmarking for cuts computation and LUDB processing
//#define LUDB_NNESTED_EVAL_BENCHMARK


// Initial seed for Random Number Generator
static int RNG_SEED = 1000;
// Random Number Generator (global)
static RNG rng(RNG_SEED);


// T-Student coefficients
#define T_STUDENT_SIZE 30
double t_student95[T_STUDENT_SIZE]={12.7062,4.3027,3.1825,2.7765,2.5706,2.4469,2.3646,2.3060,2.2622,2.2281,2.2010,2.1788,2.1604,2.1448,
	2.1315,2.1199,2.1098,2.1009,2.0930,2.0860,2.0796,2.0739,2.0687,2.0639,2.0595,2.0555,2.0518,2.0484,2.0452,2.0423};


//////////////////////////////////////////////////////////////////////////////////////////////////////////


/*	Print command line syntax
*/
int print_help()
{
	printf("Usage: deborah <config file> [options]\n\n");
	printf("A tiny Network Calculus tool to compute LUDB and Lower-bounds in tandems.\n\n");
	printf("       Options:   --ludb, --noludb     :  enable/disable computation of LUDB\n");
	printf("                  --ludb-heuristic <N> :  compute LUDB using the heuristic algorithm (max combos returned = N)\n");
	printf("                  --ludb-evaluate <tree_order> <tree_depth> <max_good_combos>\n");
	printf("                                       :  run LUDB computations on a balanced tree tandem with random provisioning\n");
	printf("                  --ludb-nnested-evaluate <N> <F>\n");
	printf("                                       :  run LUDB on a non-nested tandem with N nodes and F%% flows out of the max possible\n");
	printf("                  --ludb-nnested-sta   :  use STA algorithm to compute LUDB for non-nested tandems\n");
	printf("                  --ludb-nnested-single:  same as --ludb-nnested-sta\n");
	printf("                  --ludb-cuts-len  <N> :  limit max cuts length to <min_len+N> nodes\n");
	printf("                  --per-node           :  compute a delay bound using per-node analysis\n");
	printf("                  --lb, --nolb         :  enable/disable computation of Lower-bound\n");
	printf("                  --lb-experimental    :  turn on experimental features in Lower-bound computation (UNSTABLE)\n");
	printf("                  --lb-evaluate <tree_order> <tree_depth>\n");
	printf("                                       :  run Lower-bound computations on a balanced tree tandem with random provisioning\n");
	printf("                  --lb-random-combo <P>:  percentage of combos to be computed for Lower-bound (default 100.0)\n");
	printf("                  --tagged <N>         :  assume flow #N as tagged flow (by default it's the longest one)\n");
	printf("                  --scale-rates <F> <N>:  scale provisioned flow and node rates by F and N times, respectively\n");
	printf("                  --gen-tree <tree_order> <tree_depth> <filename> :  generate balanced-tree scenario config file\n");
	printf("                                          if <tree_order> is 1, then a sink-tree topology is generated.\n");
	printf("                  --gen-nnested <nodes> <flows percentage> <filename> :  generate non-nested tandem config file\n");
	printf("                  --rng-seed <N>       :  set seed for internal random number generator.\n");
	printf("\n");
	return 1;
}


// Terminate application immediately with an error message and code
void terminate(char *msg, int err_code)
{
	if(msg != NULL) printf("Error: %s\n",msg);
	else printf("Unspecified error, aborting.\n");
	exit(err_code);
}


/*	Dumb compare function between LUDB and LowerBound
*/
double compare_results(double ludb, double lowerbound, double ludb_pn, bool c_ludb, bool c_lowerbound)
{
	double delta;
	double rob;
	bool c_pernode = (ludb_pn > 0.0);

	if(!c_ludb && !c_lowerbound && !c_pernode) return 0.0;

	printf("DEBORAH results:\n-------------------------\n");
	delta = ludb - lowerbound;
	if(c_ludb) printf("      LUDB = %lf\n",ludb);
	if(c_lowerbound) printf("LowerBound = %lf\n",lowerbound);
	if(c_pernode) printf("   PN-LUDB = %lf\n", ludb_pn);
	if(c_ludb && c_pernode) printf("  PN-Ratio = %lf\n", ludb_pn/ludb);
	if(!c_ludb || !c_lowerbound) return 0.0;
	printf("     delta = %lf\n",delta);
	rob = 1.0 - (lowerbound / ludb);
	printf("       ROB = %lf\n", rob);
	if((delta >= 0.0) && (delta < CURVE_EPSILON))
		printf("The bounds computed are tight ;-)\n");
	else if(delta < 0.0)
		printf("Ooops, LUDB is lower than LowerBound :-[\n");
	return delta;
}


/*	This method computes the confidence interval for the 95th percentile
	Stores the mean value and the semi-width of the confidence interval (delta)
*/
bool confidence95(const double *vett, const int n, double &mean, double &delta)
{
	double var;
	double tstud_coeff;

	if(n<2) return false;
	if(n-1 >= T_STUDENT_SIZE) tstud_coeff = t_student95[T_STUDENT_SIZE-1];
	else tstud_coeff = t_student95[n-2];

	// compute mean
	mean = 0.0;
	for(int i=0; i<n; i++) mean += vett[i];
	mean /= n;

	// compute variance
	var = 0.0;
	for(int i=0; i<n; i++) var += pow(vett[i]-mean, 2.0);
	var /= n-1;

	// compute semi-interval
	delta = tstud_coeff * sqrt(var/n);

	return true;
}


/*	This method loads an array of LUDB values and simplex counts
	to avoid recomputing some stats
*/
int LUDB_LoadPrecomputedValues(char *fname, double *ludb, uint64 *simplexes, int num_trials)
{
	FILE *f;
	int trial;
	int n;
	for(int i=0; i<num_trials; i++)
	{
		ludb[i] = 0.0;
		simplexes[i] = 0;
	}
	f = fopen(fname,"rb");
	if(f==NULL) return -1;
	trial = 0;
	while(trial < num_trials)
	{
		n = fread(&ludb[trial],sizeof(double),1,f);
		if(n<=0) break;
		n = fread(&simplexes[trial],sizeof(uint64),1,f);
		if(n<=0) break;
		trial++;
	}
	fclose(f);
	return trial;
}


/*	This method dumps an array of LUDB values and simplex counts to disk, for reuse
*/
bool LUDB_SavePrecomputedValues(char *fname, double *ludb, uint64 *simplexes, int num_trials)
{
	FILE *f;
	f = fopen(fname,"wb");
	if(f == NULL) return false;
	for(int i=0; i<num_trials; i++)
	{
		fwrite(&ludb[i],sizeof(double),1,f);
		fwrite(&simplexes[i],sizeof(uint64),1,f);
	}
	fclose(f);
	return true;
}


/*	This method evaluates the accuracy of the LUDB-heuristic method compared to
	the LUDB-exact algorithm.
	A number of runs are performed on a balanced-tree topology.
*/
bool LUDB_evaluation(Tandem &tandem, int order, int depth, int max_good_combos)
{
	int num_trials = 50;
	bool compute_exact = true;
	int mgc_limit = 100;
	char precomp_fname[80];
	int ludb_exact_preload;

	uint64 simplexes_exact[num_trials];
	double ludb_exact[num_trials];
	double ludb_diff[num_trials];
	uint64 ludb_simplexes[num_trials];
	double ludb_simplexes_percent[num_trials];
	double mean, delta;

	balanced_tree_descriptor bt_desc;
	bt_desc.order = order;
	bt_desc.depth = depth;
	bt_desc.flow_burst_min = 100.0;
	bt_desc.flow_burst_max = 1000.0;
	bt_desc.flow_rate_min = 10.0;
	bt_desc.flow_rate_max = 100.0;
	bt_desc.node_rate_overprov_min = 1.01;
	bt_desc.node_rate_overprov_max = 1.50;
	bt_desc.node_latency_min = 0.0;
	bt_desc.node_latency_max = 0.0;

	tandem.CreateTree(bt_desc, rng);
	if(tandem.NumFlows() > NMAX)
	{
		printf("LUDB_Evaluation() error: too many t-nodes (%u), max %u allowed by simplex library.\n",tandem.NumFlows(),NMAX);
		return false;
	}
	tandem.BuildNestingTree();

	printf("LUDB Evaluation():  performing %i runs per simulation.\n", num_trials);
	printf("max_good_combos=%i\n",max_good_combos);
	printf("order=%i, depth=%i  --> %u flows, max %Lu simplexes\n",order,depth,tandem.NumFlows(),tandem.LUDB_TotalCombos());
	printf("flow_burst: %lf - %lf\n",bt_desc.flow_burst_min,bt_desc.flow_burst_max);
	printf("flow_rate: %lf - %lf\n",bt_desc.flow_rate_min,bt_desc.flow_rate_max);
	printf("node_latency: %lf - %lf\n",bt_desc.node_latency_min,bt_desc.node_latency_max);
	printf("node_rate_overprov: %.2lfx - %.2lfx\n",bt_desc.node_rate_overprov_min,bt_desc.node_rate_overprov_max);
	printf("\n");

	sprintf(precomp_fname,"ludb-exact-o%i-d%i-seed%i.dat",order,depth,RNG_SEED);
	ludb_exact_preload = LUDB_LoadPrecomputedValues(precomp_fname, ludb_exact, simplexes_exact, num_trials);
	if(ludb_exact_preload > 0)
		printf("Loaded precomputed LUDB values for %i runs (%s).\n",ludb_exact_preload,precomp_fname);

	if(tandem.LUDB_TotalCombos() > 1000000) compute_exact = false;

	rng = RNG(RNG_SEED);
	if(compute_exact)
	{
		printf("LUDB_exact\tSimplexes_Exact\n");
	}
	else {
		printf("LUDB_MGC%i\tSimplexes_MGC%i\n",mgc_limit,mgc_limit);
	}
	for(int trial=0; trial<num_trials; trial++)
	{
		tandem.CreateTree(bt_desc, rng);
		tandem.BuildNestingTree();
		if(simplexes_exact[trial] == 0)
		{
			if(compute_exact)
				ludb_exact[trial] = tandem.LUDB_experimental_quiet(LUDB_MODE_EXACT, simplexes_exact[trial]);
			else ludb_exact[trial] = tandem.LUDB_experimental_quiet(mgc_limit, simplexes_exact[trial]);
		}
		printf("%.12lf\t%Lu\n",ludb_exact[trial],simplexes_exact[trial]);
	}
	printf("\n");
	if(ludb_exact_preload < num_trials)
		if(LUDB_SavePrecomputedValues(precomp_fname,ludb_exact,simplexes_exact,num_trials))
			printf("Saved precomputed LUDB values for %i runs (%s)\n",num_trials,precomp_fname);

	printf("Max_GC\tMean_Err%%\tInt_semiwidth\tMax_Err%%\tMatch%%\t\tMean_Simpl%%\tInt_semiwidth\n");
	for(int k=0; k<max_good_combos+1; k++)
	{
		int nc;
		double ludb_heuristic;
		double ludb_diff_max = 0.0;
		rng = RNG(RNG_SEED);
		uint matches = 0;
		for(int trial=0; trial<num_trials; trial++)
		{
			tandem.CreateTree(bt_desc, rng);
			tandem.BuildNestingTree();
			if(compute_exact) nc = (k<max_good_combos) ? k+1 : LUDB_MODE_HEURISTIC;
			else nc = (k<max_good_combos) ? k+1 : mgc_limit;
			ludb_heuristic = tandem.LUDB_experimental_quiet(nc, ludb_simplexes[trial]);
			ludb_diff[trial] = (ludb_heuristic - ludb_exact[trial]) / ludb_exact[trial];
			ludb_diff_max = _MAX(ludb_diff_max, ludb_diff[trial]);
			ludb_simplexes_percent[trial] = (double) ludb_simplexes[trial] / (double) simplexes_exact[trial];
			if(fabs(ludb_heuristic - ludb_exact[trial]) < TNODE_LUDB_EPSILON) matches++;
		}
		confidence95(ludb_diff, num_trials, mean, delta);
		printf("%03i\t%.8lf\t%.8lf\t%.8lf\t%.4lf\t\t",nc,mean,delta,ludb_diff_max,(double) matches / (double) num_trials);
		confidence95(ludb_simplexes_percent, num_trials, mean, delta);
		printf("%.8lf\t%.8lf\n",mean,delta);
	}
	return true;
}


/*	This method evaluates the lower-bound delay using either the exhaustive or the reduced method.
	A number of runs are performed on a balanced-tree topology.
*/
bool LowerBound_evaluation(Tandem &tandem, int order, int depth, bool experimental, float percentage)
{
	int num_trials = 50;

	double lb_exact[num_trials];
	uint64 lb_combos[num_trials];
	uint64 lb_maxcombos;
	double mean, lb_max;
	char fname[80];
	uint max_flows = 8 * sizeof(uint64);
	bool dump_config = false;

	balanced_tree_descriptor bt_desc;
	bt_desc.order = order;
	bt_desc.depth = depth;
	bt_desc.flow_burst_min = 100.0;
	bt_desc.flow_burst_max = 1000.0;
	bt_desc.flow_rate_min = 10.0;
	bt_desc.flow_rate_max = 100.0;
	bt_desc.node_rate_overprov_min = 1.01;
	bt_desc.node_rate_overprov_max = 1.50;
	bt_desc.node_latency_min = 0.1;
	bt_desc.node_latency_max = 0.1;

	tandem.CreateTree(bt_desc, rng);

	if(tandem.NumFlows() >= max_flows)
	{
		printf("LowerBound_Evaluation() error: too many flows (%u), max %i allowed.\n",tandem.NumFlows(),max_flows);
		return false;
	}

	printf("LowerBound Evaluation():  performing %i run(s) per simulation using %s algorithm.\n", num_trials,(experimental?"experimental":"exhaustive"));
	printf("order=%i, depth=%i  --> %u flows\n",order,depth,tandem.NumFlows());
	printf("flow_burst: %lf - %lf\n",bt_desc.flow_burst_min,bt_desc.flow_burst_max);
	printf("flow_rate: %lf - %lf\n",bt_desc.flow_rate_min,bt_desc.flow_rate_max);
	printf("node_latency: %lf - %lf\n",bt_desc.node_latency_min,bt_desc.node_latency_max);
	printf("node_rate_overprov: %.2lfx - %.2lfx\n",bt_desc.node_rate_overprov_min,bt_desc.node_rate_overprov_max);
	printf("random combinations percentage: %.3f %%\n",percentage);
	printf("\n");

	// determine the total number of combinations that should be evaluated
	lb_maxcombos = (uint64) pow(2.0, (double)tandem.NumFlows());
	mean = 0.0;
	lb_max = 0.0;

	// initialize random number generator
	rng = RNG(RNG_SEED);

	printf("Trial\tLowerBound\tCombos\tMaxCombos\n");
	for(int trial=0; trial<num_trials; trial++)
	{
		tandem.CreateTree(bt_desc, rng);
		if(dump_config)
		{
			// save the generated tandem configuration for later debugging purposes
			sprintf(fname,"lb-o%i-d%i_%02i.conf",bt_desc.order,bt_desc.depth,trial);
			tandem.Save(fname);
		}
		// build nesting tree to initialize required data structures
		tandem.BuildNestingTree();
		if(experimental)
			lb_exact[trial] = tandem.LowerBound_experimental_quiet(lb_combos[trial], percentage);
		else {
			lb_exact[trial] = tandem.LowerBound(true, 0, 0, percentage);
			lb_combos[trial] = (uint64) pow(2.0, (double)tandem.NumFlows());
		}
		printf("%i\t%lf\t%Lu\t%Lu\n",trial,lb_exact[trial],lb_combos[trial],lb_maxcombos);

		// update some simple stats
		mean += lb_exact[trial];
		if(lb_exact[trial] - lb_max >= CURVE_EPSILON) lb_max = lb_exact[trial];
	}
	mean = mean / num_trials;
	printf("\n");
	printf("Mean=\t%lf	Max=\t%lf\n",mean,lb_max);

	return true;
}


/*	This method saves a configuration file for a tandem representing
	a balanced tree with the given parameters
*/
bool MakeTreeConfig(char *filename, const int order, const int depth)
{
	Tandem t0;
	bool res;
	balanced_tree_descriptor bt_desc;
	bt_desc.order = order;
	bt_desc.depth = depth;
	bt_desc.flow_burst_min = 100.0;
	bt_desc.flow_burst_max = 1000.0; //150.0;
	bt_desc.flow_rate_min = 10.0;
	bt_desc.flow_rate_max = 100.0; //70.0;
	bt_desc.node_rate_overprov_min = 1.01;
	bt_desc.node_rate_overprov_max = 1.50; //1.1;
	bt_desc.node_latency_min = 0.0;
	bt_desc.node_latency_max = 0.0;

	rng = RNG(RNG_SEED);

	res = t0.CreateTree(bt_desc, rng);
	if(!res) return false;
	res = t0.Save(filename);
	return res;
}


bool MakeFullTandemConfig(char *filename, const int num_nodes, const int flows_percent)
{
	Tandem tandem;
	bool res;
	balanced_tree_descriptor bt_desc;
	bt_desc.order = 0;
	bt_desc.depth = flows_percent;	// allocate nn% of all possible flows [N*(N+1)/2]
	bt_desc.flow_burst_min = 100.0;
	bt_desc.flow_burst_max = 1000.0;
	bt_desc.flow_rate_min = 10.0;
	bt_desc.flow_rate_max = 100.0;
	bt_desc.node_rate_overprov_min = 1.01;
	bt_desc.node_rate_overprov_max = 1.50;
	bt_desc.node_latency_min = 0.0;
	bt_desc.node_latency_max = 0.0;

	rng = RNG(RNG_SEED);

	bt_desc.order = num_nodes;
	int retries = 10;
	do {
		tandem.CreateFullTandem(bt_desc, rng);
		retries--;
	} while (tandem.isNested() && retries>0);
	if(retries == 0) return false;

	res = tandem.Save(filename);
	return res;
}


/*	Compares c1 with c2 and returns the index of the first cut
 *	where the two sets diverge
 */
uint compare_cutsets(CutNodelist &c1, CutNodelist &c2)
{
	uint l1 = c1.size();
	uint l2 = c2.size();
	uint len = _MIN(l1, l2);
	for(uint i=0; i<len; i++)
	{
		if(c1[i] != c2[i])
			return i;
	}
	return len;
}


//#define LUDB_NON_NESTED_DEBUG 1

/*	Computes the delay bound for the non-nested tandem 't1'
 * 	'pre_cset' points to a pre-computed set of cuts; if NULL, an exhaustive computation will be performed
 *  'num_cuts' is initialized so that cut sets longer than (min_len + num_cuts) will be skipped
 */
double LUDB_non_nested(Tandem &t1, int mode, std::vector<CutNodelist> *pre_cset, uint64 &totsimpl, uint &num_cuts)
{
	double min_ludb = TNODE_LUDB_INFINITY;
	double ludb;
	totsimpl = 0;
	int subtandems_total=0, subtandems_reused=0;
	std::vector<uint> best_results;

	// compute all possible sets of cuts
	std::vector<CutNodelist> cset;
	if(pre_cset == NULL) cset = t1.compute_cutsets();
	else cset = *pre_cset;
	//if(num_cuts != TANDEM_LUDB_MAXCUTS) t1.filter_cutsets(cset, num_cuts);
	if(num_cuts != TANDEM_LUDB_MAXCUTS) t1.filter_cutsets_v2(cset, num_cuts);
	num_cuts = cset.size();

	// remember the last set of cuts to reuse LUDB computations
	std::vector<Tandem> last_cuts;
	uint cuts_reuse = 0; // number of sub-tandem which could be reused from the previous step

	// sanity check of set of cuts found
	uint nsets = cset.size();
	if(nsets==0)
	{
		printf("Fatal error: the tandem is non-nested but no cuts have been found (bug?)\n");
		return DEBORAH_ERROR_SAFECHECK;
	}

	if(pre_cset == NULL) printf("Number of sets of cuts to try: %i\n", nsets);

	// for each set of cuts, compute end-to-end delay as the sum of intermediate LUDBs
	for(uint sel_set=0; sel_set<nsets; sel_set++)
	{
		ludb = 0.0;
		if(sel_set > 0)
		{
			// compare with the last cutset
			cuts_reuse = compare_cutsets(cset[sel_set], cset[sel_set-1]);
			if(cuts_reuse == cset[sel_set].size())
			{
				// safety check
				printf("*** CUTS_REUSE RANGE EXCEPTION ***\n");
				cuts_reuse = 0;
			}
		}
		else cuts_reuse = 0;
		std::vector<Tandem> tmp_cuts;
		tmp_cuts.clear();

		// cut the original tandem according to the i-th set of cuts
		std::vector<Tandem> cuts = t1.cut(cset[sel_set]);

		// update stats about computation savings
		subtandems_total += cuts.size();
		subtandems_reused += cuts_reuse;

#ifdef LUDB_NON_NESTED_DEBUG
		printf("\n");
		printf("Cut set #%02i contains %i sub-tandem(s)  (reusable=%i):\n", sel_set, cuts.size(), cuts_reuse);
		for(uint i=0; i<cuts.size(); i++) cuts[i].Print();
#endif

		// for each sub-tandem, compute 1) the LUDB, and 2) the relevant output arrival curves
		for(uint i=0; i<cuts.size(); i++)
		{
			double ti_ludb;
#if LUDB_NON_NESTED_DEBUG > 1
			printf("********************************** sub-tandem %i: **********************************\n",i);
#endif
			// check if we can reuse numeric results from the previous cut
			if(i < cuts_reuse)
			{
				// retrieve cached LUDB values
				Tandem ti = last_cuts[i];
				ti_ludb = ti.LUDB_getLast();
#ifdef LUDB_NON_NESTED_DEBUG
				printf("T%i.LUDB_cached = %.3lf\n", i, ti_ludb);
#endif
				ludb += ti_ludb;
				tmp_cuts.push_back(ti);
				// nothing else to do, continue with next sub-tandem
				continue;
			}
			else if((i > 0) && (i == cuts_reuse))
			{
				// for this tandem we have to recompute LUDB using the proper O.A.Curves
				Tandem &t1 = cuts[i];
				Tandem &t2 = last_cuts[i];
				// copy output arrival curves from the previous cut
				for(uint f=0; f<t1.NumFlows(); f++)
				{
					Flow *fp1 = t1.getFlow(f);
					uint fid = fp1->uid;
					if(fid == FLOW_INVALID_ID) continue;
					Flow *fp2 = t2.getFlowId(fid);
					if(fp2 == NULL) continue;
					fp1->burst = fp2->burst;
					fp1->rate = fp2->rate;
				}
			}

			// compute LUDB for the tagged flow
			Tandem ti = cuts[i];
			// save tandem
			tmp_cuts.push_back(ti);
			// aggregate flows to the maximum extent, including the tagged flow
			ti.aggregate_flows(true);
			ti.BuildNestingTree();
			//ti.Print(); ti.PrintNestingTree();
			uint64 nsimpl = 0;
			// first try with the trivial cases
			ti_ludb = ti.LUDB_trivial();
			// if not a supported case, perform normal LUDB computation
			if(ti_ludb < 0.0)
				ti_ludb = ti.LUDB_experimental_quiet(mode, nsimpl);
			tmp_cuts[i].LUDB_setLast(ti_ludb);
#ifdef LUDB_NON_NESTED_DEBUG
			printf("T%i.LUDB = %.3lf\n", i, ti_ludb);
#endif
			ludb += ti_ludb; // add the sub-tandem's delay to the global LUDB
			totsimpl += nsimpl;

			// compute output arrival curves
			if(i == cuts.size()-1) continue; // not for the last sub-tandem
#if LUDB_NON_NESTED_DEBUG > 1
			printf(">>> compute output arrival curves at exit of Tandem #%i\n",i);
#endif
			cuts[i].output_arrival_curves(cuts[i+1], mode); // recursive helper function that does the job
#if LUDB_NON_NESTED_DEBUG > 1
			printf("<<< done\n\n\n");
#endif
		}
		if(pre_cset == NULL)
		{
			printf("Cuts set %03i  ==>  LUDB= %.04lf \t cuts_len= %i \tcuts={ ", sel_set, ludb,cset[sel_set].size());
			for(uint c=0; c<cset[sel_set].size();c++)
				printf("%i ",cset[sel_set][c]);
			printf("}\n");
		}
		// keep track of the minimum LUDB found amongst all possible cut sets
		//min_ludb = _MIN(ludb, min_ludb);
		if(ludb <= min_ludb) {
			if(ludb < min_ludb) {
				best_results.clear();
				min_ludb = ludb;
			}
			best_results.push_back(sel_set);
		}

		// commit saved data for reuse
		last_cuts.clear();
		last_cuts = tmp_cuts;
	}

	if((nsets > 1) && (pre_cset == NULL))
	{
		double reusage = (double) subtandems_reused / (double) subtandems_total;
		printf("Sub-tandem computations reusage ratio = %.4lf  (%i / %i)\n", reusage, subtandems_reused, subtandems_total);
	}

    // print the IDs of the cuts that correspond to the minimum LUDB found
    printf("Best set(s) of cuts:  { ");
    for(uint c=0; c<best_results.size(); c++)
    	printf("%03i ", best_results[c]);
    printf("}\n");

	return min_ludb;
}

/*	Computes LUDB for the non-nested tandem 't1' using an alternative method:
 * 	For each set of cuts:
 * 	- compute arrival curves at cuts' boundaries
 *      - compute LUDB for the whole tandem (rather than the sum of per-cut LUDBs)
 */
double LUDB_non_nested_single(Tandem &t1, int mode, std::vector<CutNodelist> *pre_cset, uint64 &totsimpl, uint &num_cuts)
{
	double min_ludb = TNODE_LUDB_INFINITY;
	double ludb;
	totsimpl = 0;
	int subtandems_total=0, subtandems_reused=0;
	std::vector<uint> best_results;

	// compute all possible sets of cuts
	std::vector<CutNodelist> cset;
	if(pre_cset == NULL) cset = t1.compute_cutsets();
	else cset = *pre_cset;
	//if(num_cuts != TANDEM_LUDB_MAXCUTS) t1.filter_cutsets(cset, num_cuts);
	if(num_cuts != TANDEM_LUDB_MAXCUTS) t1.filter_cutsets_v2(cset, num_cuts);
	num_cuts = cset.size();

	// remember the last set of cuts to reuse LUDB computations
	std::vector<Tandem> last_cuts;
	uint cuts_reuse = 0; // number of sub-tandem which could be reused from the previous step

	// sanity check of set of cuts found
	uint nsets = cset.size();
	if(nsets==0)
	{
		printf("Fatal error: the tandem is non-nested but no cuts have been found (bug?)\n");
		return DEBORAH_ERROR_SAFECHECK;
	}

	if(pre_cset == NULL) printf("Number of sets of cuts to try: %i\n", nsets);

	// for each set of cuts, compute the global end-to-end delay
	for(uint sel_set=0; sel_set<nsets; sel_set++)
	{
		ludb = 0.0;
		if(sel_set > 0)
		{
			// compare with the last cutset
			cuts_reuse = compare_cutsets(cset[sel_set], cset[sel_set-1]);
			if(cuts_reuse == cset[sel_set].size())
			{
				// safety check
				printf("*** CUTS_REUSE RANGE EXCEPTION ***\n");
				cuts_reuse = 0;
			}
		}
		else cuts_reuse = 0;
		std::vector<Tandem> tmp_cuts;
		tmp_cuts.clear();

		//cuts_reuse = 0;

		// cut the original tandem according to the i-th set of cuts
		std::vector<Tandem> cuts = t1.cut(cset[sel_set]);

		// update stats about computation savings
		subtandems_total += cuts.size();
		subtandems_reused += cuts_reuse;

#ifdef LUDB_NON_NESTED_DEBUG
		printf("\nUsing alternative algorithm (single LUDB over joined tandem):\n");
		printf("Cut set #%02i contains %i sub-tandem(s)  (reusable=%i):\n", sel_set, cuts.size(), cuts_reuse);
		for(uint i=0; i<cuts.size(); i++) cuts[i].Print();
#endif

		// for each sub-tandem, compute the relevant output arrival curves
		for(uint i=0; i<cuts.size(); i++)
		{
#if LUDB_NON_NESTED_DEBUG > 1
			printf("********************************** sub-tandem %i: **********************************\n",i);
#endif
			// check if we can reuse numeric results from the previous cut
			if(i < cuts_reuse)
			{
				// retrieve cached LUDB values
				Tandem ti = last_cuts[i];
				tmp_cuts.push_back(ti);
				// nothing else to do, continue with next sub-tandem
				continue;
			}
			else if((i > 0) && (i == cuts_reuse))
			{
				// for this tandem we have to recompute LUDB using the proper O.A.Curves
				Tandem &t1 = cuts[i];
				Tandem &t2 = last_cuts[i];
				// copy output arrival curves from the previous cut
				for(uint f=0; f<t1.NumFlows(); f++)
				{
					Flow *fp1 = t1.getFlow(f);
					uint fid = fp1->uid;
					if(fid == FLOW_INVALID_ID) continue;
					Flow *fp2 = t2.getFlowId(fid);
					if(fp2 == NULL) continue;
					fp1->burst = fp2->burst;
					fp1->rate = fp2->rate;
				}
			}

			Tandem ti = cuts[i];
			tmp_cuts.push_back(ti);

			// compute new output arrival curves
			if(i == cuts.size()-1) continue; // not for the last sub-tandem
#if LUDB_NON_NESTED_DEBUG > 1
			printf(">>> compute output arrival curves at exit of Tandem #%i\n",i);
#endif
			cuts[i].output_arrival_curves(cuts[i+1], mode); // recursive helper function that does the job
#if LUDB_NON_NESTED_DEBUG > 1
			printf("<<< done\n\n\n");
#endif
		}
		// now re-assemble sub-tandems into a global nested tandem
		Tandem tglob;
		if(!tglob.join(tmp_cuts))
		{
			printf("Internal error: joined tandem is not nested, skipping!\n");
			last_cuts.clear();
			last_cuts = tmp_cuts;
			continue;
		}

		// compute LUDB for the global
		tglob.aggregate_flows(false);
		tglob.BuildNestingTree();
#if LUDB_NON_NESTED_DEBUG > 0
		printf("Reassembled tandem:\n");
		tglob.Print();
#endif
		uint64 nsimpl = 0;
		ludb = tglob.LUDB_trivial();
		if(ludb < 0.0) ludb = tglob.LUDB_experimental_quiet(mode, nsimpl);

		// LUDB done
		if(pre_cset == NULL)
		{
			printf("Cuts set %03i  ==>  LUDB= %.04lf \t cuts_len= %i \tcuts={ ", sel_set, ludb,cset[sel_set].size());
			for(uint c=0; c<cset[sel_set].size();c++)
				printf("%i ",cset[sel_set][c]);
			printf("}\n");
		}
		// keep track of the minimum LUDB found amongst all possible cut sets
		//min_ludb = _MIN(ludb, min_ludb);
		if(ludb <= min_ludb) {
			if(ludb < min_ludb) {
				best_results.clear();
				min_ludb = ludb;
			}
			best_results.push_back(sel_set);
		}

		// commit saved data for reuse
		last_cuts.clear();
		last_cuts = tmp_cuts;
	}

	if((nsets > 1) && (pre_cset == NULL))
	{
		double reusage = (double) subtandems_reused / (double) subtandems_total;
		printf("Sub-tandem computations reusage ratio = %.4lf  (%i / %i)\n", reusage, subtandems_reused, subtandems_total);
	}

	// print the IDs of the cuts that correspond to the minimum LUDB found
	printf("Best set(s) of cuts:  { ");
	for(uint c=0; c<best_results.size(); c++)
	  printf("%03i ", best_results[c]);
	printf("}\n");

	return min_ludb;
}


/*	Run a series of simulations to assess the performance of the LUDB analysis
 * 	in non-nested tandems.
 * 	'nnodes' is the length of the tandem (number of nodes)
 * 	'flows_percent' is the percentage (1 - 100) of how many flows will be allocated
 * 	'mode' specifies the LUDB computation mode (exact or heuristic)
 * 	'maxcuts' limits the length of the analyzed cuts to 'min_length+maxcuts'
 */
bool LUDB_non_nested_evaluation(int nnodes, int flows_percent, int mode, uint maxcuts, bool algo_classic)
{
	balanced_tree_descriptor bt_desc;
	bt_desc.order = 0;
	bt_desc.depth = flows_percent;	// allocate nn% of all possible flows [N*(N+1)/2]
	bt_desc.flow_burst_min = 100.0;
	bt_desc.flow_burst_max = 1000.0;
	bt_desc.flow_rate_min = 10.0;
	bt_desc.flow_rate_max = 100.0;
	bt_desc.node_rate_overprov_min = 1.01;
	bt_desc.node_rate_overprov_max = 1.50;
	bt_desc.node_latency_min = 0.0;
	bt_desc.node_latency_max = 0.0;

	//BasicTimer t0;
	NanoTimer t0;

	printf("LUDB_non_nested_evaluation(): tandem length = %i nodes, %i%% flow coverage, algorithm = %s\n",nnodes, bt_desc.depth, (algo_classic?"MTA":"STA"));
	printf("                              LUDB computation algorithm = ");
	if(mode == LUDB_MODE_EXACT) printf("exact\n");
	else if(mode == LUDB_MODE_HEURISTIC) printf("heuristic, all results\n");
	else printf("heuristic, best %i results\n", mode);
	//fprintf(stderr, "num_nodes,perc_flows,num_flows,num_cutsets,time,ludb,num_simplexes\n\n");
	rng = RNG(RNG_SEED);
	Tandem tandem;
	uint64 nsimpl=0;
	uint num_cuts=maxcuts;
	double ludb = 0.0;
	double t_compcuts, t_ludb;

	for(int n=nnodes; n<=nnodes; n++)
	{
		bt_desc.order = n;
		int retries = 10;
		do {
			tandem.CreateFullTandem(bt_desc, rng);
			retries--;
		} while (tandem.isNested() && retries>0);
		uint nflows = tandem.NumFlows();
		//tandem.getFlow(1)->burst *= 200;
		tandem.Print();

#ifdef LUDB_NNESTED_EVAL_BENCHMARK
		t0.mark(true);
		std::vector<CutNodelist> cuts = tandem.compute_cutsets();
		t_compcuts = t0.mark(false);
		num_cuts = _MIN(cuts.size(), num_cuts);

		t0.mark(true);
		// Comment next if() to skip LUDB and measure cuts computation time only
		//if(algo_classic) ludb = LUDB_non_nested(tandem, mode, &cuts, nsimpl, num_cuts);
		//else ludb = LUDB_non_nested_single(tandem, mode, &cuts, nsimpl, num_cuts);
		ludb=1.0;
		t0.mark(false);

		for(int i=0; i<num_cuts; i++)
		{
			printf("set #%02i: { ", i);
			for(int j=0; j<cuts[i].size(); j++)
				printf("%i ", cuts[i][j]);
			printf("}\t len=%i\n", cuts[i].size());
		}
		//fprintf(stderr, "%i,%i,%i,%i,%.3lf,%.3lf,%lf,%Lu\n", n, flows_percent, nflows, num_cuts, t_compcuts,0.0,0.0,0);
		//continue;
#else
		t_compcuts = -1.0;
		t0.mark(true);
		if(algo_classic) ludb = LUDB_non_nested(tandem, mode, NULL, nsimpl, num_cuts);
		else ludb = LUDB_non_nested_single(tandem, mode, NULL, nsimpl, num_cuts);
		t_ludb = t0.mark(false);
		std::vector<CutNodelist> cuts;
#endif

		// compare with per-node analysis
		cuts.clear();
		CutNodelist allnodes;
		uint64 nsimpl2 = 0;
		uint ncuts2 = TANDEM_LUDB_MAXCUTS;
		for(uint i=1; i<=tandem.NumNodes(); i++) allnodes.push_back(i);
		cuts.push_back(allnodes);
		double ludb_pernode = LUDB_non_nested(tandem, mode, &cuts, nsimpl2, ncuts2);
		double ratio = ludb_pernode / ludb;
		printf("LUDB= %lf \tpnLUDB= %lf \tpnLUDB/LUDB= %lf\n", ludb, ludb_pernode, ratio);

		fprintf(stderr, "%i,%i,%i,%i,%.3lf,%.3lf,%lf,%Lu\n", n, flows_percent, nflows, num_cuts, t_compcuts,t_ludb,ratio,nsimpl);
	}
	return true;
}


/*	Compute a delay bound for tandem 't1' using per-node analysis
 */
double DelayBound_per_node(Tandem &t1, int ludb_mode)
{
	std::vector<CutNodelist> cuts;
	CutNodelist allnodes;
	uint64 nsimpl = 0;
	uint ncuts = TANDEM_LUDB_MAXCUTS;
	for(uint i=1; i<=t1.NumNodes(); i++) allnodes.push_back(i);
	cuts.push_back(allnodes);
	double ludb_pernode = LUDB_non_nested(t1, ludb_mode, &cuts, nsimpl, ncuts);
	return ludb_pernode;
}


/* ******************************************************************************************************* */
/* **********************************   DEBORAH application main code   ********************************** */
/* ******************************************************************************************************* */


static Tandem t1;


int main(int argc, char *argv[])
{
	char conf_file[260]="";
	bool res;
	int ires;
	double lowerbound=0.0, ludb=0.0, ludb_pn=0.0;
	int TaggedFlow = -1;

	bool compute_ludb = false;
	bool compute_ludb_heuristic = false;
	bool compute_ludb_nnested_alt = false;
	int ludb_max_good_combos = LUDB_MODE_HEURISTIC;
	int ludb_tree_order = 2;
	int ludb_tree_depth = 3;
	bool ludb_evaluation = false;
	bool ludb_nnested_evaluation = false;
	bool compute_ludb_pernode = false;
	bool lowerbound_evaluation = false;
	bool compute_lowerbound = false;
	bool lowerbound_experimental = false;
	float lowerbound_percentage = 100.0;
	uint ludb_maxcuts = TANDEM_LUDB_MAXCUTS;
	double flowrate_mult = 1.0, noderate_mult = 1.0;

	uint64 combo_range_start = 0, combo_range_len = 0;

	cout << "--=[ DEBORAH v" << DEBORAH_VERSION << " ]=--   (C) 2008-2011 Dip. Ing. dell'Informazione, University of Pisa (Italy)   ";
	cout << "(build: " << __DATE__ << ", " << __TIME__ << ")"<< endl << endl;

	// parse command line arguments
	for(int i=1; i<argc; i++)
	{
		if(!strcasecmp(argv[i],"--noludb"))
				  compute_ludb = false;
		else if(!strcasecmp(argv[i],"--nolb"))
				  compute_lowerbound = false;
		else if(!strcasecmp(argv[i],"--lb"))
			compute_lowerbound = true;
		else if(!strcasecmp(argv[i],"--lb-experimental"))
			lowerbound_experimental = true;
		else if(!strcasecmp(argv[i],"--ludb-evaluate"))
		{
			ludb_evaluation = true;
			if(i < argc-3)
			{
				sscanf(argv[i+1],"%i",&ludb_tree_order);
				sscanf(argv[i+2],"%i",&ludb_tree_depth);
				sscanf(argv[i+3],"%i",&ludb_max_good_combos);
				i += 3;
			}
		}
		else if(!strcasecmp(argv[i], "--ludb-nnested-evaluate"))
		{
			ludb_nnested_evaluation = true;
			ludb_max_good_combos = LUDB_MODE_EXACT;
			ludb_tree_depth = 100; // default: allocate all flows (100%)
			if(i < argc-2)
			{
				sscanf(argv[i+1],"%i",&ludb_tree_order);
				sscanf(argv[i+2],"%i",&ludb_tree_depth);
				i += 2;
			}
			if(ludb_tree_order < 3) ludb_tree_order = 3;
		}
		else if(!strcasecmp(argv[i],"--lb-evaluate"))
		{
			lowerbound_evaluation = true;
			if(i < argc-2)
			{
				sscanf(argv[i+1],"%i",&ludb_tree_order);
				sscanf(argv[i+2],"%i",&ludb_tree_depth);
				i += 2;
			}
		}
		else if((!strcasecmp(argv[i],"--lb-random-combo")) && (i < argc-1))
		{
#ifdef TANDEM_LB_RANDOM_COMBOS
			sscanf(argv[i+1],"%f",&lowerbound_percentage);
#endif
			i++;
		}
		else if(!strcasecmp(argv[i],"--gen-tree"))
		{
			if(i < argc-3)
			{
				sscanf(argv[i+1],"%i",&ludb_tree_order);
				sscanf(argv[i+2],"%i",&ludb_tree_depth);
				MakeTreeConfig(argv[i+3], ludb_tree_order, ludb_tree_depth);
			}
			printf("Balanced tree (order=%i, depth=%i) written to \"%s\"\n",ludb_tree_order,ludb_tree_depth,argv[i+3]);
			return 0;
		}
		else if(!strcasecmp(argv[i],"--gen-nnested"))
		{
			if(i < argc-3)
			{
				sscanf(argv[i+1],"%i",&ludb_tree_order);
				sscanf(argv[i+2],"%i",&ludb_tree_depth);
				MakeFullTandemConfig(argv[i+3], ludb_tree_order, ludb_tree_depth);
			}
			printf("Non-nested tandem (nodes=%i, flows percentage=%i) written to \"%s\"\n",ludb_tree_order,ludb_tree_depth,argv[i+3]);
			return 0;
		}
		else if(!strcasecmp(argv[i],"--combo-range"))
		{
			if(i < argc-2)
			{
				sscanf(argv[i+1],"%Lu",&combo_range_start);
				sscanf(argv[i+2],"%Lu",&combo_range_len);
				i += 2;
			}
		}
		else if(!strcasecmp(argv[i],"--ludb"))
		{
			compute_ludb = true;
			ludb_max_good_combos = LUDB_MODE_EXACT;
		}
		else if(!strcasecmp(argv[i],"--ludb-heuristic"))
		{
			compute_ludb_heuristic = compute_ludb = true;
			ludb_max_good_combos = -1;
			if(i < argc-1)
			{
				int res = sscanf(argv[i+1],"%i",&ludb_max_good_combos);
				if(res>0) i += res;
			}
			if(ludb_max_good_combos <= 0) ludb_max_good_combos = LUDB_MODE_HEURISTIC;
		}
		else if(!strcasecmp(argv[i], "--ludb-cuts-len"))
		{
			//ludb_maxcuts = 0;
			if(i < argc-1)
			{
				sscanf(argv[i+1], "%i", &ludb_maxcuts);
				i++;
			}
		}
		else if(!strcasecmp(argv[i],"--ludb-nnested-single") || !strcasecmp(argv[i],"--ludb-nnested-sta"))
		{
			compute_ludb_nnested_alt = true;
		}
		else if(!strcasecmp(argv[i], "--per-node"))
		{
			compute_ludb_pernode = true;
			ludb_max_good_combos = LUDB_MODE_EXACT;
		}
		else if((!strcasecmp(argv[i],"--tagged")) && (i < argc-1))
		{
			sscanf(argv[i+1],"%i",&TaggedFlow);
			i++;
		}
		else if((!strcasecmp(argv[i],"--threads")) && (i < argc-1))
		{
			//sscanf(argv[i+1],"%i",&NumThreads);
			i++;
		}
		else if((!strcasecmp(argv[i],"--rng-seed")) && (i < argc-1))
		{
			sscanf(argv[i+1],"%i",&RNG_SEED);
			i++;
		}
		else if(!strcasecmp(argv[i],"--scale-rates"))
		{
			if(i < argc-2)
			{
				sscanf(argv[i+1],"%lf",&flowrate_mult);
				sscanf(argv[i+2],"%lf",&noderate_mult);
				i += 2;
				if(flowrate_mult <= 0.0) flowrate_mult = 1.0;
				if(noderate_mult <= 0.0) noderate_mult = 1.0;
			}
		}
		else if(!strcasecmp(argv[i],"--help") || !strcasecmp(argv[i],"-h"))
			return print_help();
		// last cases: either an unrecognized option, or the scenario file has been specified
		else if(argv[i][0] == '-')
		{
			printf("Error: unrecognized option \"%s\"\n",argv[i]);
			printf("       Please run \"%s --help\" for a list of valid command line options.\n",argv[0]);
			return DEBORAH_ERROR_CLI;
		}
		else strncpy(conf_file,argv[1],sizeof(conf_file)-1);
	}

	if(ludb_evaluation)
	{
		LUDB_evaluation(t1, ludb_tree_order, ludb_tree_depth, ludb_max_good_combos);
		return 0;
	}
	if(lowerbound_evaluation)
	{
		LowerBound_evaluation(t1, ludb_tree_order, ludb_tree_depth, lowerbound_experimental, lowerbound_percentage);
		return 0;
	}
	if(ludb_nnested_evaluation)
	{
		LUDB_non_nested_evaluation(ludb_tree_order, ludb_tree_depth,
				ludb_max_good_combos, ludb_maxcuts, !compute_ludb_nnested_alt);
		return 0;
	}

	// check if a configuration file has been specified
	if(!strlen(conf_file))
	{
		printf("Error: a scenario configuration file must be specified (or any -xxxx-evaluate switch).\n");
		print_help();
		//t1.Test_LowerBound_code();
		//t1.PrintAllCombinations();
		//lowerbound = t1.LowerBound_experimental();
		return -1;
	}

	// load the configuration file
	printf("Parsing scenario configuration file \"%s\"... ",conf_file);
	res = t1.Load(conf_file);
	if(!res)
	{
		printf("error!\n");
		return DEBORAH_ERROR_CONFIG;
	}
	else printf("ok.\n");

	if((flowrate_mult != 1.0) || (noderate_mult != 1.0))
	{
		double rmin,rmax;
		int nmin, nmax;
		t1.analyzeProvisioning(rmin, nmin, rmax, nmax);
		printf("Current over-provisioning status: min=%.3lf at node %i,  max=%.3lf at node %i\n", rmin, nmin+1, rmax, nmax+1);
		printf("Rescaling rates for over-provisioning: flows ratio=%.3lf,  nodes ratio=%.3lf\n", flowrate_mult, noderate_mult);
		res = t1.scaleProvisioning(flowrate_mult, noderate_mult);
		if(!res)
		{
			printf("Error: Provisioning re-scaling parameters determine violation of constraints\n");
			return DEBORAH_ERROR_PROV;
		}
		t1.analyzeProvisioning(rmin, nmin, rmax, nmax);
		printf("New over-provisioning status: min=%.3lf at node %i,  max=%.3lf at node %i\n", rmin, nmin+1, rmax, nmax+1);
	}

	// check if provisioning at the nodes is correct
	ires = t1.checkProvisioning();
	if(ires < 0)
	{
		printf("Error: Provisioning constraints not respected at node %i\n",ires);
		return DEBORAH_ERROR_PROV;
	}

	// check if user has specified a tagged flow explicitly
	if(TaggedFlow >= 0)
		t1.SetTaggedFlow(TaggedFlow);

	t1.Print();

	// check whether the tandem is nested
	bool nested = t1.isNested();

	// build the nesting tree
	t1.BuildNestingTree();

	if(nested)
	{
		printf("\nAssociated nesting tree:\n");
		printf("----------------------------------------------------\n");
		t1.PrintNestingTree();
		printf("----------------------------------------------------\n");
	}

	// check if this tandem is equivalent to a sink-tree
	if(t1.isSinkTree())
		printf("This is a sink-tree.\n");

	printf("\n");

	// check if dimensional limits of the simplex library are respected
	if(compute_ludb && (t1.NumFlows() > NMAX))
	{
		printf("LUDB computation disabled: too many t-nodes (%u), max %u allowed by simplex library.\n",t1.NumFlows(),NMAX);
		printf("(see \"simplex.h\" in source code)\n");
		compute_ludb = false;
	}

	// before proceeding, make sure we actually have some job to do
	if(!compute_ludb && !compute_lowerbound && !compute_ludb_pernode)
	{
		printf("Nothing to do, exiting...\n\n");
		return 0;
	}

	// run the LUDB algorithm
	if(compute_ludb)
	{
		if(nested) ludb = t1.LUDB_experimental(ludb_max_good_combos);
		else {
			printf("\nThis tandem is not nested, using \"%s\" algorithm.\n",
					(compute_ludb_nnested_alt ? "STA" : "MTA"));
			uint64 nsimpl = 0;
			uint ncuts = ludb_maxcuts;
			if(compute_ludb_nnested_alt)
				ludb = LUDB_non_nested_single(t1, ludb_max_good_combos, NULL, nsimpl, ncuts);
			else ludb = LUDB_non_nested(t1, ludb_max_good_combos, NULL, nsimpl, ncuts);
		}
	}

	// run the LowerBound algorithm
	if(compute_lowerbound) lowerbound = t1.LowerBound(false, combo_range_start, combo_range_len, lowerbound_percentage);
	else if(lowerbound_experimental) lowerbound = t1.LowerBound_experimental(lowerbound_percentage);

	// run the per-node delay bound analysis
	if(compute_ludb_pernode) ludb_pn = DelayBound_per_node(t1, ludb_max_good_combos);

	// compare the delay values obtained from the two algorithms
	printf("\n");
	compare_results(ludb, lowerbound, ludb_pn, compute_ludb, compute_lowerbound || lowerbound_experimental);

	printf("\n");
	return 0;
}
