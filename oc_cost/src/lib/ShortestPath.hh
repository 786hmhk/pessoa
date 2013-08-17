/**
 * @file
 * @author Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>
 * @version 0.61
 *
 * @section LICENSE
 *
 * Copyright (c) 2013, TU Delft: Delft University of Technology
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of TU Delft: Delft University of Technology nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL TU DELFT: DELFT UNIVERSITY OF TECHNOLOGY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @section DESCRIPTION
 *
 * This file contains the ShortestPath Class that is used to construct the
 * DD's used in optimal control.
 *
 * 		No details yet.
 */


#ifndef SHORTESTPATH_H_
#define SHORTESTPATH_H_


#include <cstdio>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <map>


// C Libs
#include "cudd.h"
#include "util.h"
#include "cuddInt.h"
#include "dddmp.h"   // for storing DD into files.
// C++ Lib
#include "cuddObj.hh"

//! Enables or Disables the use of the cudd cache in some operations.
//#define ENABLE_CACHE
#define ENABLE_TIME_PROFILING


#ifdef ENABLE_TIME_PROFILING
#include <sys/time.h>
#endif

// Cache Tags
//! Cache Tag for the cache used in AddOuterSumRecurTrace method.
#define DD_ADD_OUT_SUM_TRACE_TAG		0x6f
//! Cache Tag for the cache used in cuddAddApplyRecurTrace method.
#define DD_ADD_MINIMUM_TRACE_TAG		0x6d

//!
typedef std::pair<double, unsigned int> pair_double_int;
//! Priority queue to store the state number with the minimum cost. Used in the Relax method.
typedef std::priority_queue<pair_double_int, std::vector<pair_double_int>, std::greater<pair_double_int> > pq_relax;

//!
struct bdd_constNode_pos {
  int position;
  DdNode *const_node;
};

typedef struct bdd_constNode_pos min_result;

//- Create the ADD of the target set W, given the set of states of W (@ref createTargetSet).

//! This is the main class for generating the DD's used in optimal control.
//! Supports only deterministic systems at the moment.
/*!
This class is used to construct the ADD containing the shortest path values
from all states to the target set W and the ADD containing the actual path.

This class can be used to do the following:
- Create the Cost Adjacency Matrix, given the System's BDD and the ADD containing the cost
of each state (@ref createCostAdjacencyMatrix). Supports only Deterministic Systems.
- Find the all-pair shortest path and the pointer array, given the Cost Adjacency Matrix
(@ref FloydWarshall).
- Find the shortest path form all-pairs to a given set W (@ref APtoSetSP). Supports only Deterministic Systems.
*/
class ShortestPath {

private:

	// C++ manager.
	Cudd *mgr_cpp;
	// C manager.
	DdManager *mgr;
	// Systems variables.
	unsigned int no_states;
	unsigned int no_inputs;
	unsigned int no_state_vars;
	unsigned int no_input_vars;
	std::vector<BDD> bdd_x;
	std::vector<BDD> bdd_u;
	std::vector<BDD> bdd_x_;
	std::vector<ADD> add_x;
//	std::vector<ADD> add_u;
	std::vector<ADD> add_x_;
	// If the System has been analyzed yet.
	bool system_analyzed;

	// System's BDD
	BDD *system_bdd;


	unsigned int getNoBits(unsigned int number);

	void AddOuterSumTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node);
	void AddOuterSumRecurTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node);

	void cuddAddMinimumRecur(DdNode ** f, DdNode ** g, DdNode * R, DdNode **Result, int node);
	void Cudd_addApplyMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node);
	void cuddAddApplyRecurMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node);

	void minimum2(ADD *f, ADD *g, ADD *P);
	static DdNode *Cudd_addMinimumNS(DdManager * dd, DdNode ** f, DdNode ** g);
	void Cudd_addApplyMin2(DdNode * f, DdNode * g, DdNode * Pf, DdNode * Pg, DdNode **Result);
	void cuddAddApplyMin2Recur(DD_AOP op, DdNode * f, DdNode * g, DdNode * Pf, DdNode * Pg, DdNode **Result);

	inline std::vector<int> getVarsIndex(BDD *bdd);
	inline std::vector<int> getVarsIndex(ADD *bdd);

	inline BDD createMinterm(std::vector<BDD> *x, int node);
	inline BDD createMinterm(std::vector<BDD> *x, std::vector<BDD> *y, unsigned int node_x, unsigned int node_y);
	inline ADD createMinterm(std::vector<ADD> *x, int node);
	inline ADD createMinterm(std::vector<ADD> *x, std::vector<ADD> *y, int x_node, int y_node);

	inline ADD createColumn(std::vector<ADD> *y, unsigned int column);

	inline void createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<BDD> *x, std::vector<BDD> *u, std::vector<BDD> *x_);
	inline void createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<BDD> *x, std::vector<BDD> *x_);
	inline void createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<ADD> *x, std::vector<ADD> *x_);

	inline BDD createXstates(int no_states);

	//!
	inline unsigned int findSequentNode(ADD *APSP_PA, unsigned int *target_node, std::vector<ADD> *x_);

	/* APtoSetSP helper functions */


	inline BDD operatorXUsz(BDD *W, BDD *W_swapped, BDD *Q, BDD *Z, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_);
	inline void relax(BDD *XUz, ADD *APSP_W, BDD *PA_W, ADD *SC, pq_relax *pq_mincost, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_, std::vector<ADD> *add_x, std::vector<ADD> *add_x_);


#ifdef ENABLE_TIME_PROFILING
	long long get_usec(void);
#endif


public:
	//! ShortestPath Constructor.
	/**
	 * It assumes the CUDD manager has already been initialized.
	 * @param mgr_cpp the pointer to the CUDD manager's object.
	 */
	ShortestPath(Cudd *mgr_cpp);

	//! ShortestPath Constructor. Analyzes the System first.
	/**
	 * It assumes the CUDD manager has already been initialized. It takes also as argument the
	 * system in order to analyze it and creates the system variables that are going to be needed
	 * in the methods of this class. The purpose of this is to speed up all methods that need the
	 * system variables. Use this constructor if know in advance that you are going to use more than
	 * one method for one particular system.
	 *
	 * @param mgr_cpp the pointer to the CUDD manager's object.
	 * @param system is the pointer to the BDD of the system.
	 * @param no_states is the number of states the system has.
	 * @param no_inputs is the number of inputs the system has.
	 * @see ShortestPath(Cudd *mgr_cpp)
	 */
	ShortestPath(Cudd *mgr_cpp, BDD *system, unsigned int no_states, unsigned int no_inputs);

	//! ShortestPath De-Constructor.
	/**
	 * Nothing special yet. :)
	 */
	virtual ~ShortestPath();

	//! Create the ADD Target set W given the set of states of W as a vector<int>.
	/**
	 * Given the System's BDD and the cost of each state, as an ADD, this method constructs the
	 * Cost Adjacency Matrix of the System. This is done, by first extracting a valid transition
	 * (x,u,x') and then attaching the corresponding cost of that transition c(x,u,x').
	 * @param system is the pointer to the System's BDD.
	 * @param state_cost is the pointer to the ADD, describing the cost of each state of the system.
	 * @param no_states is the number of states of the System.
	 * @param no_inputs is is the number of the inputs of the System.
	 * @return The ADD of the Cost Adjacency Matrix.
	 */
	ADD createCostAdjacencyMatrix(BDD *system, ADD *state_cost, int no_states = 0, int no_inputs = 0);

	 //! Given the Cost Adjacency Matrix of a DD, get the all-pair shortest path values and the pointer array used to trace back the desired shortest path.
	/**
	 * Method takes as input the Cost Adjacency Matrix of the System as an ADD and returns the all-pairs shortest
	 * path values and the pointer array. To receive the return values, two ADD objects have to be passed as arguments.
	 * Implements the well-known Floyd-Warshall algorithm.\n
	 * __Important__: This method assumes that the arguments (ASPS, PA), which are used to pass the result,
	 * have been already allocated. Empty pointers of these will result to an error!
	 * @param AG is the pointer to the System's BDD.
	 * @param APSP is the _allocated_ pointer to the ADD for returning the APSP cost values.
	 * @param PA is the _allocated_ pointer to the ADD for returning the pointer array of the APSP.
	 */
	void FloydWarshall(ADD *AG, ADD *APSP, ADD *PA);

	//! Finds the shortest path from all pairs to a given target set W. Returns the vector containing the shortest path values and the pointer vector. Supports only deterministic transitions.
	/**
	 * This method takes as input the DDs containing the all-pairs shortest path values, the pointer array of the all-pairs shortest path and a target set W,
	 * for which we want to find the shortest path from all the pairs to that set. The set is given as vector of integers, denoting the states. The method
	 * returns the vector containing the shortest path value as ADD and the pointer vector that shows which node to follow to achieve the shortest path value.\n
	 *
	 * __Important Notice__:
	 * - Memory should have been already allocated for the results (@param APSP_W and @param PA_W).
	 * - The pointer does not contain the "map" to achieve the shortest path value from all pairs not the set. It only contains the intermediate node
	 * to be followed in order to achieve the minimum path. Therefore this result should be used together with the initial pointer array (@param PA).
	 *
	 * @param APSP is the ADD containing the all-pairs shortest path values.
	 * @param PA is the ADD containing the pointer array of the APSP.
	 * @param W	is the desired target set given as a BDD.
	 * @param APSP_W is the returned ADD, containing the all-pairs shortest path values to the set W.
	 * @param PA_W is the returned ADD, containing the intermediate nodes to be followed to achieve the all-pairs shortest path to the set W. This result should be used
	 *        together with the PA ADD.
	 * @see   FloydWarshall, APtoSetSP(ADD *APSP, ADD *PA, std::vector<int> W, ADD *APSP_W, ADD *PA_W)
	 */
	void APtoSetSP(ADD *APSP, ADD *PA, BDD *W, ADD *APSP_W, ADD *PA_W);

	//! Finds the shortest path from all pairs to a given target set W. Returns the vector containing the shortest path values and the pointer vector. Supports also non-deterministic transitions.
	void APtoSetSP(BDD *S, ADD *SC, BDD *W, ADD *APSP_W, BDD *PA_W, unsigned int no_states, unsigned int no_inputs);

	//!
	BDD createControllerBDD(BDD *S, ADD *APSP_PA, ADD *APSP_PA_W);

	//! Dumps the argument BDD to file.
	/** Dumping is done through Dddmp_cuddBddArrayStore. A dummy array of 1 BDD root is used for this purpose.
	 *
	 * @param f is the BDD to be dumped.
	 * @param fname is the file name, containing the dumped BDD.
	 * @param ddname is the DD name (or NULL).
	 * @param varnames is the array of variable names (or NULL)
	 * @param auxids is the array of converted var ids.
	 * @param mode is the storing mode selector.
	 * @param varinfo is used for extra info for the variables in text mode.
	 * @param fp is the file pointer to the store file
	 * @return returns whether the dump was successful or not.
	 * @see dddmp.h in the CUDD Library for more info.
	 */
	bool Dddmp_cuddStore(BDD *f, char *fname, char *ddname = NULL, char **varnames = NULL, int *auxids = NULL, int mode = DDDMP_MODE_TEXT, Dddmp_VarInfoType varinfo = DDDMP_VARIDS, FILE *fp = NULL);

	//! Dumps the argument ADD to file.
	/** Dumping is done through Dddmp_cuddBddArrayStore. A dummy array of 1 BDD root is used for this purpose.
	 *
	 * @param f is the ADD to be dumped.
	 * @param fname is the file name, containing the dumped BDD.
	 * @param ddname is the DD name (or NULL).
	 * @param varnames is the array of variable names (or NULL)
	 * @param auxids is the array of converted var ids.
	 * @param mode is the storing mode selector.
	 * @param varinfo is used for extra info for the variables in text mode.
	 * @param fp is the file pointer to the store file
	 * @return returns whether the dump was successful or not.
	 * @see dddmp.h in the CUDD Library for more info.
	 */
	bool Dddmp_cuddStore(ADD *f, char *fname, char *ddname = NULL, char **varnames = NULL, int *auxids = NULL, int mode = DDDMP_MODE_TEXT, Dddmp_VarInfoType varinfo = DDDMP_VARIDS, FILE *fp = NULL);

	bool checkControllerDom(BDD *contrl, BDD *dom);

};

#endif /* SHORTESTPATH_H_ */
