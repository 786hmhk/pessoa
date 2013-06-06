/**
 * @file
 * @author Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>
 * @version 0.5
 *
 * @section LICENSE
 *
 * Copyright (c) <2013>, <TU Delft: Delft University of Technology>
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

// C Libs
#include "cudd.h"
#include "util.h"
#include "cuddInt.h"
// C++ Lib
#include "cuddObj.hh"

//! Enables or Disables the use of the cudd cache in some operations.
//#define ENABLE_CACHE

// Cache Tags
//! Cache Tag for the cache used in AddOuterSumRecurTrace method.
#define DD_ADD_OUT_SUM_TRACE_TAG		0x6f
//! Cache Tag for the cache used in cuddAddApplyRecurTrace method.
#define DD_ADD_MINIMUM_TRACE_TAG		0x6d

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

	void AddOuterSumTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node);
	void AddOuterSumRecurTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node);

	void cuddAddMinimumRecur(DdNode ** f, DdNode ** g, DdNode * R, DdNode **Result, int node);
	void Cudd_addApplyMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node);
	void cuddAddApplyRecurMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node);

	void minimum2(ADD *f, ADD *g, ADD *P);
	static DdNode *Cudd_addMinimumNS(DdManager * dd, DdNode ** f, DdNode ** g);
	void Cudd_addApplyMin2(DdNode * f0, DdNode * g0, DdNode * Pf, DdNode * Pg, DdNode **Result);
	void cuddAddApplyMin2Recur(DD_AOP op, DdNode * f0, DdNode * g0, DdNode * Pf, DdNode * Pg, DdNode **Result);

	inline std::vector<int> getVarsIndex(BDD *bdd);
	inline std::vector<int> getVarsIndex(ADD *bdd);

	inline BDD createMinterm(std::vector<BDD> *x, int node);
	inline ADD createMinterm(std::vector<ADD> *x, int node);
	inline ADD createMinterm(std::vector<ADD> *x, std::vector<ADD> *y, int x_node, int y_node);

	inline ADD createColumn(std::vector<ADD> *y, unsigned int column);

	inline void createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<BDD> *x, std::vector<BDD> *u, std::vector<BDD> *x_);
	inline void createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<ADD> *x, std::vector<ADD> *x_);

	//! Create the ADD Target set W given the set of states of W as a vector<int>.
	/**
	 * Method takes as input the Systems ADD or any other ADD that contains the System's variables,
	 * the vector containing the states of the target set and returns the ADD of the target set.
	 * @param system is the pointer to the System's ADD or any other ADD that contains the System's variables.
	 * @param no_states is the number of states of the System.
	 * @param target_set is the vector containing the states of the target set W.
	 * @return The ADD of the target set W.
	 */
	ADD createTargetSet(ADD *system, int no_states, int no_inputs, std::vector<int> target_set);

	//! Create the ADD Target set W given the set of states of W as a vector<int>.
	/**
	 * Method takes as input the Systems BDD, the vector containing the states
	 * of the target set and returns the ADD of the target set.
	 * @param system is the pointer to the System's BDD.
	 * @param no_states is the number of states of the System.
	 * @param target_set is the vector containing the states of the target set W.
	 * @return The ADD of the target set W.
	 */
	ADD createTargetSet(BDD *system, int no_states, int no_inputs, std::vector<int> target_set);




public:
	//! ShortestPath Constructor.
	/**
	 * It assumes the CUDD manager has already been initialized.
	 * @param mgr_cpp the pointer to the CUDD manager's object.
	 */
	ShortestPath(Cudd *mgr_cpp);

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
	ADD createCostAdjacencyMatrix(BDD *system, ADD *state_cost, int no_states, int no_inputs);

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

	//! To be implemented. (if needed).
	void APtoSetSP(ADD *APSP, ADD *PA, ADD *W, ADD *APSP_W, ADD *PA_W);

	//! Finds the shortest path from all pairs to a given target set @param W. Returns the vector containing the shortest path values and the pointer vector.
	/**
	 * This method takes as input the DDs containing the all-pairs shortest path values, the pointer array of the all-pairs shortest path and a target set W,
	 * for which we want to find the shortest path from all the pairs to that set. The set is given as vector of integers, denoting the states. The method
	 * returns the vector containing the shortest path value as ADD and the pointer vector that shows which node to follow to achieve the shortest path value.\n
	 *
	 * __Important__:
	 * - Memory should have been already allocated for the results (@param APSP_W and @param PA_W).
	 * - The pointer does not contain the "map" to achieve the shortest path value from all pairs not the set. It only contains the intermediate node
	 * to be followed in order to achieve the minimum path. Therefore this result should be used together with the initial pointer array (@param PA).
	 *
	 * @param APSP is the ADD containing the all-pairs shortest path values.
	 * @param PA is the ADD containing the pointer array of the APSP.
	 * @param W	is the desired target set. Nodes are expressed as integers.
	 * @param APSP_W is the returned ADD, containing the all-pairs shortest path values to the set W.
	 * @param PA_W is the returned ADD, containing the intermediate nodes to be followed to achieve the all-pairs shortest path to the set W. This result should be used
	 *        together with the PA ADD.
	 */
	void APtoSetSP(ADD *APSP, ADD *PA, std::vector<int> W, ADD *APSP_W, ADD *PA_W);
};

#endif /* SHORTESTPATH_H_ */
