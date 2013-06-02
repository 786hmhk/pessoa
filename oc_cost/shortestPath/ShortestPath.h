/*
 * ShortestPath.h
 *
 *  Created on: May 30, 2013
 *      Author: tanasaki
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


class ShortestPath {

private:

	Cudd *mgr_cpp;

	void
	AddOuterSumTrace(
	  DdNode *M,
	  DdNode *r,
	  DdNode *c,
	  DdNode **Result,
	  unsigned int node);

	void
	AddOuterSumRecurTrace(
	  DdNode *M,
	  DdNode *r,
	  DdNode *c,
	  DdNode **Result,
	  unsigned int node);

	inline std::vector<int> getVarsIndex(BDD *bdd);
	inline std::vector<int> getVarsIndex(ADD *bdd);
	bool createBDDVariables(BDD *system, std::vector<int> vars_index, int no_states, int no_inputs, BDD *x, BDD *u, BDD *x_);
	bool createADDVariables(BDD *system, std::vector<int> vars_index, int no_states, int no_inputs, ADD *x, ADD *x_);
	BDD BDDTransition(BDD *x, BDD *u, BDD *x_, int no_states, int no_inputs, int xi, int ui, int xi_);

public:
	ShortestPath(Cudd *mgr_cpp);
	virtual ~ShortestPath();

	ADD getCostAdjacencyMatrix(BDD *system, ADD *cost, int no_states, int no_inputs);


	bool FloydWarshall(ADD *AG);
};

#endif /* SHORTESTPATH_H_ */
