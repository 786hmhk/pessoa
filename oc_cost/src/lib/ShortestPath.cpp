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

#include "ShortestPath.hh"


//! ShortestPath Constructor.
/**
 * It assumes the CUDD manager has already been initialized.
 * @param mgr_cpp the pointer to the CUDD manager's object.
 */
ShortestPath::ShortestPath(Cudd *mgr_cpp) {

	// Cudd Managers
	this->mgr_cpp = mgr_cpp;
	this->mgr     = mgr_cpp->getManager();

	// System not analyzed this way.
	this->system_analyzed = false;

	this->no_states       = 0;
	this->no_inputs       = 0;
	this->no_state_vars   = 0;
	this->no_input_vars   = 0;

	system_bdd = NULL;

} /* ShortestPath */

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
ShortestPath::ShortestPath(Cudd *mgr_cpp, BDD *system, unsigned int no_states, unsigned int no_inputs){

	// Cudd Managers
	this->mgr_cpp = mgr_cpp;
	this->mgr     = mgr_cpp->getManager();

	// Get the number of states and inputs;
	this->no_states = no_states;
	this->no_inputs = no_inputs;

	// Get the number of variables
	this->no_state_vars = getNoBits(no_states - 1);
	this->no_input_vars = getNoBits(no_inputs - 1);

	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(system);

	if (vars_index.size() != (2 * no_state_vars + no_input_vars)){
		printf("***Critical error! Number of states/inputs mismatch! Revise! (%d:%d)\n", vars_index.size(), (2 * no_state_vars + no_input_vars));
		exit(-1);
	}

	// Create the variables of the System.
	createVariables(vars_index, no_state_vars, no_input_vars, &bdd_x, &bdd_u, &bdd_x_);
	createVariables(vars_index, no_state_vars, no_input_vars, &add_x, &add_x_);

	// System is being analyzed this way.
	this->system_analyzed = true;

	// Get the System's BDD.
	system_bdd = system;


	printf("Self test:\n");
	printf("No States Vars: %3d - No Inputs Vars: %3d\n", this->no_state_vars, this->no_input_vars);
	printf("x  index begin: %3d - end: %3d\n", bdd_x[0].getNode()->index,  bdd_x[bdd_x.size()-1].getNode()->index);
	printf("u  index begin: %3d - end: %3d\n", bdd_u[0].getNode()->index,  bdd_u[bdd_u.size()-1].getNode()->index);
	printf("x' index begin: %3d - end: %3d\n", bdd_x_[0].getNode()->index, bdd_x_[bdd_x_.size()-1].getNode()->index);

} /* ShortestPath */


//! ShortestPath De-Constructor.
/**
 * Nothing special yet. :)
 */
ShortestPath::~ShortestPath() {

} /* ~ShortestPath */


//! Create the Cost Adjacency Matrix of a given System. (Deterministic System)
/**
 * Given the System's BDD and the cost of each state, as an ADD, this method constructs the
 * Cost Adjacency Matrix of the System. This is done, by first extracting a valid transition
 * (x,u,x') and then attaching the corresponding cost of that transition c(x,u,x').\n
 *
 * __Important Notice__: This method assumes only deterministic systems. If for one input, more than
 * one end nodes exist, then only the first one encountered is being taken into account.
 *
 * @param system is the pointer to the System's BDD.
 * @param state_cost is the pointer to the ADD, describing the cost of each state of the system.
 * @param no_states is the number of states of the System.
 * @param no_inputs is is the number of the inputs of the System.
 * @return The ADD of the Cost Adjacency Matrix.
 */
ADD ShortestPath::createCostAdjacencyMatrix(BDD *system, ADD *state_cost, int no_states, int no_inputs){

	ADD AG;

	// BDD system variables
	std::vector<BDD> *bdd_x;
	std::vector<BDD> *bdd_u;
	std::vector<BDD> *bdd_x_;
	// ADD system variables
	std::vector<ADD> *add_x;
	std::vector<ADD> *add_x_;

	ADD add_minterm;
	BDD bdd_minterm;

#ifdef ENABLE_TIME_PROFILING
#ifdef C_CPP_PRINT
	long long start_time = get_usec();
#endif
#endif

	if (!system_analyzed){

		// Allocate memory for the vectors.
		bdd_x  = new std::vector<BDD>;
		bdd_u  = new std::vector<BDD>;
		bdd_x_ = new std::vector<BDD>;
		add_x  = new std::vector<ADD>;
		add_x_ = new std::vector<ADD>;

		// Get the number of variables
		int no_state_vars = getNoBits(no_states - 1);
		int no_input_vars = getNoBits(no_inputs - 1);

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(system);

		// Create the variables of the System.
		createVariables(vars_index, no_state_vars, no_input_vars, bdd_x, bdd_u, bdd_x_);
		createVariables(vars_index, no_state_vars, no_input_vars, add_x, add_x_);
	}
	else{
		bdd_x  = &this->bdd_x;
		bdd_u  = &this->bdd_u;
		bdd_x_ = &this->bdd_x_;

		add_x  = &this->add_x;
		add_x_ = &this->add_x_;
	}

#ifdef C_CPP_PRINT
	printf("ShortestPath::getCostAdjacencyMatrix: No states: %d - No inputs: %d \n", no_states, no_inputs);
#endif


	ADD SC_swapped = state_cost->SwapVariables(*add_x, *add_x_);

	// Get rid of the input u. Keep only x,x'.
	int existental[mgr_cpp->ReadSize()];
	for (unsigned k = 0; k < (*bdd_u)[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = (*bdd_u)[0].getNode()->index; k < (*bdd_x_)[0].getNode()->index; k++){
		existental[k] = 1;
	}
	for (unsigned k = (*bdd_x_)[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 2;
	}

	DdNode *cube_array = Cudd_CubeArrayToBdd(mgr,existental);
	BDD S = system->ExistAbstract(BDD(*mgr_cpp, cube_array), 0);

	// Create the Cost Adjacency Matrix.
	AG = ((mgr_cpp->background() * ((~S).Add())) + mgr_cpp->addOne()) * SC_swapped;
//	AG = system->Add() * SC_swapped;

	/* Iterate over all states. Create the costs of the self-transitions. */
	for (int i = 0; i < no_states; i++){
		/* Add the cost of each state to zero where applicable. */ //TODO: This is if we count the cost of the self transition...
		// Create the minterm.
		bdd_minterm = createMinterm(bdd_x, bdd_x_, i, i);

		if (S.Restrict(bdd_minterm).IsZero()){
			continue;
		}

//		printf("Adding self-cost for state: %d\n", i);

		add_minterm = createMinterm(add_x, add_x_, i, i);
		// Create the constant node. Zero node
		AG = add_minterm.Ite(mgr_cpp->constant(0), AG);
	}

#ifdef ENABLE_TIME_PROFILING
#ifdef C_CPP_PRINT
	/* Print execution time. */
	printf("ShortestPath::getCostAdjacencyMatrix: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
#endif
	return AG;
}/* createCostAdjacencyMatrix */

//! Creates a minterm as an BDD based on the @param node and the number of x variables available.
/**
 * Takes as input the initialized boolean variables and the minterm expressed by an integer and
 * creates the minterm that is expressed by the integer as an BDD.\n
 *
 * __Important Notice__: The maximum value of @param node is 2^(number of x variables). If @param node
 * is larger than this value, then the bits after the max value are not taken into account.
 *
 * @param x has the created boolean variables.
 * @param node is the desired minterm expressed by an integer.
 * @return the created minterm as BDD.
 */
inline BDD ShortestPath::createMinterm(std::vector<BDD> *x, int node){

	unsigned int j;

	BDD one  = mgr_cpp->bddOne();
	BDD zero = mgr_cpp->bddZero();

	BDD minterm, temp;


	minterm = one;

	for (j = 0; j < (*x).size(); j++){
		if (node & 0x01){
			temp = (*x)[j].Ite(minterm, zero);
		}
		else{
			// row
			temp = (*x)[j].Ite(zero, minterm);
		}
		minterm = temp;
		node >>= 1;
	}

	return minterm;
} /* createMinterm */

//!
inline BDD ShortestPath::createMinterm(std::vector<BDD> *x, std::vector<BDD> *y, unsigned int node_x, unsigned int node_y){
	unsigned int j;

	BDD one  = mgr_cpp->bddOne();
	BDD zero = mgr_cpp->bddZero();

	BDD minterm, temp;


	minterm = one;

	for (j = 0; j < (*x).size(); j++){
		if (node_x & 0x01){
			temp = (*x)[j].Ite(minterm, zero);
		}
		else{
			// row
			temp = (*x)[j].Ite(zero, minterm);
		}

		if (node_y & 0x01){
			temp = (*y)[j].Ite(temp, zero);
		}
		else{
			// row
			temp = (*y)[j].Ite(zero, temp);
		}

		minterm = temp;
		node_x >>= 1;
		node_y >>= 1;
	}

	return minterm;
} /* createMinterm */

//! Creates a minterm as an ADD based on the @param node and the number of x variables available.
/**
 * Takes as input the initialized boolean variables and the minterm expressed by an integer and
 * creates the minterm that is expressed by the integer as an ADD.\n
 *
 * __Important Notice__: The maximum value of @param node is 2^(number of x variables). If @param node
 * is larger than this value, then the bits after the max value are not taken into account.
 *
 * @param x has the created boolean variables.
 * @param node is the desired minterm expressed by an integer.
 * @return the created minterm as BDD.
 */
inline ADD ShortestPath::createMinterm(std::vector<ADD> *x, int node){

	unsigned int j;

	ADD one  = mgr_cpp->addOne();
	ADD minterm, temp;

	minterm = one;

	for (j = 0; j < (*x).size(); j++){
		if (node & 0x01){
			temp = minterm * (*x)[j];
		}
		else{
			// row
			temp = minterm * (~(*x)[j]);
		}
		minterm = temp;
		node >>= 1;
	}

	return minterm;
} /* createMinterm */


//! Creates a minterm as an ADD based on the @param node and the number of x variables available.
/**
 * Takes as input the initialized boolean variables x, y and the two sub-minterms expressed
 * by an integer and creates the minterm, which is the logical AND between the two sub-minterms.
 * The resulted minterm is expressed as an ADD. \n
 *
 * __Important Notice__: The maximum value of @param x_node is 2^(number of x variables). If
 * @param x_node is larger than this value, then the bits after the max value are not taken into account.
 * The same applies also for @param y_node.
 *
 * @param x has the created x boolean variables.
 * @param y has the created y boolean variables.
 * @param node is the desired minterm expressed by an integer.
 * @return the created minterm as BDD.
 */
inline ADD ShortestPath::createMinterm(std::vector<ADD> *x, std::vector<ADD> *y, int x_node, int y_node){
	unsigned int j;

	ADD one  = mgr_cpp->addOne();

	ADD minterm, temp;

	minterm = one;

	for (j = 0; j < (*x).size(); j++){
		if (x_node & 0x01){
			temp = minterm * (*x)[j];
		}
		else{
			// row
			temp = minterm * (~(*x)[j]);
		}
		minterm = temp;
		x_node >>= 1;
	}
	for (j = 0; j < (*y).size(); j++){
		if (y_node & 0x01){
			temp = minterm * (*y)[j];
		}
		else{
			// row
			temp = minterm * (~(*y)[j]);
		}
		minterm = temp;
		y_node >>= 1;
	}

	return minterm;
} /* createMinterm */

//! Creates all the variables x, u and x_ given the corresponding indexes.
/**
 * This method works under the assumption that no reordering has been applied to the ADD/BDD to which the indexes correspond.
 * Furthermore it assumes the following order: x-->u-->x_. So:\n
 *
 * - the index for @param x is from 0 to @param no_state_vars - 1
 * - the index for @param u is from @param no_state_vars to @param no_state_vars + @param no_input_vars -1
 * - the index for @param x_ is from @param no_state_vars + @param no_input_vars to 2 *@param no_state_vars + @param no_input_vars -1 \n
 *
 * The result is returned using the pointers to the vectors @param x, @param u and @param x_.
 *
 * __Imporatant Notice__: The memory for the vectors @param x, @param u and @param x_ is assumed to be already allocated.
 *
 * @param vars_index contains all available indexes.
 * @param no_state_vars is the number of state variables (or x) variables being used.
 * @param no_input_vars is the number of input variables (or u) variables being used.
 * @param x is used to return the created x variables.
 * @param u is used to return the created u variables.
 * @param x_ is used to return the created x_ variables.
 * @see createVariables
 */
void ShortestPath::createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<BDD> *x, std::vector<BDD> *u, std::vector<BDD> *x_){

	int i;

	// Create the Variables.
	for (i = 0; i < no_state_vars; i++){
		(*x).push_back(mgr_cpp->bddVar(vars_index[i]));
		(*x_).push_back(mgr_cpp->bddVar(vars_index[i + (no_state_vars + no_input_vars)]));
	}

	for (i = 0; i <  no_input_vars; i++){
		(*u).push_back(mgr_cpp->bddVar(vars_index[no_state_vars + i]));
	}
} /* createVariables */


//! Creates all the variables x, u and x_ given the corresponding indexes.
/**
 * This method works under the assumption that no reordering has been applied to the ADD/BDD to which the indexes correspond.
 * Furthermore it assumes the following order: x-->u-->x_. So:\n
 *
 * - the index for @param x is from 0 to @param no_state_vars - 1
 * - the index for @param x_ is from @param no_state_vars + @param no_input_vars to 2 *@param no_state_vars + @param no_input_vars -1 \n
 *
 * The result is returned using the pointers to the vectors @param x, @param u and @param x_.
 *
 * __Imporatant Notice__: The memory for the vectors @param x, @param u and @param x_ is assumed to be already allocated.
 *
 * @param vars_index contains all available indexes.
 * @param no_state_vars is the number of state variables (or x) variables being used.
 * @param no_input_vars is the number of input variables (or u) variables being used.
 * @param x is used to return the created x variables as an BDD.
 * @param x_ is used to return the created x_ variables as an BDD.
 * @see createVariables
 */
inline void ShortestPath::createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<BDD> *x, std::vector<BDD> *x_){
	int i;

	// Create the Variables.
	for (i = 0; i < no_state_vars; i++){
		(*x).push_back(mgr_cpp->bddVar(vars_index[i]));
		(*x_).push_back(mgr_cpp->bddVar(vars_index[i + (no_state_vars + no_input_vars)]));
	}
}

//! Creates all the variables x and x_ given the corresponding indexes.
/**
 * This method works under the assumption that no reordering has been applied to the ADD/BDD to which the indexes correspond.
 * Furthermore it assumes the following order: x-->u-->x_. So:\n
 *
 * - the index for @param x is from 0 to @param no_state_vars - 1
 * - the index for @param x_ is from @param no_state_vars + @param no_input_vars to 2 *@param no_state_vars + @param no_input_vars -1 \n
 *
 * The result is returned using the pointers to the vectors @param x, @param u and @param x_.
 *
 * __Imporatant Notice__: The memory for the vectors @param x, @param u and @param x_ is assumed to be already allocated.
 *
 * @param vars_index contains all available indexes.
 * @param no_state_vars is the number of state variables (or x) variables being used.
 * @param no_input_vars is the number of input variables (or u) variables being used.
 * @param x is used to return the created x variables as an ADD.
 * @param x_ is used to return the created x_ variables as an ADD.
 * @see createVariables
 */
void ShortestPath::createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<ADD> *x, std::vector<ADD> *x_){

	int i;

	ADD one  = mgr_cpp->addOne();
	ADD zero = mgr_cpp->addZero();

	ADD tempX, tempX_;

	// Create the Variables.
	for (i = 0; i < no_state_vars; i++){

		if (x != NULL){
			tempX  = mgr_cpp->addVar(vars_index[i]);
			tempX  = tempX.Ite(one, zero);
			(*x).push_back(tempX);
		}

		if (x_ != NULL){
			tempX_ = mgr_cpp->addVar(vars_index[i + (no_state_vars + no_input_vars)]);
			tempX_ = tempX_.Ite(one, zero);
			(*x_).push_back(tempX_);
		}
	}
} /* createVariables */


//! Gets all the indexes of the variables being used, given an BDD. The total number of indexes is the total number of variables being used.
inline std::vector<int> ShortestPath::getVarsIndex(BDD *bdd){

	// Vector containing all possible indexes of the BDD's variables.
	std::vector<int> vars_index;

	// Get the BDD.
	std::vector<BDD> system;
	system.push_back(*bdd);

	// Take the union of the supports of each output function.
	BDD support = mgr_cpp->VectorSupport(system);

	DdNode *scan = support.getNode();

	// Get the indexes of the variables from the BDD.
	while (!cuddIsConstant(scan)){
		vars_index.push_back(scan->index);
		scan = cuddT(scan);
	}

	Cudd_RecursiveDeref(mgr_cpp->getManager(),scan);

//	printf("ShortestPath::getVarsIndex: Getting the indexes of the variables (Total: %d).\n", vars_index.size());

//	for (std::vector<int>::iterator i = vars_index.begin(); 	i != vars_index.end(); ++i) {
//		printf("%d - ", *i);
//	}
//	printf("\n");
//	for(i = 0; i < vars_index.size(); i++){
//		printf("%d - ", vars_index[i]);
//	}

	return vars_index;
} /* getVarsIndex */


//! Gets all the indexes of the variables being used, given an ADD. The total number of indexes is the total number of variables being used.
inline std::vector<int> ShortestPath::getVarsIndex(ADD *add){

	// Vector containing all possible indexes of the ADD's variables.
	std::vector<int> vars_index;

	// Get the ADD.
	std::vector<ADD> system;
	system.push_back(*add);

	// Take the union of the supports of each output function.
	BDD support = mgr_cpp->VectorSupport(system);

	DdNode *scan = support.getNode();

	// Get the indexes of the variables from the ADD.
	while (!cuddIsConstant(scan)){
		vars_index.push_back(scan->index);
		scan = cuddT(scan);
	}

	Cudd_RecursiveDeref(mgr_cpp->getManager(),scan);

	return vars_index;
} /* getVarsIndex */

//! Finds the shortest path from all pairs to a given target set W. Returns the vector containing the shortest path values and the pointer vector.
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
void ShortestPath::APtoSetSP(ADD *APSP, ADD *PA, BDD *W, ADD *APSP_W, ADD *PA_W){

	printf("ShortestPath::APtoSetSP\n");

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif
	int i;
	// The x and y variables of the system.
	std::vector<ADD> *add_x;
	std::vector<ADD> *add_y;
	// BDD system variables
	std::vector<BDD> *bdd_x;
	std::vector<BDD> *bdd_y;
	// The total number of x and y variables. (or x and x').
	int no_states;


	if (!system_analyzed){

		// Allocate Memory for the vectors.
		add_x = new std::vector<ADD>;
		add_y = new std::vector<ADD>;
		bdd_x = new std::vector<BDD>;
		bdd_y = new std::vector<BDD>;

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(APSP);

		// Create the variables of the System.
		createVariables(vars_index, vars_index.size()/2, 0, add_x, add_y);
		createVariables(vars_index, vars_index.size()/2, 0, bdd_x, bdd_y);

		// Number of x and y variables.
		no_states = vars_index.size();	// not optimal since it is in the power of 2.
										// Example. If no states 4 then state vars 3.
										// which means no_states = 2*3 = 6 => 2 extra states!!
	}
	else{
		// x and y variables.
		add_x = &this->add_x;
		add_y = &this->add_x_;

		bdd_x = &this->bdd_x;
		bdd_y = &this->bdd_x_;

		// Number of x and y variables.
		no_states = this->no_states;  // this is optimal
	}

	/* Find the shortest path and the pointer array, given the target set W. */
	BDD system_rstct_x, bdd_minterm;
	ADD add_minterm;
	ADD op1 = mgr_cpp->plusInfinity();
	ADD P1  = mgr_cpp->addZero();
	ADD op2;
	ADD P2;
	ADD temp;
	DdNode *result[2];

	int no_states_w = 0;

	/* Iterate over all y (=x') states. */
	// TODO: Check ddPrintMintermAux() in cuddUtil.c
	// it might a better/faster implementation.
	for (i = 0; i < no_states; i++){

//		printf("Checking: %4d ", i);

		// not a valid input.
		if (W->Restrict(createMinterm(bdd_x,i)).IsZero()){
//			printf(" ...not found!\n");
			continue;
		}

		// valid state
//		printf("found target!(%d)\n", i);
		add_minterm = createMinterm(add_y, i);

		op2 = APSP->Cofactor(add_minterm);
		P2  = PA->Cofactor(add_minterm);

		// A little trick to get rid of the Zero in the ADD. :) Instead of zero, now I am getting the node number (+1).
		temp = (P2 * mgr_cpp->minusInfinity()) + mgr_cpp->constant(i+1);
		temp = temp.Maximum(P2);
		//
		Cudd_addApplyMin2(op1.getNode(),op2.getNode(),P1.getNode(),temp.getNode(),result);
		op1 = ADD(*mgr_cpp, result[0]);
		P1  = ADD(*mgr_cpp, result[1]);

		no_states_w++;
	}

	/* Get the result */
	*APSP_W = op1;
	*PA_W   = P1;

	// De-Allocate memory for the vectors.
	if (!system_analyzed){
		delete bdd_x;
		delete bdd_y;
		delete add_x;
		delete add_y;
	}

	printf("Total States in the Target Set: %d - Total States: %d\n", no_states_w, no_states);

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::APtoSetSP: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

} /* APtoSetSP */

//! Creates a BDD with only the x states. I believe this is faster than using Cudd_bddExistAbstract(). //TODO: Check if true.
inline BDD ShortestPath::createXstates(int no_states){
	//
	BDD minterm;
	BDD one  = mgr_cpp->bddOne();
	BDD zero = mgr_cpp->bddZero();
	BDD cofactor;
	BDD temp;

	int i;
	unsigned int j;
	int node;

	/* Create the target set */
	cofactor = zero;

	for (i = 0; i < no_states; i++){

		minterm = one;
		node    = i;

		for (j = 0; j < no_state_vars; j++){
			if (node & 0x01){
				temp = bdd_x[j].Ite(minterm, zero);
			}
			else{
				// row
				temp = bdd_x[j].Ite(zero, minterm);
			}
			minterm = temp;
			node >>= 1;
		}

		// Create the constant node.
		temp     = minterm.Ite(one, cofactor);
		cofactor = temp;
	}

	return cofactor;
}


//!
inline void ShortestPath::relax(BDD *XUz, ADD *APSP_W, BDD *PA_W, ADD *SC, pq_relax *pq_mincost, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_, std::vector<ADD> *add_x, std::vector<ADD> *add_x_){
	printf("ShortestPath::relax\n");

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif


	//
	BDD bdd_minterm_x, bdd_minterm_u;
	BDD XUz_rstct_x, XUz_rstct_u, XUz_rstct_x_;
	BDD PA_W_temp;
	BDD bdd_zero;
	ADD APSP_W_rstct, SC_rstct;
	ADD APSP_W_wo_state;
	ADD SC_rstct_x;
	ADD end_states_x;
	ADD add_minterm_x;

	unsigned int state_sp = 0xFFFF;
	std::vector<unsigned int> input_sp;
	double current_lowest_sp;

	double no_transitions, no_endStates;

	bool found_sp = false;

	bdd_zero = mgr_cpp->bddZero();

//	// Create .dot file
//	std::vector<BDD> nodes_bdd;
//	std::vector<ADD> nodes_add;
//	FILE *outfile;


	/* Iterate over all states. */
	for (unsigned int i = 0; i < no_states; i++){

		// reset current_lowest_sp
		current_lowest_sp = 0xFFFF;

		XUz_rstct_x = XUz->Restrict(createMinterm(bdd_x, i));

		if (XUz_rstct_x.IsZero()){
			continue;
		}

		// If no_transitions == 1, then we have definitely a deterministic transition,
		// else we have to check the inputs. This may speed things up. // TODO: check if it optimizes things.
		no_transitions = XUz_rstct_x.CountMinterm(no_state_vars+no_input_vars);
		printf("-State: %d. No transitions: %f\n", i, no_transitions);

		/* Get Cost Dw(x). */
		APSP_W_rstct = APSP_W->Restrict(createMinterm(add_x, i));

		add_minterm_x = createMinterm(add_x, i);
		bdd_minterm_x = createMinterm(bdd_x, i);

		// Check the inputs
		for (unsigned int j = 0; j < no_inputs; j++){

			bdd_minterm_u = createMinterm(bdd_u, j);
			XUz_rstct_u   = XUz_rstct_x.Restrict(bdd_minterm_u);

			if (XUz_rstct_u.IsZero()){
				continue;
			}


			/* Get Dw(x') + c(x,u,x'). Covers also non-determinism. */
			end_states_x = XUz_rstct_u.SwapVariables(*bdd_x, *bdd_x_).Add();
//			SC_rstct_x = ((*SC) & end_states_x) + ((*APSP_W) & end_states_x);
			SC_rstct_x = SC->Restrict(end_states_x) + APSP_W->Restrict(end_states_x);

			// In case it is non-deterministic: max(Dw(x') + c(x,u,x')).
			SC_rstct_x = SC_rstct_x.FindMax();
			no_endStates = XUz_rstct_u.CountMinterm(no_state_vars);
			printf(" Input: %d. No end states: %f\n", j, no_endStates);

			// Found a shorter path...update, i.e. relax the current state.
			if (APSP_W_rstct > SC_rstct_x){
				printf("  Found shorter path!\n");
				found_sp = true;
				// Update Dw(x)
				APSP_W_wo_state = (*APSP_W) & (~add_minterm_x);
				*APSP_W = APSP_W_wo_state + (add_minterm_x * SC_rstct_x);


				// Shorter path found!
				if (SC_rstct_x.getNode()->type.value < current_lowest_sp){
					current_lowest_sp = SC_rstct_x.getNode()->type.value;
					state_sp = i;
//					input_sp.clear();
//					input_sp.push_back(j);

					printf("Adding state %d with cost %f to the priority queue. Number of inputs: %d\n", state_sp, current_lowest_sp, input_sp.size());
					// Update Priority Queue.
					pq_mincost->push(std::make_pair(current_lowest_sp, state_sp));

					// Creating temp transition for the PA_W
					PA_W_temp = bdd_minterm_x.Ite(bdd_minterm_u, bdd_zero);

				}
				else if (SC_rstct_x.getNode()->type.value == current_lowest_sp){
					// Several inputs can yield the same shortest path.
					if (state_sp == i){
//						input_sp.push_back(j);
						PA_W_temp += bdd_minterm_x.Ite(bdd_minterm_u, bdd_zero);
					}
				}

			}

			// If the number of transitions is one, then we have definitely a deterministic transition.
			if (no_transitions == 1.0){
				break;
			}
		}

		if (found_sp){
			// Update Pw(x). Remember that multiple inputs can yield the same shortest path (sp).
//			PA_W_temp = mgr_cpp->bddZero();
//			for (m = 0; m < input_sp.size(); i++){
//				PA_W_temp  += bdd_minterm_x.Ite(createMinterm(bdd_u, input_sp[i]), bdd_zero);
//			}
			// Remove older pointer entry
			*PA_W &= (~bdd_minterm_x);
			// Add new entry
			*PA_W += PA_W_temp;

		}

	}


//	if (found_sp){
//		// Update Pw(x). Remember that multiple inputs can yield the same shortest path (sp).
//		bdd_minterm_x = createMinterm(bdd_x, state_sp);
//		for (i = 0; i < input_sp.size(); i++){
//			*PA_W  += bdd_minterm_x.Ite(createMinterm(bdd_u, input_sp[i]), bdd_zero);
//		}
//	}
//	else{
//		printf("No shorter path found in this iteration!\n");
//	}



//	nodes_add.push_back(*APSP_W);
//	outfile = fopen("APSP_W_new.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	// Create .dot file
//	nodes_bdd.push_back(*PA_W);
//	outfile = fopen("PA_W_new.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();





#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::relax (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
} /* relax */



//!
inline BDD ShortestPath::operatorXUsz(BDD *W, BDD *W_swapped, BDD *Q, BDD *Z, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_){
	printf("ShortestPath::operatorXUsz\n");

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	BDD Q_swapped, Z_swapped;
	BDD system_rstct_Z, system_rstct_ZW, system_rstct_ZWu;
	BDD rg_ns_xux_;


	//
	Q_swapped = Q->SwapVariables(*bdd_x, *bdd_x_);
	Z_swapped = Z->SwapVariables(*bdd_x, *bdd_x_);


	// Filter-out the states that do not have a transition to Z.
	system_rstct_Z  = system_bdd->Restrict(Z_swapped);

	// Get rid of x'. Keep only x and u.
	int existental[mgr_cpp->ReadSize()];
	for (unsigned k = 0; k < (*bdd_u)[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = (*bdd_u)[0].getNode()->index; k < (*bdd_x_)[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = (*bdd_x_)[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 1;
	}

	DdNode *cube_array  = Cudd_CubeArrayToBdd(mgr,existental);
	BDD system_rstct_Zxu = system_rstct_Z.ExistAbstract(BDD(*mgr_cpp, cube_array), 0);


	// Filter-out the states that belong to the W set. (They are already in the W set :) )
	system_rstct_ZW = system_rstct_Zxu & (~(*W));
//	system_rstct_ZW = system_rstct_Zxu;
	// With this trick, we are going to find out if the states that we have so far all satisfy
	// the reachability game. Namely if all the non-deterministic transitions end up in the Z set
	// This found by the intersection of the system with the above found states: System \cap system_rstct_ZW
	// So, in this case we are going to expose all inputs and then...
	system_rstct_ZWu = (*system_bdd) & system_rstct_ZW;

	// We are going to see which states dot not satisfy the reachability game.
	// These are the intersection of the system_rstct_ZWu with all the states that are not in the union of Z and W:
	// system_rstct_ZWu /cap (Z /cup W)
//	rg_ns_xux_ = system_rstct_ZWu & (~(Z_swapped + (*W_swapped)));
	rg_ns_xux_ = system_rstct_ZWu & (~Z_swapped);


	// Now we know which state/input does not satisfy the reachability game. We only have to filter-out the x' of the BDD.
	// Get rid of x'. Keep only "bad" x and u.
	BDD rg_ns_x = rg_ns_xux_.ExistAbstract(BDD(*mgr_cpp, cube_array), 0);


	// Now, remove undesired states. Of course states in the Z set have also to be removed.
	// This is basically the BDD with the desired (x,u). :)
	BDD XUsz = system_rstct_ZW & (~rg_ns_x) & (~(*Z));
	// We need also the x' states. So final form: (x,u,x').
	XUsz &= (*system_bdd);


	// Testing - delete
	printf("Valid states - inputs:\n");
	for (unsigned int i = 0; i < no_states; i++){

		BDD XUsz_rstr_x = XUsz.Restrict(createMinterm(bdd_x, i));

		if (XUsz_rstr_x.IsZero()){
			continue;
		}

		for (unsigned j = 0; j < no_inputs; j++){
			if(XUsz_rstr_x.Restrict(createMinterm(bdd_u, j)).IsZero()){
				continue;
			}
			printf("(%d,%d)\n",i,j);
		}
	}


//	// Create .dot file
//	std::vector<BDD> nodes_bdd;
//	FILE *outfile;
//
//	nodes_bdd.push_back(Z_swapped);
//	outfile = fopen("Z_swapped.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//
//	nodes_bdd.push_back(system_rstct_Z);
//	outfile = fopen("system_rstct_Z.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	nodes_bdd.push_back(system_rstct_ZW);
//	outfile = fopen("system_rstct_ZW.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	nodes_bdd.push_back(system_rstct_ZWu);
//	outfile = fopen("system_rstct_ZWu.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	nodes_bdd.push_back(rg_ns_xux_);
//	outfile = fopen("rg_ns_xux_.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	nodes_bdd.push_back(rg_ns_x);
//	outfile = fopen("rg_ns_x.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	nodes_bdd.push_back(XUsz);
//	outfile = fopen("XUsz.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::operatorXUsz (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif


	return XUsz;
}

//!
void ShortestPath::initPA_W(BDD *W, BDD *W_swapped, BDD *PA_W, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_){
	printf("ShortestPath::initPA_W\n");

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	BDD PA_W_temp = W->Ite(*W_swapped, mgr_cpp->bddZero());

	*PA_W = PA_W_temp & (*system_bdd);

//	// Create .dot file
//	std::vector<BDD> nodes_bdd;
//	FILE *outfile;
//	nodes_bdd.push_back(*PA_W);
//	outfile = fopen("PA_W_init.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::initPA_W (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
}


//! Finds the shortest path from all pairs to a given target set W. Returns the vector containing the shortest path values and the pointer vector. Supports also non-deterministic transitions.
void ShortestPath::APtoSetSP(BDD *S, ADD *SC, BDD *W, ADD *APSP_W, BDD *PA_W, unsigned int no_states, unsigned int no_inputs){
	printf("ShortestPath::APtoSetSP (ND)\n");


	// BDD system variables
	std::vector<BDD> *bdd_x;
	std::vector<BDD> *bdd_u;
	std::vector<BDD> *bdd_x_;
	// ADD system variables
	std::vector<ADD> *add_x;
	std::vector<ADD> *add_x_;

	unsigned int left_states;

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	if (!system_analyzed){

		// Allocate memory for the vectors.
		bdd_x  = new std::vector<BDD>;
		bdd_u  = new std::vector<BDD>;
		bdd_x_ = new std::vector<BDD>;
		add_x  = new std::vector<ADD>;
		add_x_ = new std::vector<ADD>;

		// Get the number of variables
		int no_state_vars = getNoBits(no_states - 1);
		int no_input_vars = getNoBits(no_inputs - 1);

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(S);

		// Create the variables of the System.
		createVariables(vars_index, no_state_vars, no_input_vars, bdd_x, bdd_u, bdd_x_);
//		createVariables(vars_index, no_state_vars, no_input_vars, add_x, add_x_);
	}
	else{
		bdd_x  = &this->bdd_x;
		bdd_u  = &this->bdd_u;
		bdd_x_ = &this->bdd_x_;

		add_x  = &this->add_x;
		add_x_ = &this->add_x_;
	}

	// The set containing the states of the system
	BDD X;
	// The target set Z. Contains the resolved states. (States where a solution exists)
	BDD Z;
	// Temporary/helper DD
	BDD Q, x_temp;
	// The set Xz \subseteq X that contains the states that guarantee a transitions to Z.
	BDD XUz;

	BDD W_swapped = W->SwapVariables(*bdd_x, *bdd_x_);


	// Priority queue to store the state number with the minimum cost.
	pq_relax pq_mincost;

	/*
	 * Initiliaze the Algorithm
	 */

	left_states = no_states - 1;
	int iteration = 0;

	X = createXstates(no_states);
	// Initialize the shortest path cost function. (Dw)
	*APSP_W = (~(*W)).Add() * mgr_cpp->background();
	// Initialize the pointer function. (Pw)
	*PA_W   = mgr_cpp->bddZero();
	initPA_W(W, &W_swapped, PA_W, bdd_x, bdd_u, bdd_x_);
	// Q = X
	Q = X;
	// Z = 0
	Z = mgr_cpp->bddZero();

	//
	X = X.SwapVariables(*bdd_x, *bdd_x_);

	// Create the priority queue.
	for (unsigned int i = 0; i < no_states; i++){

		if (W->Restrict(createMinterm(bdd_x,i)).IsZero()){
//			pq_mincost.push(std::make_pair(0xFFFF, i));  // TODO: is this needed? Consider it....
			continue;
		}

		pq_mincost.push(std::make_pair(0.0, i));
	}


	/*
	 * Main Loop. Break when Q != empty.
	 */
	while(left_states){

		iteration++;
		printf("Iteration %d\n ", iteration);



		// Get the state with the minimum cost: x = min{Dw(x) | x \in Q}
		x_temp = createMinterm(bdd_x, pq_mincost.top().second);
		printf("Adding state %d to Z. Cost: %f\n", pq_mincost.top().second, pq_mincost.top().first);
		pq_mincost.pop(); // Remove lowest priority state
		// If it is already resolved... continue. (Priority Queue unpleasant property)
		if (!(Z.Restrict(x_temp).IsZero())){
			continue;
		}

		// Z = Z \cup x
		Z = Z + x_temp;

		// Q = Q\x (set minus)
		Q = Q - x_temp;




		// Testing - delete
		printf("States in the Z set: (-");
		for (unsigned int i = 0; i < no_states; i++){
			if (Z.Restrict(createMinterm(bdd_x, i)).IsZero()){
				continue;
			}
			printf("%d-",i);
		}
		printf(")\n");
		// Testing - delete
		printf("States in the Q set: (-");
		for (unsigned int i = 0; i < no_states; i++){
			if (Q.Restrict(createMinterm(bdd_x, i)).IsZero()){
				continue;
			}
			printf("%d-",i);
		}
		printf(")\n");





		// Operator XUsz. Implements both the Xsz and Usz operators.
		XUz = operatorXUsz(W, &W_swapped, &Q, &Z, bdd_x, bdd_u, bdd_x_);

		// Relax
		relax(&XUz, APSP_W, PA_W, SC, &pq_mincost, bdd_x, bdd_u, bdd_x_, add_x, add_x_);




		// Testing - delete
		printf("Current System (x,u):\n");
		for (unsigned int i = 0; i < no_states; i++){

			BDD PA_W_rstr_x = PA_W->Restrict(createMinterm(bdd_x, i));

			if (PA_W_rstr_x.IsZero()){
				continue;
			}

			for (unsigned j = 0; j < no_inputs; j++){
				if(PA_W_rstr_x.Restrict(createMinterm(bdd_u, j)).IsZero()){
					continue;
				}
				printf("(%d,%d)\n",i,j);
			}
		}



		// Some critical error that should not happen dude!
		if (XUz.IsZero()){
			printf("***WTF?? Critical ERROR!*** (System is not satisfying the reachability game)\n");
			return;
		}

		left_states--;

//		// delete
//		if (left_states == (no_states - 4)){
//			break;
//		}
	}





	// Create .dot file
	std::vector<ADD> nodes_add;
	FILE *outfile;
	nodes_add.push_back(*APSP_W);
	outfile = fopen("System_APSP_W.dot", "w");
	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	nodes_bdd.push_back(*PA_W);
	outfile = fopen("System_PA_W.dot", "w");
	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();






#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::APtoSetSP (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

} /* APtoSetSP */

//!
inline unsigned int ShortestPath::findSequentNode(ADD *APSP_PA, unsigned int *target_node, std::vector<ADD> *x_){

	unsigned int sq_node = (unsigned int)APSP_PA->Restrict(createMinterm(x_, *target_node)).getNode()->type.value;

	if (sq_node){
		*target_node = sq_node - 1;
		findSequentNode(APSP_PA, target_node, x_);
	}
	else{
		return *target_node;
	}

	return *target_node;
}

//!
/**
 *
 * @param S
 * @param APSP_PA
 * @param APSP_PA_W
 * @return
 */
BDD ShortestPath::createControllerBDD(BDD *S, ADD *APSP_PA, ADD *APSP_PA_W){
	printf("ShortestPath::createControllerBDD\n");
#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	int i;
	BDD controller;
	BDD restrct_x, restrct_x_;
	BDD x;

	ADD add_restrct_x;

	// BDD system variables
	std::vector<BDD> *bdd_x;
	std::vector<BDD> *bdd_x_;
	std::vector<ADD> *add_x;
	std::vector<ADD> *add_x_;
	// The total number of x and y variables. (or x and x').
	int no_states;
	unsigned int target_node;
	unsigned int fw_node;


	if (!system_analyzed){

		// Allocate Memory for the vectors.
		bdd_x  = new std::vector<BDD>;
		bdd_x_ = new std::vector<BDD>;
		add_x  = new std::vector<ADD>;
		add_x_ = new std::vector<ADD>;

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(S);

		// Create the variables of the System.
		createVariables(vars_index, vars_index.size()/2, 0, bdd_x, bdd_x_);
		createVariables(vars_index, vars_index.size()/2, 0, add_x, add_x_);

		// Number of x and y variables.
		no_states = vars_index.size();	// not optimal since it is in the power of 2.
										// Example. If no states 4 then state vars 3.
										// which means no_states = 2*3 = 6 => 2 extra states!!
	}
	else{

		bdd_x  = &this->bdd_x;
		bdd_x_ = &this->bdd_x_;
		add_x  = &this->add_x;
		add_x_ = &this->add_x_;

		// Number of x and y variables.
		no_states = this->no_states;  // this is optimal
	}

	controller = mgr_cpp->bddZero();

	// Iterate over all states of the system.
	for (i = 0; i < no_states; i++){

		x = createMinterm(bdd_x, i);

		restrct_x   = S->Restrict(x);

		if (restrct_x.IsZero())
			continue;

		target_node = (unsigned int)APSP_PA_W->Restrict(createMinterm(add_x, i)).getNode()->type.value - 1; // because all add entries are +1.

		// Now check how to get to the target node by looking at the FW pointer array. If it is zero, then we have a direct link and nothing has to be done.
//		printf("Checking...: (%d,%d) -> ", i, target_node);
		add_restrct_x = APSP_PA->Restrict(createMinterm(add_x, i));

		fw_node = findSequentNode(&add_restrct_x, &target_node, add_x_);

		// Get the valid inputs. (if there is one/some).
		restrct_x_ = restrct_x.Restrict(createMinterm(bdd_x_, fw_node));

		if (!restrct_x_.IsZero()){
//			printf("(%d,u,%d)\n", i, fw_node);
			controller += x.Ite(restrct_x_, mgr_cpp->bddZero());

//			// TODO: DELETE THIS.
//			for (unsigned int j =0; j <no_inputs; j++){
//				BDD u = createMinterm(&bdd_u, j);
//				if (!(restrct_x_.Restrict(u).IsZero()))
//					printf("(%d,%d)\n", i, j);
//			}
		}
	}

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::createControllerBDD: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

	// De-Allocate memory for the vectors.
	if (!system_analyzed){
		delete bdd_x;
		delete bdd_x_;
		delete add_x;
		delete add_x_;
	}


	return controller;
}

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
bool ShortestPath::Dddmp_cuddStore(BDD *f, char *fname, char *ddname, char **varnames, int *auxids, int mode, Dddmp_VarInfoType varinfo, FILE *fp){
	bool ok;

	DdNode *fn = f->getNode();

	ok = Dddmp_cuddBddStore(mgr_cpp->getManager(), ddname, fn, varnames, auxids, mode, varinfo, fname, fp);

	if (ok) return true;
	else return false;
}


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
bool ShortestPath::Dddmp_cuddStore(ADD *f, char *fname, char *ddname, char **varnames, int *auxids, int mode, Dddmp_VarInfoType varinfo, FILE *fp){
	bool ok;

	DdNode *fn = f->getNode();

	ok = Dddmp_cuddAddStore(mgr_cpp->getManager(), ddname, fn, varnames, auxids, mode, varinfo, fname, fp);

	if (ok) return true;
	else return false;
}


//! Given the Cost Adjacency Matrix of a DD, get the all-pair shortest path values and the pointer array used to trace back the desired shortest path.
/**
* Method takes as input the Cost Adjacency Matrix of the System as an ADD and returns the all-pairs shortest
* path values and the pointer array. To receive the return values, two ADD objects have to be passed as arguments.
* Implements the well-known Floyd-Warshall algorithm.\n
* __Important Notice__:
* - This method assumes that the arguments (ASPS, PA), which are used to pass the result,
* have been already allocated. Empty pointers of these will result to an error!
* - In the pointer array the value of the node is incremented by one. This is done because we use zero to denote that no intermediate node exist. So, when
* for example node 0 is being added to the pointer array for some minterm, then it will show up as 1.
* @param AG is the pointer to the System's BDD.
* @param APSP is the _allocated_ pointer to the ADD for returning the APSP cost values.
* @param PA is the _allocated_ pointer to the ADD for returning the pointer array of the APSP.
*/
void ShortestPath::FloydWarshall(ADD *AG, ADD *APSP, ADD *PA) {

	printf("ShortestPath::FloydWarshall.\n");

	// C++ to C
	DdNode *AG_ = AG->getNode();

	DdNode *S, *P;
	DdNode *one, *zero, *xminterm, *yminterm, *temp_node;
	DdNode *R, *C;
	DdNode **xx;
	DdNode **yy;

	int i, j;
	int matrix_elements;
	int element;
	int no_vars;

	DdNode *Result[2];
	DdNode *TR, *TR_temp;

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	if (!system_analyzed){
		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(AG);
		no_vars = vars_index.size()/2; // spliting into x and y.

		/* Memory allocation */
		xx = (DdNode **) malloc(sizeof(DdNode *) * no_vars);
		yy = (DdNode **) malloc(sizeof(DdNode *) * no_vars);

		/* Create the ADD variables. */
		// It is important here that they have the same index
		// with the ADD they are going to be used with.
		for (i = 0; i < no_vars; i++) {
			xx[i] = Cudd_addIthVar(mgr, vars_index[i]);
			yy[i] = Cudd_addIthVar(mgr, vars_index[i + no_vars]);
			Cudd_Ref(xx[i]);
			Cudd_Ref(yy[i]);
		}

		/* Iterate over all matrix elements, i.e. nodes of the initial DD. */
		matrix_elements = 1 << no_vars;
	}
	else{

		/* Memory allocation */
		xx = (DdNode **) malloc(sizeof(DdNode *) * no_state_vars);
		yy = (DdNode **) malloc(sizeof(DdNode *) * no_state_vars);

		/* Create the ADD variables. */
		for (i = 0; (unsigned int)i < no_state_vars; i++){
			xx[i] = add_x [i].getNode();
			yy[i] = add_x_[i].getNode();
		}

		/* Iterate over all matrix elements, i.e. nodes of the initial DD. */
		matrix_elements = no_states;
		no_vars         = no_state_vars;
	}

	/* Zero and One (constant) nodes. Used to create the minterms.*/
	one  = Cudd_ReadOne(mgr);
	zero = Cudd_ReadZero(mgr);

	/* "Copy" the AG matrix. */
	S = AG_;
	Cudd_Ref(S);
	/* Initialize the Pointer array. */
	TR = Cudd_addConst(mgr, 0);
	Cudd_Ref(TR);


	for (i = 0; i < matrix_elements; i++){
//		printf("Node (%d)\n", i);

		element = i;
		xminterm = one;
		yminterm = one;
		Cudd_Ref(xminterm);
		Cudd_Ref(yminterm);

		/* Creating the minterms. */
		// LSB -- MSB
		for (j = 0; j < no_vars; j++){
			if (element & 1){
				// row
				temp_node = Cudd_addIte(mgr, xx[j], xminterm, zero);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr,xminterm);
				xminterm = temp_node;
				//column
				temp_node = Cudd_addIte(mgr, yy[j], yminterm, zero);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr,yminterm);
				yminterm = temp_node;
			}
			else{
				// row
				temp_node = Cudd_addIte(mgr, xx[j], zero, xminterm);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr,xminterm);
				xminterm = temp_node;
				//column
				temp_node = Cudd_addIte(mgr, yy[j], zero, yminterm);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr,yminterm);
				yminterm = temp_node;
			}

			element >>= 1;
		}

		/* Co-factor the matrix. */
		// row
		R = Cudd_Cofactor(mgr, S, xminterm);
		Cudd_Ref(R);
		Cudd_RecursiveDeref(mgr, xminterm);
		// column
		C = Cudd_Cofactor(mgr, S, yminterm);
		Cudd_Ref(C);
		Cudd_RecursiveDeref(mgr, yminterm);



		/* Compute the outer sum. */
//		P = Cudd_addOuterSum(mgr,S,R,C); // if you want only the APSP and not the PA.
		AddOuterSumTrace(S,R,C,Result,i+1);
		P = Result[0];
		Cudd_Ref(P);
		Cudd_Ref(Result[1]);
		// Caution here! It might happen that in previous steps the c(i,k) + c(k,j)
		// to have been already less than c(i,j), meaning that we have already recorded
		// a pointer in the pointer array for that path. Now we find out that a new
		// c(i,k') + c(k',j) is even less that the previous c(i,k) + c(k,j), so the pointer,
		// has to be updated. Now, since we are iterating over all nodes, starting from 0,
		// we know that the new updated pointer has to be larger than the old one. This is
		// why we use maximum here.
		TR_temp = Cudd_addApply(mgr,Cudd_addMaximum,TR,Result[1]);
		Cudd_Ref(TR_temp);

		/* Keep the new matrix P. Delete others. */
		Cudd_RecursiveDeref(mgr,R);
		Cudd_RecursiveDeref(mgr,C);
		Cudd_RecursiveDeref(mgr,S);
		Cudd_RecursiveDeref(mgr,TR);
		Cudd_RecursiveDeref(mgr, Result[1]);

		S  = P;
		TR = TR_temp;
	}

	Cudd_Deref(S);
	Cudd_Deref(TR);

	/* Memory De-allocation. */
	free(xx);
	free(yy);

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::FloydWarshall: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

	/* Create the C++ objects. */
	*APSP = ADD(*mgr_cpp, S);
	*PA   = ADD(*mgr_cpp, TR);
} /* FloydWarshall */


//! Re-implemented method of the AddOuterSum function of the cudd library, to be used for the @ref FloydWarshall method.
/**
 * Takes the minimum of a matrix and the outer sum of two vectors. It also keeps track the node (@param node), if the
 * outer sum of the two vectors is minimum than the matrix. This outcome is used to construct the pointer array of the
 * all-pair shortest path problem. So, Result[0] contains the minumum and Result[1] contains the minterm (based on the node
 * that is being processed) of the pointer array.
 *
 * __Important__: The memory for the result (@param Result) should be already allocated, before passing it as an argument.
 *
 * @param M
 * @param r
 * @param c
 * @param Result
 * @param node
 * @see FloydWarshall, AddOuterSumRecurTrace
 */
void ShortestPath::AddOuterSumTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node){
    do {
		mgr_cpp->getManager()->reordered = 0;
		AddOuterSumRecurTrace(M, r, c, Result, node);
	} while (mgr_cpp->getManager()->reordered == 1);
} /* AddOuterSumTrace */



//! Re-implemented method of the AddOuterSumRecur function of the cudd library, to be used for the @ref AddOuterSumTrace method.
/**
 * Implements the recursive step of the @ref AddOuterSumTrace.
 */
void ShortestPath::AddOuterSumRecurTrace(DdNode *M, DdNode *r, DdNode *c, DdNode **Result, unsigned int node){

//	printf("cuddAddOuterSumRecurTrace\n");

	DdManager *mgr = mgr_cpp->getManager();
	DdNode *R, *T, *Mt, *Me, *Tt, *Te, *rt, *re, *ct, *ce, *Rt, *Re;
	int topM, topc, topr;
	int v, index;

	statLine(mgr_cpp->getManager());
	/* Check special cases. */
	if (r == DD_PLUS_INFINITY(mgr) || c == DD_PLUS_INFINITY(mgr)) {
//		printf("cuddAddOuterSumRecurTrace: return PLUS INFINITY\n");
		Result[0] = M;
		Result[1] = cuddUniqueConst(mgr, 0);
		return;
	}

	if (cuddIsConstant(c) && cuddIsConstant(r)) {
		R = cuddUniqueConst(mgr, Cudd_V(c) + Cudd_V(r) );
		cuddRef(R);
		if (cuddIsConstant(M)) {
			if (cuddV(R) < cuddV(M)) {
				cuddDeref(R);
				Result[0] = R;
				Result[1] = cuddUniqueConst(mgr, node);
				return;
			} else {
				Cudd_RecursiveDeref(mgr,R);
				Result[0] = M;
				Result[1] = cuddUniqueConst(mgr, 0);
				return;
			}
		}
		else {
			DdNode *min_result[2];
			Cudd_addApplyMinTrace(Cudd_addMinimum, R, M, R, min_result, node);
			cuddRef(min_result[0]);
			cuddRef(min_result[1]);

			Result[0] = min_result[0];
			// Caution!
			// It might be the case where d(i,j) = d(i,k) + d(k,j). So we can add the k node
			// to the pointer matrix, but it is better not to. :)
			if (M == min_result[0]){
				Result[1] = cuddUniqueConst(mgr, 0);
			}
			else{
				Result[1] = min_result[1];
			}
			Cudd_RecursiveDeref(mgr, R);
			Cudd_Deref(min_result[0]);
			Cudd_Deref(min_result[1]);
			return;
		}
	}


#ifdef ENABLE_CACHE
	/* Check the cache. */ // TODO: Need to check the cache. It is causing errors.
//	R = cuddCacheLookup(mgr, DD_ADD_OUT_SUM_TAG, M, r, c);
//	if (R != NULL) {
//		printf("cuddAddOuterSumRecurTrace: return cache (%d)\n", node);
//		Result[0] = R;
//		R = cuddCacheLookup(mgr, DD_ADD_OUT_SUM_TRACE_TAG, M, r, c);
//		if (R != NULL){
//			Result[1] = R;
//		}
//		else{
//			printf("ShortestPath::AddOuterSumRecurTrace. Cache Miss might cause error!");
//			Result[1] = NULL;
//		}
//		return;
//	}
#endif

	// Find the top variable.
	topM = cuddI(mgr,M->index);
	topr = cuddI(mgr,r->index);
	topc = cuddI(mgr,c->index);
	v = ddMin(topM,ddMin(topr,topc));

	/* Compute cofactors. */
	if (topM == v) {
		Mt = cuddT(M);
		Me = cuddE(M);
	} else {
		Mt = Me = M;
	}
	if (topr == v) {
		rt = cuddT(r);
		re = cuddE(r);
	} else {
		rt = re = r;
	}
	if (topc == v) {
		ct = cuddT(c);
		ce = cuddE(c);
	} else {
		ct = ce = c;
	}

	/* Recursively solve. */
	// If
	AddOuterSumRecurTrace(Mt, rt, ct, Result, node);
	Rt = Result[0];
	Tt = Result[1];
	if (Rt == NULL) {
		Result[0] = NULL;
		if (Tt == NULL)
			Result[1] = NULL;
		return;
	}
	cuddRef(Rt);
	cuddRef(Tt);

	// Else
	AddOuterSumRecurTrace(Me, re, ce, Result, node);
	Re = Result[0];
	Te = Result[1];
	if (Re == NULL) {
		Cudd_RecursiveDeref(mgr, Rt);
		Result[0] = NULL;
		if (Te == NULL)
			Result[1] = NULL;
		return;
	}
	cuddRef(Re);
	cuddRef(Te);

	index = mgr->invperm[v];
	R = (Rt == Re) ? Rt : cuddUniqueInter(mgr, index, Rt, Re);
	T = (Tt == Te) ? Tt : cuddUniqueInter(mgr, index, Tt, Te);

	if (R == NULL) {
		Cudd_RecursiveDeref(mgr, Rt);
		Cudd_RecursiveDeref(mgr, Re);
		Result[0] = NULL;
		if (T == NULL) {
			Cudd_RecursiveDeref(mgr, Tt);
			Cudd_RecursiveDeref(mgr, Te);
			Result[1] = NULL;
		}
		return;
	}
	cuddDeref(Rt);
	cuddDeref(Re);
	cuddDeref(Tt);
	cuddDeref(Te);

#ifdef ENABLE_CACHE
	/* Store the result in the cache. */
//	cuddCacheInsert(mgr, DD_ADD_OUT_SUM_TAG, M, r, c, R);
//	cuddCacheInsert(mgr, DD_ADD_OUT_SUM_TRACE_TAG, M, r, c, T);
#endif

	Result[0] = R;
	Result[1] = T;

} /* AddOuterSumRecurTrace */




//! Re-implemented method of the Cudd_addApply function of the cudd library, to be used for the @ref AddOuterSumTrace method.
/**
 * This method applies the min(f,g), atm, but also returns the minterm of the pointer array need in the all-pairs shortest path
 * problem. A minterm of the pointer array is created for each matrix element (i.e. each node of the DD). So, This applies the
 * minimum to the corresponding discriminants of f and g. Based on the outcome of each end-node of the @param f and @param g, the
 * either the current matrix element (@param node) is recorded by @param R or the zero value is used.\n
 *
 * __Important__: The memory for the result (@param Result) should be already allocated, before passing it as an argument.
 *
 * @param op is the desired operation to be applied between @param f and @param g. Currently only the minimum is needed/supported.
 * @param f	is the first operand.
 * @param g is the second operand.
 * @param R the x-@param node minterm of the pointer array.
 * @param Result is the result of the operation, where Result[0] contains the min(f,g) and Result[1] the minterm of the pointer array.
 * @param node is the value of the current matrix element (DD node) to be processed. This is needed to record the current node for the
 *        pointer array.
 */
void ShortestPath::Cudd_addApplyMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node){

    do {
		mgr_cpp->getManager()->reordered = 0;
		cuddAddApplyRecurMinTrace(op,f,g,R,Result,node);
    } while (mgr_cpp->getManager()->reordered == 1);
} /* Cudd_addApplyTrace */


//! Re-implemented method of the Cudd_addApply function of the cudd library, to be used for the @ref Cudd_addApplyMinTrace method.
/**
 * Implements the recursive step of the @ref Cudd_addApplyMinTrace.
 */
void ShortestPath::cuddAddApplyRecurMinTrace(DD_AOP op, DdNode * f, DdNode * g, DdNode * R, DdNode **Result, int node)
{
    DdNode *res,
	   *fv, *fvn, *gv, *gvn,
	   *T0, *T1, *E0, *E1;
    unsigned int ford, gord;
    unsigned int index;


    /* Check terminal cases. Op may swap f and g to increase the
     * cache hit rate.
     */
    statLine(dd);
    res = (*op)(mgr_cpp->getManager(),&f,&g);
    if (res != NULL) {
    	Result[0] = res;
    	if (cuddV(R) == cuddV(res)){
    		Result[1] = cuddUniqueConst(mgr_cpp->getManager(), node);
    	}
    	else{
    		Result[1] = cuddUniqueConst(mgr_cpp->getManager(), 0);
    	}
    	return;
    }

#ifdef ENABLE_CACHE
    /* Check cache. */
    DD_CTFP cacheOp = (DD_CTFP) op;
    res = cuddCacheLookup2(mgr_cpp->getManager(),cacheOp,f,g);
    if (res != NULL) {
    	printf("***ShortestPath::cuddAddApplyRecurTrace: Cache hit!\n");
    	Result[0] = res;
    	Result[1] = cuddCacheLookup(mgr_cpp->getManager(), DD_ADD_MINIMUM_TRACE_TAG, f, g, R); // TODO: check this!
    	return;
    }
#endif

    /* Recursive step. */
	ford = cuddI(mgr_cpp->getManager(),f->index);
	gord = cuddI(mgr_cpp->getManager(),g->index);
	if (ford <= gord) {
		index = f->index;
		fv = cuddT(f);
		fvn = cuddE(f);
	} else {
		index = g->index;
		fv = fvn = f;
	}
	if (gord <= ford) {
		gv = cuddT(g);
		gvn = cuddE(g);
	} else {
		gv = gvn = g;
	}

    cuddAddApplyRecurMinTrace(op,fv,gv,R,Result,node);
    T0 = Result[0];
    T1 = Result[1];
    if (T0 == NULL) {
    	Result[0] = NULL;
    	Result[1] = NULL;
    	return;
    }
    cuddRef(T0);
    cuddRef(T1);

    cuddAddApplyRecurMinTrace(op,fvn,gvn,R,Result,node);
    E0 = Result[0];
    E1 = Result[1];
    if (E0 == NULL) {
    	Result[0] = NULL;
    	Result[1] = NULL;
    	Cudd_RecursiveDeref(mgr_cpp->getManager(),T0);
    	Cudd_RecursiveDeref(mgr_cpp->getManager(),T1);
    	return;
    }
    cuddRef(E0);
    cuddRef(E1);


    res = (T0 == E0) ? T0 : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T0,E0);
    if (res == NULL) {
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T0);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E0);
    	Result[0] = NULL;
    	Result[1] = NULL;
		return;
    }
    cuddDeref(T0);
    cuddDeref(E0);

    Result[0] = res;

    res = (T1 == E1) ? T1 : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T1,E1);
    if (res == NULL) {
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T1);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E1);
    	Result[0] = NULL;
    	Result[1] = NULL;
		return;
    }
    cuddDeref(T1);
    cuddDeref(E1);

    Result[1] = res;

#ifdef ENABLE_CACHE
    /* Store result. */
    cuddCacheInsert2(mgr_cpp->getManager(),cacheOp,f,g,Result[0]);
    cuddCacheInsert(mgr_cpp->getManager(), DD_ADD_MINIMUM_TRACE_TAG, f, g, R, Result[1]);
#endif
} /* cuddAddApplyRecurTrace */


//! ADD Integer and floating point minimum without node swapping. Based on the implementation of the cudd library.
/**
 * Integer and floating point min for Cudd_addApply*. Returns NULL if not a terminal case; min(f,g) otherwise.
 *
 * @param dd is the cudd manager. (cudd library)
 * @param f is the first operand. (cudd library)
 * @param g is the second operand. (cudd library)
 * @return returns the minimum of f,g.
 * @see Cudd_addMinimum, Cudd_addApply
 */
DdNode *ShortestPath::Cudd_addMinimumNS(DdManager * dd, DdNode ** f, DdNode ** g){
	DdNode *F, *G;

	F = *f;
	G = *g;
	if (F == DD_PLUS_INFINITY(dd))
		return (G);
	if (G == DD_PLUS_INFINITY(dd))
		return (F);
	if (F == G){
//		printf("EQUAL!!!!\n");
		return (F);
	}


	if (cuddIsConstant(F) && cuddIsConstant(G)) {
		if (cuddV(F) <= cuddV(G)) {
			return (F);
		} else {
			return (G);
		}
	}
	return (NULL);

} /* end of Cudd_addMinimum */


//! Re-implemented method of the Cudd_addApply function of the cudd library, to be used for the @ref APtoSetSP method.
/**
 * This is implemented for the @ref APtoSetSP method. This Applies the minimum to the corresponding discriminants of f and g. Based on the outcome of
 * each end node of the @param f and @param g, the correct pointer (value) of the two "pointer"-discriminants @param Pf and @param Pg is beeing
 * chosen.\n
 *
 * Returns a pointer to the result if successful; NULL otherwise.\n
 *
 * __Important__: The memory for the Result should be allocated before passing it as an argument.
 *
 * @param f is the first operand.
 * @param g is the second operand.
 * @param Pf is the first discriminants of the pointer array
 * @param Pg is the second discriminants of the pointer array
 * @param Result is the result of min(f,g) in Result[0] and the updated minterm of pointer array in Result[1].
 * @see APtoSetSP, Cudd_addMonadicApply Cudd_addPlus Cudd_addTimes Cudd_addThreshold Cudd_addSetNZ Cudd_addDivide Cudd_addMinus
 *  Cudd_addMinimum Cudd_addMaximum Cudd_addOneZeroMaximum Cudd_addDiff Cudd_addAgreement Cudd_addOr Cudd_addNand Cudd_addNor
 *  Cudd_addXor Cudd_addXnor
 */

void ShortestPath::Cudd_addApplyMin2(DdNode * f, DdNode * g, DdNode * Pf, DdNode * Pg, DdNode **Result){

    do {
		mgr_cpp->getManager()->reordered = 0;
		cuddAddApplyMin2Recur(Cudd_addMinimumNS,f,g,Pf,Pg,Result);
    } while (mgr_cpp->getManager()->reordered == 1);
} /* Cudd_addApplyTrace */

//! Re-implemented method of the Cudd_addApply function of the cudd library, to be used for the @ref Cudd_addApplyMin2 method.
/**
 * Implements the recursive step of the @ref Cudd_addApplyMin2.
 */
void ShortestPath::cuddAddApplyMin2Recur(DD_AOP op, DdNode * f, DdNode * g, DdNode * Pf, DdNode * Pg, DdNode **Result)
{
    DdNode *res,
	   *fv, *fvn, *gv, *gvn,
	   *Pfv, *Pfvn, *Pgv, *Pgvn,
	   *T0, *T1, *E0, *E1;
    unsigned int index;
    unsigned int topf, topg, topPf, topPg, v;
//    DD_CTFP cacheOp;

    /* Check terminal cases. Op may swap f and g to increase the
     * cache hit rate.
     */
    statLine(dd);

    res = (*op)(mgr_cpp->getManager(),&f,&g);
    if (res != NULL) {
//    	printf("Res! \n");
    	Result[0] = res;
    	if(cuddIsConstant(Pf) && cuddIsConstant(Pg)){
//    		printf("Res 1. Is constant! \n");
    		if(cuddV(f) < cuddV(g)){
    			Result[1] = Pf;
    		}
    		else {
    			Result[1] = Pg;
    		}
    		return;
    	}
    }


#ifdef ENABLE_CACHE
    /* Check cache. */	// TODO: Need to check if the Cache is working like that. Until then, it is disabled.
//    cacheOp = (DD_CTFP) op;
//    res = cuddCacheLookup2(mgr_cpp->getManager(),cacheOp,f,g);
//    if (res != NULL) {
//    	Result[0] = res;
//    	res = cuddCacheLookup2(mgr_cpp->getManager(),cacheOp,Pf,Pg);
//    	if (res != NULL){
//	    	printf("***ShortestPath::cuddAddApplyRecurTrace: Double Cache hit!\n");
//    		Result[1] = res;
//    		return;
//    	}
//    }
#endif




	// Find the top variable.
	topf  = cuddI(mgr,f->index);
	topg  = cuddI(mgr,g->index);
	topPf = cuddI(mgr,Pf->index);
	topPg = cuddI(mgr,Pg->index);
	v = ddMin(topf,ddMin(topg,ddMin(topPf,topPg)));

	/* Compute cofactors. */
	if (topf == v) {
		fv  = cuddT(f);
		fvn = cuddE(f);
	} else {
		fv = fvn = f;
	}
	if (topg == v) {
		gv  = cuddT(g);
		gvn = cuddE(g);
	} else {
		gv = gvn = g;
	}
	if (topPf == v) {
		Pfv  = cuddT(Pf);
		Pfvn = cuddE(Pf);
	} else {
		Pfv = Pfvn = Pf;
	}
	if (topPg == v) {
		Pgv  = cuddT(Pg);
		Pgvn = cuddE(Pg);
	} else {
		Pgv = Pgvn = Pg;
	}

	index = mgr->invperm[v];

	// True
	cuddAddApplyMin2Recur(op,fv,gv,Pfv,Pgv,Result);
    T0 = Result[0];
    T1 = Result[1];
    if (T0 == NULL) {
    	Result[0] = NULL;
    	Result[1] = NULL;
    	return;
    }
    cuddRef(T0);
    cuddRef(T1);


    // Else
    cuddAddApplyMin2Recur(op,fvn,gvn,Pfvn,Pgvn,Result);
    E0 = Result[0];
    E1 = Result[1];
    if (E0 == NULL) {
    	Result[0] = NULL;
    	Result[1] = NULL;
    	Cudd_RecursiveDeref(mgr_cpp->getManager(),T0);
    	Cudd_RecursiveDeref(mgr_cpp->getManager(),T1);
    	return;
    }
	cuddRef(E0);
	cuddRef(E1);


	// If it exists return it, otherwise create the node.
	res = (T0 == E0) ? T0 : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T0,E0);
	if (res == NULL) {
		printf("Res0 cuddUniqueInter FAIL!\n");
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T0);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E0);
		Result[0] = NULL;
		Result[1] = NULL;
		return;
	}
	cuddDeref(T0);
	cuddDeref(E0);
	// Result
	Result[0] = res;

	// If it exists return it, otherwise create the node.
	res = (T1 == E1) ? T1 : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T1,E1); // this is normally index1.
	if (res == NULL) {
		printf("Res1 cuddUniqueInter FAIL!\n");
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T0);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E0);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T1);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E1);
		Result[0] = NULL;
		Result[1] = NULL;
		return;
	}
	cuddDeref(T1);
	cuddDeref(E1);
	// Result
	Result[1] = res;


    /* Store result. */
#ifdef ENABLE_CACHE
//    cuddCacheInsert2(mgr_cpp->getManager(),cacheOp,f,g,Result[0]);	//TODO: need to check the cache if it works ok.
//    cuddCacheInsert2(mgr_cpp->getManager(),cacheOp,Pf,Pg,Result[0]);
#endif

} /* cuddAddApplyRecurTrace */


//! Get the number of bits of an integer.
unsigned int ShortestPath::getNoBits(unsigned int number){

	unsigned int no_bits = 0;

	for(;;){
		number >>= 1;
		no_bits++;
		if (number == 0)
			break;
	}
	return no_bits;
}



#ifdef ENABLE_TIME_PROFILING
long long ShortestPath::get_usec(void){
	long long r;
	struct timeval t;


	gettimeofday(&t, NULL);
	r = t.tv_sec * 1000000 + t.tv_usec;
	return r;
}
#endif




bool ShortestPath::checkControllerDom(BDD *contrl, BDD *dom){

	printf("ShortestPath::checkControllerDom\n");
#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	// BDD system variables
	std::vector<BDD> *bdd_x;
	std::vector<BDD> *bdd_y;
	// The total number of x and y variables. (or x and x').
	int no_states;
	unsigned int no_states_w = 0;


	if (!system_analyzed){

		// Allocate Memory for the vectors.
		bdd_x = new std::vector<BDD>;
		bdd_y = new std::vector<BDD>;

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(contrl);

		// Create the variables of the System.
		createVariables(vars_index, vars_index.size()/2, 0, bdd_x, bdd_y);

		// Number of x and y variables.
		no_states = vars_index.size();	// not optimal since it is in the power of 2.
										// Example. If no states 4 then state vars 3.
										// which means no_states = 2*3 = 6 => 2 extra states!!
	}
	else{

		bdd_x = &this->bdd_x;
		bdd_y = &this->bdd_x_;

		// Number of x and y variables.
		no_states = this->no_states;  // this is optimal
	}

	BDD x_restrct;

	for (int i = 0; i < no_states; i++){

		x_restrct = createMinterm(bdd_x, i);

		if((contrl->Restrict(x_restrct).IsZero() && dom->Restrict(x_restrct).IsZero()) || (contrl->Restrict(x_restrct).IsOne() && dom->Restrict(x_restrct).IsOne())){
			no_states_w++;
		}
		else{
			printf("ShortestPath::checkControllerDom: Controller domain does not match the actual Controller domain! State: %d\n", i);
#ifdef ENABLE_TIME_PROFILING
			/* Print execution time. */
			printf("ShortestPath::checkControllerDom: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
			return false;
		}
	}

	printf("Total States in the Controller Set: %d - Total States: %d\n", no_states_w, no_states);

	// De-Allocate memory for the vectors.
	if (!system_analyzed){
		delete bdd_x;
		delete bdd_y;
	}

#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::checkControllerDom: Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

	return true;
}

