/**
 * @file
 * @author Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>
 * @version 0.9
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

	optimized = false;

} /* end of ShortestPath */

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
ShortestPath::ShortestPath(Cudd *mgr_cpp, BDD *system, unsigned int no_states, unsigned int no_inputs, bool optimized){

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

	// Create the Existental BDD(s) and ADD(s).
	createExistentalBDD();
	createExistentalADD();

	// System is being analyzed this way.
	this->system_analyzed = true;

	// If optimized
	if (optimized){
		this->optimized = false;
		initMinterms();
		this->optimized = true;
	}
	else{
		this->optimized = false;
	}

	// Get the System's BDD.
	system_bdd = system;


	printf("-\n");
	printf("Self test:\n");
	printf("No States Vars: %3d - No Inputs Vars: %3d\n", this->no_state_vars, this->no_input_vars);
	printf("x  index begin: %3d - end: %3d\n", bdd_x[0].getNode()->index,  bdd_x[bdd_x.size()-1].getNode()->index);
	printf("u  index begin: %3d - end: %3d\n", bdd_u[0].getNode()->index,  bdd_u[bdd_u.size()-1].getNode()->index);
	printf("x' index begin: %3d - end: %3d\n", bdd_x_[0].getNode()->index, bdd_x_[bdd_x_.size()-1].getNode()->index);
	printf("-\n");

} /* end of ShortestPath */


//! ShortestPath De-Constructor.
/**
 * Nothing special yet. :)
 */
ShortestPath::~ShortestPath() {

} /* end of ~ShortestPath */


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
	BDD S = system->ExistAbstract(bddExistental_xx_, 0);


//	ADD test = SC_swapped & S.Add();
//
//	// Create .dot file
//	std::vector<ADD> nodes_add;
//	FILE *outfile;
//	nodes_add.push_back(test);
//	outfile = fopen("Costs_.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();


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
}/* end of createCostAdjacencyMatrix */

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


	if (optimized){

		if ((*x)[0].getNode()->index == bdd_x[0].getNode()->index){
			return bdd_mterm_x[node];
		}
		else{
			return bdd_mterm_x_[node];
		}
	}
	else{
		unsigned int j;

		BDD one  = mgr_cpp->bddOne();
		BDD zero = mgr_cpp->bddZero();

		BDD minterm, temp;


		minterm = one;

#ifdef LSB_MSB
		for (j = 0; j < (*x).size(); j++){
#endif
#ifdef MSB_LSB
		for (j = (*x).size();;){
			if (j==0) break; j--;
#endif
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
	}
} /* end of createMinterm */

//!
inline BDD ShortestPath::createMinterm(std::vector<BDD> *x, std::vector<BDD> *y, unsigned int node_x, unsigned int node_y){
	unsigned int j;

	BDD one  = mgr_cpp->bddOne();
	BDD zero = mgr_cpp->bddZero();

	BDD minterm, temp;


	minterm = one;

#ifdef LSB_MSB
	for (j = 0; j < (*x).size(); j++){
#endif
#ifdef MSB_LSB
	for (j = (*x).size();;){
		if (j==0) break; j--;
#endif
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
} /* end of createMinterm */

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


	if (optimized){

		if ((*x)[0].getNode()->index == add_x[0].getNode()->index){
			return add_mterm_x[node];
		}
		else{
			return add_mterm_x_[node];
		}

	}
	else{

		unsigned int j;

		ADD one  = mgr_cpp->addOne();
		ADD minterm, temp;

		minterm = one;

#ifdef LSB_MSB
		for (j = 0; j < (*x).size(); j++){
#endif
#ifdef MSB_LSB
		for (j = (*x).size();;){
			if (j==0) break; j--;
#endif
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
	}
} /* end of createMinterm */


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

#ifdef LSB_MSB
	for (j = 0; j < (*x).size(); j++){
#endif
#ifdef MSB_LSB
	for (j = (*x).size();;){
		if (j==0) break; j--;
#endif
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
#ifdef LSB_MSB
	for (j = 0; j < (*y).size(); j++){
#endif
#ifdef MSB_LSB
	for (j = (*y).size();;){
		if (j==0) break; j--;
#endif
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
} /* end of createMinterm */

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
} /* end of createVariables */


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
} /* end of createVariables */

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
} /* end of createVariables */


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
} /* end of getVarsIndex */


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
} /* end of getVarsIndex */


//! Initializes the necessary BDDs needed to perform various "Existental" operations.
void ShortestPath::createExistentalBDD(){

	// Get rid of x'. Keep only x and u.
	int existental[mgr_cpp->ReadSize()];
	for (unsigned k = 0; k < bdd_u[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_u[0].getNode()->index; k < bdd_x_[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_x_[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 1;
	}
	DdNode *cube_array  = Cudd_CubeArrayToBdd(mgr,existental);
	bddExistental_xu = BDD(*mgr_cpp, cube_array);

	// Get rid of u. Keep only x and x'.
	for (unsigned k = 0; k < bdd_u[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_u[0].getNode()->index; k < bdd_x_[0].getNode()->index; k++){
		existental[k] = 1;
	}
	for (unsigned k = bdd_x_[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 2;
	}
	DdNode *cube_array1  = Cudd_CubeArrayToBdd(mgr,existental);
	bddExistental_xx_ = BDD(*mgr_cpp, cube_array1);

	// Get rid of u and x'. Keep only x.
	for (unsigned k = 0; k < bdd_u[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_u[0].getNode()->index; k < bdd_x_[0].getNode()->index; k++){
		existental[k] = 1;
	}
	for (unsigned k = bdd_x_[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 1;
	}
	DdNode *cube_array2  = Cudd_CubeArrayToBdd(mgr,existental);
	bddExistental_x = BDD(*mgr_cpp, cube_array2);

} /* end of createExistentalBDD */

//! Initializes the necessary ADDs needed to perform various "Existental" operations.
void ShortestPath::createExistentalADD(){

	// Get rid of x'. Keep only x and u.
	int existental[mgr_cpp->ReadSize()];
	for (unsigned k = 0; k < bdd_u[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_u[0].getNode()->index; k < bdd_x_[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_x_[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 1;
	}
	DdNode *cube_array  = Cudd_CubeArrayToBdd(mgr,existental);
	Cudd_Ref(cube_array);
	DdNode *result = Cudd_BddToAdd(mgr, cube_array);
	addExistental_xu = ADD(*mgr_cpp, result);
	Cudd_RecursiveDeref(mgr, cube_array);

	// Get rid of u and x'. Keep only x.
	for (unsigned k = 0; k < bdd_u[0].getNode()->index; k++){
		existental[k] = 2;
	}
	for (unsigned k = bdd_u[0].getNode()->index; k < bdd_x_[0].getNode()->index; k++){
		existental[k] = 1;
	}
	for (unsigned k = bdd_x_[0].getNode()->index; k < (unsigned int)(mgr_cpp->ReadSize()); k++){
		existental[k] = 1;
	}
	DdNode *cube_array1  = Cudd_CubeArrayToBdd(mgr,existental);
	Cudd_Ref(cube_array1);
	DdNode *result1 = Cudd_BddToAdd(mgr, cube_array1);
	addExistental_x = ADD(*mgr_cpp, result1);
	Cudd_RecursiveDeref(mgr, cube_array1);

} /* end of createExistentalBDD */


void ShortestPath::initMinterms(){

	printf("ShortestPath::initMinterms().\n");

	for(unsigned int i = 0; i < no_states; i++){
		bdd_mterm_x.push_back(createMinterm(&bdd_x, i));
		bdd_mterm_x_.push_back(createMinterm(&bdd_x_, i));

		add_mterm_x.push_back(createMinterm(&add_x, i));
		add_mterm_x_.push_back(createMinterm(&add_x_, i));
	}

	for(unsigned int i = 0; i < no_inputs; i++){
		bdd_mterm_u.push_back(createMinterm(&bdd_u, i));
	}
} /* end of initMinterms */

//! Finds the shortest path from all pairs to a given target set W. Returns the vector containing the shortest path values and the pointer vector.
/**
 * This method takes as input the DDs containing the all-pairs shortest path values, the pointer array of the all-pairs shortest path and a target set W,
 * for which we want to find the shortest path from all the pairs to that set. The set is given as vector of integers, denoting the states. The method
 * returns the vector containing the shortest path value as ADD and the pointer vector that shows which node to follow to achieve the shortest path value.\n
 *
 * __Important Notice__:
 * - Memory should have been already allocated for the results (<strong class="paramname">APSP_W</strong> and <strong class="paramname">PA_W</strong>).
 * - The pointer does not contain the "map" to achieve the shortest path value from all pairs not the set. It only contains the intermediate node
 * to be followed in order to achieve the minimum path. Therefore this result should be used together with the initial pointer array (<strong class="paramname">PA</strong>).
 *
 * @param APSP is the ADD containing the all-pairs shortest path values.
 * @param PA is the ADD containing the pointer array of the APSP.
 * @param W	is the desired target set given as a BDD.
 * @param APSP_W is the returned ADD, containing the all-pairs shortest path values to the set W.
 * @param PA_W is the returned ADD, containing the intermediate nodes to be followed to achieve the all-pairs shortest path to the set W. This result should be used
 *        together with the PA ADD.
 * @see   FloydWarshall, APtoSetSP, Cudd_addApplyMin2, cuddAddApplyMin2Recur
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

} /* end of APtoSetSP */

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

#ifdef LSB_MSB
	for (j = 0; j < no_state_vars; j++){
#endif
#ifdef MSB_LSB
	for (j = no_state_vars;;){
			if (j==0) break; j--;
#endif
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
} /* end of createXstates */


//! This method "relaxes" the states, i.e. updates the shortest path value and the pointer index of the states that are consideres as candidates for the Z set, i.e. the resolved set.
inline void ShortestPath::relax(BDD *XUz, ADD *APSP_W, BDD *PA_W, ADD *SC, pq_relax *pq_mincost, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_, std::vector<ADD> *add_x, std::vector<ADD> *add_x_){
#ifdef APtoSetSP_DEBUG
	printf("ShortestPath::relax\n");
#endif

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif


	//
	BDD XUz_swapped;
	BDD validInput, validTransition;
	BDD bdd_minterm_x;
	BDD bdd_zero;

	ADD XUz_add;
	ADD Dwx, Dwx_, DwC, c_xux_, c_xux_swapped;
	ADD Cdet_u;
	ADD add_minterm_x;
	ADD Cmin, Cndet_max, Cndet_u_max, Cdet_min;

	//
	DdNode *validInput_;

	bool non_deterministic = false;
	bool deterministic = false;
	bool add_nondet_u  = false;


	// Create .dot file
	std::vector<BDD> nodes_bdd;
	std::vector<ADD> nodes_add;
	FILE *outfile;

	bdd_zero = mgr_cpp->bddZero();

	XUz_swapped = XUz->SwapVariables(*bdd_x, *bdd_x_);

	// Get Dw(x).
	Dwx = XUz->ExistAbstract(bddExistental_x).Add() & (*APSP_W);

	// Get Dw(x')
	Dwx_ = XUz_swapped.Add() & (*APSP_W);

	// Get c(x,u,x'). This is basically the cost of the x'. Since SC is in the form of x,
	// we have to swap XUz.
	c_xux_swapped = XUz_swapped.Add() & (*SC);
	// Now, we have the c(x,u,x') in the form of x'->u->x so we have to swap x and x'.
	c_xux_ = c_xux_swapped.SwapVariables(*add_x, *add_x_);



	// With this trick we are going to which transitions are non-deterministic,
	// as in this case the ExistAbstract() function sums their cost.
	ADD Dw_plus_C_swapped = c_xux_swapped + Dwx_;
	ADD Dw_plus_C         = Dw_plus_C_swapped.SwapVariables(*add_x, *add_x_);
	ADD c_xux_e           = Dw_plus_C.ExistAbstract(addExistental_xu);

	// So, we know have an ADD in the form of (x,u), where the costs for the non-deterministic
	// transitions have been summed up and the costs of the remaining (deterministic) transitions
	// remain as they were. So, using Xnor we are going to see where the initial ADD (c_xux_), that
	// holds all transitions differs with the c_xux_e where the ExistAbstract() has been applied.
	// This difference can be used to break down the initial ADD (c_xux_) into non-deterministic and
	// deterministic transitions.
	ADD c_xux_ee   = c_xux_e.Xnor(Dw_plus_C);
	ADD c_xux_det  = c_xux_ee * Dw_plus_C;
	ADD c_xux_ndet = Dw_plus_C - c_xux_det;

	// Again we are getting rid of the x' from both the on-deterministic and deterministic c(x,u,x').
	// But in the case of the non-deterministic we are taking the maximum. This is done by a modified
	// ExistAbstract() function, that does not sum the non-deterministic transitions but takes their maximum.
	ADD c_xux_det_xu  = c_xux_det.ExistAbstract(addExistental_xu);
	ADD c_xux_ndet_xu = ADD(*mgr_cpp, CuddaddExistAbstract(mgr, c_xux_ndet.getNode(), addExistental_xu.getNode()));


	// With this trick we are creating an ADD that has two values: 1 and inf. 1 is for costs that exists, i.e.
	// valid given costs, and inf for the rest, i.e. invalid states. This is used to filter the deterministic
	// c(x,u,x'), which is in the form of (x,u). We are doing this, because for the invalid inputs, i.e. inputs
	// that do not yield a transition to the Z set, their corresponding values is zero (0). But this will cause
	// problems if we want to find the minimum of all (valid) transitions. So, we have to convert this zero to inf,
	// to apply the Minimum function.
	DdNode *covert = covertSomeValue(&c_xux_det_xu, mgr_cpp->addZero(), mgr_cpp->plusInfinity());
	c_xux_det_xu = ADD(*mgr_cpp, covert);



	// This might speed things up.
	BDD xStates = XUz->ExistAbstract(bddExistental_x);
	double avail_states = xStates.CountMinterm(no_state_vars);
	double count_states = 0.0;

	/* Iterate over all states. */
#ifdef APtoSetSP_DEBUG
	printf("Number of states to be relaxed: %d\n", (int)avail_states);
#endif
	for (unsigned int i = no_states;;){
		if (i==0) break; i--;

		add_minterm_x = createMinterm(add_x, i);

		// Deterministic c(x,u,x').
		Cdet_u = c_xux_det_xu.Restrict(add_minterm_x);
		// Non-deterministic c(x,u,x').
		Cndet_u_max = c_xux_ndet_xu.Restrict(add_minterm_x);

		// Check if we have deterministic transitions...
		if (Cdet_u == mgr_cpp->plusInfinity()){
			deterministic = false;
		}
		else deterministic = true;
		// how about non-deterministic?
		if (Cndet_u_max.IsZero()){
			non_deterministic = false;
		}else non_deterministic = true;

		if (!non_deterministic && !deterministic){
			continue;
		}

#ifdef APtoSetSP_DEBUG
		printf("-Checking State: %d\n", i);
#endif


		if (deterministic && non_deterministic){
			// Deterministic transitions.
			DdNode *result = Cudd_addFindMin(mgr, Cdet_u.getNode()); 	// Cdet_min  = Cdet_u.FindMin();
			Cdet_min = ADD(*mgr_cpp, result);
			// Non-deterministic transition.
			Cndet_max = Cndet_u_max.FindMax();
#ifdef APtoSetSP_DEBUG
			printf("   All types of transitions available. Cdet_min = %.2f and Cndet_max = %.2f\n", Cdet_min.getNode()->type.value, Cndet_max.getNode()->type.value);
#endif
			if (Cdet_min < Cndet_max){
				Cmin = Cdet_min;
			}
			else if (Cdet_min == Cndet_max){
#ifdef APtoSetSP_DEBUG
				printf("      Deterministic and non-deterministic transitions have the same cost.\n");
#endif
				Cmin = Cndet_max;
				add_nondet_u = true;
			}
			else{
				add_nondet_u = true;
				deterministic = false;
				Cmin = Cndet_max;
			}
		}
		else if (non_deterministic){
			add_nondet_u = true;
			Cndet_max = Cndet_u_max.FindMax();
			Cmin = Cndet_max;
#ifdef APtoSetSP_DEBUG
			printf("   Only non-deterministic transitions available. Cndet_max = %.2f\n", Cndet_max.getNode()->type.value);
#endif
		}
		else{
			// Deterministic transitions.
			DdNode *result = Cudd_addFindMin(mgr, Cdet_u.getNode()); 	// Cdet_min  = Cdet_u.FindMin();
			Cdet_min = ADD(*mgr_cpp, result);
			Cmin = Cdet_min;
#ifdef APtoSetSP_DEBUG
			printf("   Only deterministic transitions available. Cdet_min = %.2f\n", Cdet_min.getNode()->type.value);
#endif
		}


//		/* Dw(x) > Dw(x') + c(x,u,x') */
//		// Get Dw(x')
//		Dwx_ = XUz->Restrict(add_minterm_x).Add() & (*APSP_W);


//		DwC = Dwx_.Restrict(add_minterm_x.SwapVariables(*add_x, *add_x_)).FindMax() + Cmin;
		DwC = Cmin;
		if (Dwx.Restrict(add_minterm_x) >= DwC){
#ifdef APtoSetSP_DEBUG
			printf(" Found shorter path!\n");
			printf(" Adding state %d with cost %.2f to the priority queue. \n", i, DwC.getNode()->type.value);
#endif
			bdd_minterm_x = createMinterm(bdd_x, i);
			// Update Priority Queue.
			pq_mincost->push(std::make_pair(DwC.getNode()->type.value, i));

			if (deterministic){
				validInput_ = getTransitionFromValue(&Cdet_u, &Cmin);
				validInput  = BDD(*mgr_cpp, validInput_);

				// Creating temp transition for the PA_W
				validTransition = bdd_minterm_x.Ite(validInput, bdd_zero);

				// Update Dw(x)
				*APSP_W &= (~add_minterm_x);
				*APSP_W += (add_minterm_x * DwC);

				/* Update PA_W. */
				// Remove older pointer entry
				*PA_W &= (~bdd_minterm_x);
				// Add new entry
				*PA_W += validTransition;
			}

			if(non_deterministic && add_nondet_u){
				validInput_ = getTransitionFromValue(&Cndet_u_max, &Cndet_max);
				validInput  = BDD(*mgr_cpp, validInput_);

				// Creating temp transition for the PA_W
				validTransition = bdd_minterm_x.Ite(validInput, bdd_zero);

				// Remove older entries from APSP_W and PA_W.
				if (!deterministic){
					*APSP_W &= (~add_minterm_x);
					*PA_W   &= (~bdd_minterm_x);
					// Update Dw(x)
					*APSP_W += (add_minterm_x * DwC);
				}

				// Add new entry
				*PA_W   += validTransition;
			}

			deterministic = false;
			non_deterministic = false;
			add_nondet_u = false;

		}

		count_states++;
		// This might speed things up.
		if (count_states == avail_states){
			break;
		}
	}


//	// Create .dot file
//	nodes_bdd.push_back(testa);
//	outfile = fopen("Ztest.dot", "w");
//	mgr_cpp->DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//	break;
//
//
//	nodes_add.push_back(Dwx);
//	outfile = fopen("Dwx.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(Dwx_);
//	outfile = fopen("Dwx_.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_);
//	outfile = fopen("c_xux_.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_e);
//	outfile = fopen("c_xux_e.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_det);
//	outfile = fopen("c_xux_det.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_ndet);
//	outfile = fopen("c_xux_ndet.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_det_xu);
//	outfile = fopen("c_xux_det_xu.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//
//	nodes_add.push_back(c_xux_ndet_xu);
//	outfile = fopen("c_xux_ndet_xu.dot", "w");
//	mgr_cpp->DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();

//	nodes_add.push_back(Dw_c_xu);
//	outfile = fopen("Dw_c_xu.dot", "w");
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


#ifdef APtoSetSP_DEBUG
#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::relax (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
#endif
} /* end of relax */


//! This method both implements the Xsz and Usz operators as defined in the theory.
inline BDD ShortestPath::operatorXUsz(BDD *W, BDD *W_swapped, BDD *Q, BDD *Z, std::vector<BDD> *bdd_x, std::vector<BDD> *bdd_u, std::vector<BDD> *bdd_x_){
#ifdef APtoSetSP_DEBUG
	printf("ShortestPath::operatorXUsz\n");
#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif
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
	BDD system_rstct_Zxu = system_rstct_Z.ExistAbstract(bddExistental_xu, 0);


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
	BDD rg_ns_x = rg_ns_xux_.ExistAbstract(bddExistental_xu, 0);


	// Now, remove undesired states. Of course states in the Z set have also to be removed.
	// This is basically the BDD with the desired (x,u). :)
	BDD XUsz = system_rstct_ZW & (~rg_ns_x) & (~(*Z));
	// We need also the x' states. So final form: (x,u,x').
	XUsz &= (*system_bdd);


#ifdef APtoSetSP_DEBUG
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
#endif


	// Create .dot file
	std::vector<BDD> nodes_bdd;
	FILE *outfile;
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

#ifdef APtoSetSP_DEBUG
#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::operatorXUsz (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif
#endif


	return XUsz;
} /* end of operatorXUsz */

//! This method initialized the PA_W (BDD) to include the states that belong to the target set W and their inputs. PA_W represents the refined system/controller.
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
/**
 * This method is used mainly to solve the non-deterministic shortest path problem. It is based on the reachability game that is performed by the operators \f$X_{S_{Z}}\f$ and \f$U_{S_{Z}}\f$ (operatorXUsz).
 * These operators point out the next valid state-input(s), i.e. the next candidate (with its valid inputs) for the Z set. That is the set that holds the states where the shortest path has been computed.
 * For each of these states the relax() function is called, which updates the shortest path value towards the target set and adds these states to a priority queue. Based on that queue the state with the
 * lowest cost value is added to the Z set and is marked as resolved. The process ends when all states have been resolved. \n
 * This method supports also deterministic systems.
 *
 *__Important Notice__:
 * - This algorithm assumes that the liveness constraints of the system/controller have been already solved. That is, it is guaranteed that the target set <strong class="paramname">W</strong> can been reached by all
 *	 states of the system/controller.
 * - This method assumes that the arguments (<strong class="paramname">APSP_W</strong>, <strong class="paramname">PA_W</strong>), which are used to pass the result, have been already allocated. Empty pointers of these will result to an error!
 *
 * @param S is the pointer to the System's BDD.
 * @param SC is the pointer to the System's costs ADD.
 * @param W is the pointer to the target set.
 * @param APSP_W is the _allocated_ pointer that will hold / will store the all-pairs to a target set shortest path values. (Is used to return the result)
 * @param PA_W is the _allocated_ pointer the will hold / will store the new refined system. (Is used to return the result)
 * @param no_states is the number of states of the system.
 * @param no_inputs is the number of inputs of the system.
 * @see operatorXUsz, relax
 */
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
		createVariables(vars_index, no_state_vars, no_input_vars, add_x, add_x_);

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
#ifdef APtoSetSP_DEBUG
	int iteration = 0;
	int db_state;
	double db_cost;
#endif

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
	 * Main Loop. Breaks when Q != empty.
	 */
	while(left_states){

#ifdef APtoSetSP_DEBUG
		iteration++;
		printf("*Iteration %d\n", iteration);
#endif

//		printf("%.2f%% Queue size: %d. Top state: %d\n", (1.0 - (double)left_states/no_states)*100.0, pq_mincost.size(), pq_mincost.top().second);
		if (pq_mincost.size()==0){
			break;
		}

		// Get the state with the minimum cost: x = min{Dw(x) | x \in Q}
		x_temp = createMinterm(bdd_x, pq_mincost.top().second);
#ifdef APtoSetSP_DEBUG
		db_state = pq_mincost.top().second;
		db_cost  = pq_mincost.top().first;
#endif
		pq_mincost.pop(); // Remove lowest priority state
		// If it is already resolved... continue. (Priority Queue unpleasant property)
		if (!(Z.Restrict(x_temp).IsZero())){
			continue;
		}
		left_states--;

#ifdef APtoSetSP_DEBUG
		printf("==>Adding state %d to Z. Cost: %.3f<==\n", db_state, db_cost);
#endif

		// Z = Z \cup x
		Z = Z + x_temp;

		// Q = Q\x (set minus)
		Q = Q - x_temp;


#ifdef APtoSetSP_DEBUG
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
#endif



		// Operator XUsz. Implements both the Xsz and Usz operators.
		XUz = operatorXUsz(W, &W_swapped, &Q, &Z, bdd_x, bdd_u, bdd_x_);

		if (XUz.IsZero()){
			continue;
		}

		// Relax
		relax(&XUz, APSP_W, PA_W, SC, &pq_mincost, bdd_x, bdd_u, bdd_x_, add_x, add_x_);


#ifdef APtoSetSP_DEBUG
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
#endif
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



	// De-Allocate memory for the vectors.
	if (!system_analyzed){
		delete bdd_x;
		delete bdd_u;
		delete bdd_x_;
		delete add_x;
		delete add_x_;
	}



#ifdef ENABLE_TIME_PROFILING
	/* Print execution time. */
	printf("ShortestPath::APtoSetSP (ND): Execution Time: %ds (%dms) (%dus)\n",  (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#endif

} /* end of APtoSetSP */

//! This method finds the next state to be followed to achieve the shortest path. That is, if there is no direct link.
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
} /* end of findSequentNode */


//!
//! Creates a BBD of the new refined controller in form of (x,u), based on the old one and the results of the Deterministic Shortest Path (@ref APtoSetSP).
/**
 * This method first checks how to get to the target node, by looking at the FW pointer array. If it is zero, it means that we have a direct link and nothing has to be done.
 * Otherwise it finds the next (subsequent) node towards the target set, by calling the findSequentNode() method. In any case, as soon as it finds a valid states it records also
 * the valid inputs. The final result is a BDD of the refined system in the form of (x,u).
 *
 * @param S a pointer to the BDD of the initial system, i.e. the controller that needs to be refined.
 * @param APSP_PA a pointer to the ADD representing the all-pairs to a target set shortest path values.
 * @param APSP_PA_W a pointer to the ADD representing the all-pairs to a target set shortest path pointer array.
 * @return a BDD of the refined system in the form of (x,u).
 * @see APtoSetSP(ADD *APSP, ADD *PA, BDD *W, ADD *APSP_W, ADD *PA_W), findSequentNode
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
} /* end of createControllerBDD */


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
} /* end of Dddmp_cuddStore */


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
} /* end of Dddmp_cuddStore */


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
* @see createCostAdjacencyMatrix, AddOuterSumTrace, AddOuterSumRecurTrace
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
#ifdef LSB_MSB
		for (j = 0; j < no_vars; j++){
#endif
#ifdef MSB_LSB
	for (j = no_vars;;){
		if (j==0) break; j--;
#endif
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
} /* end of FloydWarshall */


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
} /* end of AddOuterSumTrace */



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

} /* end of AddOuterSumRecurTrace */




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
} /* end of Cudd_addApplyTrace */


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
} /* end of cuddAddApplyRecurTrace */


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
} /* end of Cudd_addApplyTrace */

//! Re-implemented method of the Cudd_addApply function of the cudd library, to be used for the @ref Cudd_addApplyMin2 method.

//!Implements the recursive step of the @ref Cudd_addApplyMin2.
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

} /* end of cuddAddApplyRecurTrace */



//!
DdNode *ShortestPath::covertSomeValue(ADD *f, ADD from, ADD to){
	DdNode *res;

	do {
		mgr_cpp->getManager()->reordered = 0;
		res = covertSomeValueRecur(f->getNode(), from.getNode(), to.getNode());
    } while (mgr_cpp->getManager()->reordered == 1);

	return res;
} /* end of covertSomeValue */


//!
DdNode *ShortestPath::covertSomeValueRecur(DdNode *f, DdNode *from, DdNode *to) {

	 DdNode *res,
		    *fv, *fvn,
		    *T, *E;

	 unsigned int index;

	statLine(dd);


	// Change the values.
    if (cuddIsConstant(f)){

    	if (f == from){
    		return to;
    	}

    	return f;
    }


	fv  = cuddT(f);
	fvn = cuddE(f);

	index = cuddI(mgr,f->index);

	// True
	T = covertSomeValueRecur(fv, from, to);
	if (T == NULL ){
		return NULL;
	}
	cuddRef(T);


	// Else
	E = covertSomeValueRecur(fvn, from, to);
	if (E == NULL ){
		Cudd_RecursiveDeref(mgr_cpp->getManager(),T);
		return NULL;
	}
	cuddRef(E);


	// If it exists return it, otherwise create the node.
	res = (T == E) ? T : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T,E);
	if (res == NULL) {
		printf("cuddUniqueInter FAIL!\n");
		Cudd_RecursiveDeref(mgr_cpp->getManager(), T);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), E);
		return NULL;
	}
	cuddDeref(T);
	cuddDeref(E);

	return res;

} /* end of cuddAddApplyRecurTrace */


//!
DdNode *ShortestPath::getTransitionFromValue(ADD *Cu, ADD *value){

	DdNode *res;
	bool found = false;

	do {
		mgr_cpp->getManager()->reordered = 0;
		res = getTransitionFromValueRecur(Cu->getNode(), value->getNode(), &found);
    } while (mgr_cpp->getManager()->reordered == 1);

//	printf("returning from getTransitionFromValue\n");

	if (res == NULL){
		return DD_ONE(mgr);
	}
//	printf("valid bdd\n");
	return res;
} /* end of getInputFromValue */


//!
DdNode *ShortestPath::getTransitionFromValueRecur(DdNode *f, DdNode *value, bool *found){


    DdNode *fv, *fvn,
	       *T, *E,
	       *res;
    unsigned int index;
//    unsigned int topf;

    DdNode *one, *zero;

    one  = DD_ONE(mgr_cpp->getManager());
    zero = Cudd_Not(one);

//    printf("IN\n");


    statLine(dd);
//    res = NULL;
    // Checking if desired value has been found.
    if (cuddIsConstant(f)){

//    	printf("f is constant\n");
    	if (f == value){

//    		printf("Found the value! Returning!\n");
    		*found = true;
    		return one;
    	}
    	return zero;
    }



	// Find the top variable.
//    F = Cudd_Regular(f);
//    printf("index: %d - %d\n", f->index, F->index);
//    topf = mgr_cpp->getManager()->perm[F->index];
//    topf  = cuddI(mgr,f->index);


//	printf("topf = %d\n", topf);


	/* Compute cofactors. */
	fv  = cuddT(f);
	fvn = cuddE(f);
	index = cuddI(mgr,f->index);


//	printf("TRUE\n");
	// True
	T = getTransitionFromValueRecur(fv, value, found);
//	printf("AFTER TRUE\n");
	if (T == NULL ){
		return NULL;
	}
	cuddRef(T);

	if (!(*found)){
//		printf("FALSE\n");
		// Else
		E = getTransitionFromValueRecur(fvn, value, found);
//		printf("AFTER FALSE\n");
		if (E == NULL ){
			Cudd_RecursiveDeref(mgr_cpp->getManager(),T);
			return NULL;
		}
		cuddRef(E);

		if (*found){
//			Cudd_RecursiveDeref(mgr_cpp->getManager(),T);
			T = zero;
			cuddRef(T);
		}
	}
	else{
		E = zero;
		cuddRef(E);
	}

	if (*found){

//		printf("creating node at %d ...\n", index);
		// If it exists return it, otherwise create the node.
//		res = (T == E) ? T : cuddUniqueInter(mgr_cpp->getManager(),(int)index,T,E);
//		if (res == NULL) {
//			printf("cuddUniqueInter FAIL!\n");
//			Cudd_RecursiveDeref(mgr_cpp->getManager(), T);
//			Cudd_RecursiveDeref(mgr_cpp->getManager(), E);
//			return NULL;
//		}
		DdNode *x = Cudd_bddIthVar(mgr, index);
		cuddRef(x);
		res = Cudd_bddIte(mgr, x, T, E);
		cuddDeref(T);
		cuddDeref(E);
		cuddDeref(x);

		return res;
	}

	return NULL;

} /* end of getInputFromValue */


/*---------------------------------------------------------------------------*/
/* Definition of exported functions                                          */
/*---------------------------------------------------------------------------*/

/**Function********************************************************************

  Synopsis    [Existentially Abstracts all the variables in cube from f.]

  Description [Abstracts all the variables in cube from f by summing
  over all possible values taken by the variables. Returns the
  abstracted ADD.]

  SideEffects [None]

  SeeAlso     [Cudd_addUnivAbstract Cudd_bddExistAbstract
  Cudd_addOrAbstract]

******************************************************************************/
DdNode *ShortestPath::CuddaddExistAbstract(DdManager * manager, DdNode * f, DdNode * cube)
{
    DdNode *res;

    DdNode *two = cuddUniqueConst(manager,(CUDD_VALUE_TYPE) 2);
    if (two == NULL) return(NULL);
    cuddRef(two);

    if (CuddaddCheckPositiveCube(manager, cube) == 0) {
        (void) fprintf(manager->err,"Error: Can only abstract cubes");
        return(NULL);
    }

    do {
	manager->reordered = 0;
	res = cuddAddExistAbstractRecur(manager, f, cube, two);
    } while (manager->reordered == 1);

    if (res == NULL) {
	Cudd_RecursiveDeref(manager,two);
	return(NULL);
    }
    cuddRef(res);
    Cudd_RecursiveDeref(manager,two);
    cuddDeref(res);

    return(res);

} /* end of Cudd_addExistAbstract */


/**Function********************************************************************

  Synopsis    [Checks whether cube is an ADD representing the product
  of positive literals.]

  Description [Checks whether cube is an ADD representing the product of
  positive literals. Returns 1 in case of success; 0 otherwise.]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
int ShortestPath::CuddaddCheckPositiveCube(DdManager * manager, DdNode * cube)
{
    if (Cudd_IsComplement(cube)) return(0);
    if (cube == DD_ONE(manager)) return(1);
    if (cuddIsConstant(cube)) return(0);
    if (cuddE(cube) == DD_ZERO(manager)) {
        return(CuddaddCheckPositiveCube(manager, cuddT(cube)));
    }
    return(0);

} /* end of addCheckPositiveCube */


/**Function********************************************************************

  Synopsis    [Performs the recursive step of Cudd_addExistAbstract.]

  Description [Performs the recursive step of Cudd_addExistAbstract.
  Returns the ADD obtained by abstracting the variables of cube from f,
  if successful; NULL otherwise.]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
DdNode *ShortestPath::cuddAddExistAbstractRecur(DdManager * manager, DdNode * f, DdNode * cube, DdNode *two)
{
    DdNode	*T, *E, *res, *res1, *res2, *zero;

    statLine(manager);
    zero = DD_ZERO(manager);

    /* Cube is guaranteed to be a cube at this point. */
    if (f == zero || cuddIsConstant(cube)) {
        return(f);
    }

    /* Abstract a variable that does not appear in f => multiply by 2. */
    if (cuddI(manager,f->index) > cuddI(manager,cube->index)) {
	res1 = cuddAddExistAbstractRecur(manager, f, cuddT(cube), two);
	if (res1 == NULL) return(NULL);
	cuddRef(res1);
	/* Use the "internal" procedure to be alerted in case of
	** dynamic reordering. If dynamic reordering occurs, we
	** have to abort the entire abstraction.
	*/
	res = cuddAddApplyRecur(manager,Cudd_addTimes,res1,two);
	if (res == NULL) {
	    Cudd_RecursiveDeref(manager,res1);
	    return(NULL);
	}
	cuddRef(res);
	Cudd_RecursiveDeref(manager,res1);
	cuddDeref(res);
        return(res);
    }

//    if ((res = cuddCacheLookup2(manager, Cudd_addExistAbstract, f, cube)) != NULL) {
//	return(res);
//    }

    T = cuddT(f);
    E = cuddE(f);

    /* If the two indices are the same, so are their levels. */
    if (f->index == cube->index) {
	res1 = cuddAddExistAbstractRecur(manager, T, cuddT(cube), two);
	if (res1 == NULL) return(NULL);
        cuddRef(res1);
	res2 = cuddAddExistAbstractRecur(manager, E, cuddT(cube), two);
	if (res2 == NULL) {
	    Cudd_RecursiveDeref(manager,res1);
	    return(NULL);
	}
        cuddRef(res2);
	res = cuddAddApplyRecur(manager, Cudd_addMaximum, res1, res2);
	if (res == NULL) {
	    Cudd_RecursiveDeref(manager,res1);
	    Cudd_RecursiveDeref(manager,res2);
	    return(NULL);
	}
	cuddRef(res);
	Cudd_RecursiveDeref(manager,res1);
	Cudd_RecursiveDeref(manager,res2);
//	cuddCacheInsert2(manager, Cudd_addExistAbstract, f, cube, res);
	cuddDeref(res);
        return(res);
    } else { /* if (cuddI(manager,f->index) < cuddI(manager,cube->index)) */
	res1 = cuddAddExistAbstractRecur(manager, T, cube, two);
	if (res1 == NULL) return(NULL);
        cuddRef(res1);
	res2 = cuddAddExistAbstractRecur(manager, E, cube, two);
	if (res2 == NULL) {
	    Cudd_RecursiveDeref(manager,res1);
	    return(NULL);
	}
        cuddRef(res2);
	res = (res1 == res2) ? res1 :
	    cuddUniqueInter(manager, (int) f->index, res1, res2);
	if (res == NULL) {
	    Cudd_RecursiveDeref(manager,res1);
	    Cudd_RecursiveDeref(manager,res2);
	    return(NULL);
	}
	cuddDeref(res1);
	cuddDeref(res2);
//	cuddCacheInsert2(manager, Cudd_addExistAbstract, f, cube, res);
        return(res);
    }

} /* end of cuddAddExistAbstractRecur */


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



//! Experimental. No description.
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



//! Experimental. Analyzes given system's BDD and system's cost ADD.
void ShortestPath::diagnostics(BDD *S, ADD *costs){

	printf("\n***\nDiagnostics initiated...\n");

	unsigned int count = 0;

	for (unsigned int i = 0; i < no_states; i++){
		if (S->Restrict(createMinterm(&bdd_x,i)).IsZero()){


			if (S->Restrict(createMinterm(&bdd_x_,i)).IsZero()){
//				printf("***Warning! Invalid State: %d. ", i);
				if (costs->Restrict(createMinterm(&add_x,i)) == mgr_cpp->plusInfinity()){
//					printf("Cost infinity.\n");
				}
				count++;
				continue;
			}

			if (costs->Restrict(createMinterm(&add_x,i)) == mgr_cpp->plusInfinity()){
				printf("***Warning! Cost infinity.  State: %d\n", i);
			}
		}
	}

	printf("Diagnostics end!\nNumber of invalid states: %d of %d (%.2f%%). Only %d valid.\n***\n", count, no_states, 100*((double)(count)/(double)no_states), no_states - count);
}


//! Experimental. Filters out (assigns infinity) to states that are invalid.
ADD ShortestPath::filterCosts(BDD *S, ADD *costs, int mode){

	ADD valid_states  = mgr_cpp->addZero();

	if (mode == KEEP_VALID_STATES) {
		for (unsigned int i = 0; i < no_states; i++) {
			if (S->Restrict(createMinterm(&bdd_x, i)).IsZero()) {

				if (!S->Restrict(createMinterm(&bdd_x_, i)).IsZero()) {
					valid_states += createMinterm(&add_x, i);
				}
			} else {
				valid_states += createMinterm(&add_x, i);
			}
		}
	}
	else if (mode == KEEP_VALID_TRANSITIONS) {
		for (unsigned int i = 0; i < no_states; i++) {
			if (!S->Restrict(createMinterm(&bdd_x_, i)).IsZero()) {
				valid_states += createMinterm(&add_x, i);
			}
		}
	}

	return ((*costs) & valid_states);
}





