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

#include "ShortestPath.h"


//! ShortestPath Constructor.
/**
 * It assumes the CUDD manager has already been initialized.
 * @param mgr_cpp the pointer to the CUDD manager's object.
 */
ShortestPath::ShortestPath(Cudd *mgr_cpp) {

	this->mgr_cpp = mgr_cpp;
} /* ShortestPath */

//! ShortestPath De-Constructor.
/**
 * Nothing special yet. :)
 */
ShortestPath::~ShortestPath() {

} /* ~ShortestPath */

//! Create the ADD Target set W given the set of states of W as a vector<int>.
/**
 * Method takes as input the Systems ADD or any other ADD that contains the System's variables,
 * the vector containing the states of the target set and returns the ADD of the target set.
 * @param system is the pointer to the System's ADD or any other ADD that contains the System's variables.
 * @param no_states is the number of states of the System.
 * @param target_set is the vector containing the states of the target set W.
 * @return The ADD of the target set W.
 */
ADD ShortestPath::createTargetSet(ADD *system, int no_states, int no_inputs, std::vector<int> target_set){
	//
	std::vector<ADD> x;
	std::vector<ADD> y;
	ADD minterm;
	ADD one  = mgr_cpp->addOne();
	ADD zero = mgr_cpp->addZero();
	ADD cofactor;
	ADD temp;

	unsigned int i;

	printf("ShortestPath::createTargetSet: Number of states in W: %d\n", target_set.size());

	// Get the number of variables
	int no_state_vars = no_states/2 + no_states % 2;

	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(system);

	// Create the variables of the System.
	createVariables(vars_index, no_state_vars, no_inputs, &x, &y);

	cofactor = zero;
//	cofactor = mgr_cpp->background();

	for (i = 0; i < target_set.size(); i++){

		minterm = createMinterm(&y, target_set[i]);

		// Create the constant node.
		temp     = minterm.Ite(one, cofactor);
		cofactor = temp;
	}

	return cofactor;
} /* createTargetSet */

//! Create the ADD Target set W given the set of states of W as a vector<int>.
/**
 * Method takes as input the Systems BDD, the vector containing the states
 * of the target set and returns the ADD of the target set.
 * @param system is the pointer to the System's BDD.
 * @param no_states is the number of states of the System.
 * @param target_set is the vector containing the states of the target set W.
 * @return The ADD of the target set W.
 */
ADD ShortestPath::createTargetSet(BDD *system, int no_states, int no_inputs, std::vector<int> target_set){

	//
	std::vector<ADD> x;
	std::vector<ADD> y;
	ADD minterm;
	ADD one  = mgr_cpp->addOne();
	ADD zero = mgr_cpp->addZero();
	ADD cofactor;
	ADD temp;

	unsigned int i;

//	printf("ShortestPath::createTargetSet: No states: %d\n", target_set.size());

	// Get the number of variables
	int no_state_vars = no_states/2 + no_states % 2;

	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(system);

	// Create the variables of the System.
	createVariables(vars_index, no_state_vars, no_inputs, &x, &y);

	cofactor = zero;
//	cofactor = mgr_cpp->background();

	for (i = 0; i < target_set.size(); i++){

		minterm = createMinterm(&y, target_set[i]);

		// Create the constant node.
		temp     = minterm.Ite(one, cofactor);
		cofactor = temp;
	}

	return cofactor;
} /* createTargetSet */


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

	int i,j,k;
	ADD AG;

	//
	std::vector<BDD> x;
	std::vector<BDD> u;
	std::vector<BDD> x_;
	//
	std::vector<ADD> xx;
	std::vector<ADD> xx_;

	// Get the number of variables
	int no_state_vars = no_states/2 + no_states % 2;
	int no_input_vars = no_inputs/2 + no_inputs % 2;

	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(system);

	// Create the variables of the System.
	createVariables(vars_index, no_state_vars, no_input_vars, &x, &u, &x_);
	createVariables(vars_index, no_state_vars, no_input_vars, &xx, &xx_);

	// temporary variables
	BDD system_rstct_x;
	BDD system_rstct_u;
	BDD bdd_minterm;
	ADD add_minterm, add_temp;

	int total_iter = 0; // for statistical reasons. TODO: delete this.

	printf("ShortestPath::getCostAdjacencyMatrix: No states: %d - No inputs: %d \n", no_states, no_inputs);

	/* Creating the Cost Adjacency Matrix
	 *
	 * First check: if (x,u,x') -> 1
	 * then create: (x,x') = cost
	 *
	 */

	/* Initialize the AG matrix */
	AG = mgr_cpp->background();

	/* Iterate over all states. */
	for (i = no_states - 1; i >= 0; i--){

		// Create the x minterm.
		bdd_minterm = createMinterm(&x,i);

		// Restrict the System's BDD given that state.
		system_rstct_x = system->Restrict(bdd_minterm);
		total_iter++;

		/* Add the cost of each state to zero. */
		// Create the minterm.
		add_minterm = createMinterm(&xx, &xx_, i, i);
		// Create the constant node. Zero node
		AG = add_minterm.Ite(mgr_cpp->constant(0), AG);

		/* Iterate over all possible inputs */
		for (j = no_inputs - 1; j >=0 ; j--){
			total_iter++;

			// Create the u minterm.
			bdd_minterm = createMinterm(&u,j);

			// not a valid input.
			if ((system_rstct_x.Restrict(bdd_minterm)).IsZero()){
				continue;
			}
			// valid input
			else{

				system_rstct_u = system_rstct_x.Restrict(bdd_minterm);

				/* Iterate over all possible end states. */
				for (k = no_states - 1; k >= 0; k--){
					total_iter++;

					// Create the x' minterm.
					bdd_minterm = createMinterm(&x_,k);

					// not a valid input.
					if ((system_rstct_u.Restrict(bdd_minterm)).IsZero()){
						continue;
					}
					// valid input
					/* Creating (x,x') = cost */
					else{
						printf("%d,%d,%d\n", i,j,k);

						// Create the minterm.
						add_minterm = createMinterm(&xx, &xx_, i, k);

						// Create the constant node.
						AG = add_minterm.Ite(state_cost->Restrict(createMinterm(&xx_, k)), AG);

						break; // TODO: If it is deterministic then we can break here.
					}
				}
			}
		}
	}
	printf("ShortestPath::getCostAdjacencyMatrix: Total iterations: %d\n", total_iter);
	return AG;
} /* createCostAdjacencyMatrix */

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

	int j;

	BDD one  = mgr_cpp->bddOne();
	BDD zero = mgr_cpp->bddZero();

	BDD minterm, temp;


	minterm = one;

	for (j = (*x).size() - 1; j >=0; j--){
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

	int j;

	ADD one  = mgr_cpp->addOne();
	ADD minterm, temp;

	minterm = one;

	for (j = (*x).size() - 1; j >=0; j--){
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
	int j;

	ADD one  = mgr_cpp->addOne();
	ADD zero = mgr_cpp->addZero();

	ADD minterm, temp;

	minterm = one;

	for (j = (*x).size() - 1; j >=0; j--){
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
	for (j = (*y).size() - 1; j >=0; j--){
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
 * __Imporatant__: The memory for the vectors @param x, @param u and @param x_ is assumed to be already allocated.
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
//		printf("x: %d - x': %d\n", vars_index[i], vars_index[i + (no_state_vars + no_input_vars)]);
		(*x).push_back(mgr_cpp->bddVar(vars_index[i]));
		(*x_).push_back(mgr_cpp->bddVar(vars_index[i + (no_state_vars + no_input_vars)]));
	}

	for (i = 0; i <  no_input_vars; i++){
//		printf("u: %d\n", vars_index[i + no_state_vars]);
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
 * __Imporatant__: The memory for the vectors @param x, @param u and @param x_ is assumed to be already allocated.
 *
 * @param vars_index contains all available indexes.
 * @param no_state_vars is the number of state variables (or x) variables being used.
 * @param no_input_vars is the number of input variables (or u) variables being used.
 * @param x is used to return the created x variables.
 * @param x_ is used to return the created x_ variables.
 * @see createVariables
 */
void ShortestPath::createVariables(std::vector<int> vars_index, int no_state_vars, int no_input_vars, std::vector<ADD> *x, std::vector<ADD> *x_){

	int i;

	ADD one  = mgr_cpp->addOne();
	ADD zero = mgr_cpp->addZero();

	ADD tempX, tempX_;

	// Create the Variables.
	for (i = 0; i < no_state_vars; i++){
//		printf("x: %d - x': %d\n", vars_index[i], vars_index[i + (no_states + no_inputs)]);

		tempX  = mgr_cpp->addVar(vars_index[i]);
		tempX_ = mgr_cpp->addVar(vars_index[i + (no_state_vars + no_input_vars)]);

		tempX  = tempX.Ite(one, zero);
		tempX_ = tempX_.Ite(one, zero);

		(*x).push_back(tempX);
		(*x_).push_back(tempX_);
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



//! To be implemented. (if needed).
void ShortestPath::APtoSetSP(ADD *APSP, ADD *PA, ADD *W, ADD *APSP_W, ADD *PA_W){


} /* APtoSetSP */

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
void ShortestPath::APtoSetSP(ADD *APSP, ADD *PA, std::vector<int> W, ADD *APSP_W, ADD *PA_W){

	/* No point passing only one state since we have the APSP, but you never know :P */
	if (W.size() == 1){
		APSP_W = APSP;
		PA_W   = PA;
	}
	else{
		unsigned int i;
		//
		std::vector<ADD> x;
		std::vector<ADD> y;
		//
		ADD APSP_rstrct, PA_rstrct;

		// Get the index of the variables of the BDD, representing the system.
		std::vector<int> vars_index = getVarsIndex(APSP);

		//
		createVariables(vars_index, vars_index.size()/2, 0, &x, &y);

		// Create the ADD of the Target set W.
		ADD add_W = createTargetSet(APSP, vars_index.size(), 0, W);

		// Get the columns of the APSP and the PA ADD, pointed out by W.
		APSP_rstrct = APSP->Restrict(add_W);
		PA_rstrct   = PA->Restrict(add_W);

		ADD op1 = APSP_rstrct.Restrict(createMinterm(&y,W[0]));
		ADD P1  = PA_rstrct.Restrict(createMinterm(&y,W[0]));
		ADD op2;
		ADD P2;
		DdNode *result[2];

		/* Iterate over all states of the set W. */
		for (i = 1; i < W.size(); i++){
			op2 = APSP_rstrct.Restrict(createMinterm(&y,W[i]));
			P2  = PA_rstrct.Restrict(createMinterm(&y,W[i]));
			Cudd_addApplyMin2(op1.getNode(),op2.getNode(),P1.getNode(),P2.getNode(),result);
			op1 = ADD(*mgr_cpp, result[0]);
			P1  = ADD(*mgr_cpp, result[1]);
		}

		// Get the result
		*APSP_W = ADD(*mgr_cpp, result[0]);
		*PA_W   = ADD(*mgr_cpp, result[1]);
	}




} /* APtoSetSP */

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
void ShortestPath::FloydWarshall(ADD *AG, ADD *APSP, ADD *PA) {

	printf("Floyd Warshall Algorithm.\n");

//	printf("Number of variables: %d\n", param->vars);

	// C++ to C
	DdNode *AG_    = AG->getNode();
	DdManager *mgr = mgr_cpp->getManager();
//	cf_minterm * minterms = createMinterms(mgr, param);

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


//	long start_time = util_cpu_time();
	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(AG); // TODO: this is done twice if we call the getCostAdjacencyMatrix function
	no_vars = vars_index.size()/2;

	/* Memory allocation */
	xx = (DdNode **) malloc(sizeof(DdNode *) * no_vars);
	yy = (DdNode **) malloc(sizeof(DdNode *) * no_vars);

	/* "Copy" the AG matrix. */
	S = AG_;
	/* Initialize the Pointer array. */
	TR = Cudd_addConst(mgr, 0);

	/* Zero and One (constant) nodes. Used to create the minterms.*/
	one  = Cudd_ReadOne(mgr);
	zero = Cudd_ReadZero(mgr);
	Cudd_Ref(one);
	Cudd_Ref(zero);

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
//	printf("Iterating over all matrix elements (%d).\n", matrix_elements);

	for (i = 0; i < matrix_elements; i++){
//		printf("Node (%d).\n", i);

		element = i;
		xminterm = one;
		yminterm = one;

		/* Creating the minterms. */
		// always construct the ADD from the bottom to the top.
		for (j = no_vars - 1; j >=0; j--){
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

		S = P;
		TR = TR_temp;;
	}

	Cudd_Deref(S);

	/* Memory De-allocation. */
	free(xx);
	free(yy);

	/* Print execution time. */
//	printf("Floyd Warshall execution time: %s\n", util_print_time(util_cpu_time() - start_time));

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
//				printf("cuddAddOuterSumRecurTrace: return constant/			less %d\n", node);
				cuddDeref(R);
				Result[0] = R;
				Result[1] = cuddUniqueConst(mgr, node);
				return;
			} else {
//				printf("cuddAddOuterSumRecurTrace: return constant/nf\n");
				Cudd_RecursiveDeref(mgr, R);
				Result[0] = M;
				Result[1] = cuddUniqueConst(mgr, 0);
				return;
			}
		} else {
			printf("IN (node %d)\n", node);

			DdNode *min_result[2];
			Cudd_addApplyMinTrace(Cudd_addMinimum, R, M, R, min_result, node);
			cuddDeref(min_result[0]);
			cuddDeref(min_result[1]);

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
	//if (topM == v) {      T = (Tt == Te) ? Tt : Cudd_addIte(mgr,M,Tt,Te); }
	//else if (topr == v) { T = (Tt == Te) ? Tt : Cudd_addIte(mgr,r,Tt,Te); }
	//else if (topc == v) { T = (Tt == Te) ? Tt : Cudd_addIte(mgr,c,Tt,Te);}

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
//    	printf("***ShortestPath::cuddAddApplyRecurTrace: Cache hit!\n");
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
	if (F == G)
		return (F);

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
 * __Important__: The memory for the Result should be allocated before passing it as an argoument.
 *
 * @param f is the first operand.
 * @param g is the second opernd.
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
    unsigned int ford, gord, Pford, Pgord;
    unsigned int index0, index1;
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
    		if(cuddV(f) < cuddV(g)){
    			Result[1] = Pf;
    		}
    		else{
    			Result[1] = Pg;
    		}
    	}
    	return;
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

    /* Recursive step for f,g. */
	ford = cuddI(mgr_cpp->getManager(),f->index);
	gord = cuddI(mgr_cpp->getManager(),g->index);
	if (ford <= gord) {
		index0 = f->index;
		// It might happen that f is going to reach
		// a constant node first. In that case we
		// actually can go further down the DD.
		if(!cuddIsConstant(f)){
			fv  = cuddT(f);
			fvn = cuddE(f);
		}
		else{
			index0 = g->index;
			fv = fvn = f;
		}
	} else {
		index0 = g->index;
		fv = fvn = f;
	}
	if (gord <= ford) {
		// It might happen that g is going to reach
		// a constant node first. In that case we
		// actually can go further down the DD.
		if(!cuddIsConstant(g)){
			gv  = cuddT(g);
			gvn = cuddE(g);
		}
		else{
			gv = gvn = g;
		}
	} else {
		gv = gvn = g;
	}


    /* Recursive step for Pf, Pg. */
	Pford = cuddI(mgr_cpp->getManager(),Pf->index);
	Pgord = cuddI(mgr_cpp->getManager(),Pg->index);

	if (Pford <= Pgord) {
		index1 = Pf->index;
		// It might happen that Pf is going to reach
		// a constant node first. In that case we
		// actually can go further down the DD.
		if(!cuddIsConstant(Pf)){
			Pfv    = cuddT(Pf);
			Pfvn   = cuddE(Pf);
		}
		else{
			Pfv = Pfvn = Pf;
		}
	} else {
		index1 = Pg->index;
		Pfv = Pfvn = Pf;
	}
	if (Pgord <= Pford) {
		// It might happen that Pg is going to reach
		// a constant node first. In that case we
		// actually can go further down the DD.
		if(!cuddIsConstant(Pg)){
			Pgv  = cuddT(Pg);
			Pgvn = cuddE(Pg);
		}
		else{
			Pgv = Pgvn = Pg;
		}
	} else {
		Pgv = Pgvn = Pg;
	}


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
    res = (T0 == E0) ? T0 : cuddUniqueInter(mgr_cpp->getManager(),(int)index0,T0,E0);
    if (res == NULL) {
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
    res = (T1 == E1) ? T1 : cuddUniqueInter(mgr_cpp->getManager(),(int)index1,T1,E1);
    if (res == NULL) {
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


