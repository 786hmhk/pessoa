/*
 * ShortestPath.cpp
 *
 *  Created on: May 30, 2013
 *      Author: tanasaki
 */

#include "ShortestPath.h"

ShortestPath::ShortestPath(Cudd *mgr_cpp) {

	this->mgr_cpp = mgr_cpp;
}

ShortestPath::~ShortestPath() {

}

/**/
ADD ShortestPath::getCostAdjacencyMatrix(BDD *system, ADD *cost, int no_states, int no_inputs){

	int i,j,k,l;
	ADD AG;

	//
	BDD x[no_states];
	BDD u[no_inputs];
	BDD x_[no_states];
	//
	ADD xx[no_states];
	ADD xx_[no_states];

	// Get the index of the variables of the BDD, representing the system.
	std::vector<int> vars_index = getVarsIndex(system);

	// Create the variables of the System.
	createBDDVariables(system, vars_index, no_states, no_inputs, x, u, x_);
	createADDVariables(system, vars_index, no_states, no_inputs, xx, xx_);

	// temporary variables
	BDD system_rstct_x;
	BDD system_rstct_u;
	ADD minterm;
	bool init_AG = false;

	int total_iter = 0; // for statistical reasons. TODO: delete this.



	/* Creating the Cost Adjacency Matrix
	 *
	 * First check: if (x,u,x') -> 1
	 * then create: (x,x') = cost
	 *
	 * */
	// Iterate over all states.
	for (i = no_states - 1; i >= 0; i--){
		system_rstct_x = system->Restrict(x[i]);
		total_iter++;

		// iterate over all possible inputs
		for (j = no_inputs - 1; j >=0 ; j--){
			total_iter++;

			// not a valid input.
			if ((system_rstct_x.Restrict(u[j])).IsZero()){
				continue;
			}
			// valid input
			else{
				system_rstct_u = system_rstct_x.Restrict(u[j]);

				// iterate over all possible end states.
				for (k = no_states - 1; k >= 0; k--){
					total_iter++;
					// not a valid input.
					if ((system_rstct_u.Restrict(x_[k])).IsZero()){
						continue;
					}
					// valid input

					/* Creating (x,x') = cost */
					else{
						printf("%d,%d,%d\n", i,j,k);

						// Source States
						if (k == no_states - 1){
							minterm = xx_[no_states - 1].Ite(cost->Restrict(xx_[k]), mgr_cpp->plusInfinity());
						}
						else{
							minterm = xx_[no_states - 1].Ite(mgr_cpp->plusInfinity(), cost->Restrict(xx_[k]));
						}

						for (l = no_states - 2; l >= 0; l--){

							if (l != k){
								minterm = xx_[l].Ite(mgr_cpp->plusInfinity(), minterm);
							}
							else{
								minterm = xx_[l].Ite(minterm, mgr_cpp->plusInfinity());
							}
						}

						// Target States
						for (l = no_states - 1; l >= 0; l--){

							if (l != i){
								minterm = xx[l].Ite(mgr_cpp->plusInfinity(), minterm);
							}
							else{
								minterm = xx[l].Ite(minterm, mgr_cpp->plusInfinity());
							}
						}

						break; // TODO: If it is deterministic then we can break here.
					}
				}

				// Sum of products. Adding (:= or) all midterms.
				if (!init_AG){
					AG = minterm;
					init_AG = true;
				}
				else{
					AG |= minterm;
				}
			}
		}
	}
	printf("ShortestPath::getCostAdjacencyMatrix: Total iterations: %d\n", total_iter);
	return AG;
}





/**/
bool ShortestPath::createBDDVariables(BDD *system, std::vector<int> vars_index, int no_states, int no_inputs, BDD *x, BDD *u, BDD *x_){

	int i;
	bool status = false;

	printf("ShortestPath::createVariables: Getting the indexes of the variables (Total: %d).\n", vars_index.size());

		for (std::vector<int>::iterator i = vars_index.begin(); i!= vars_index.end(); ++i){
			printf("%d - ", *i);
		}
		printf("\n");
	//	for(i = 0; i < vars_index.size(); i++){
	//		printf("%d - ", vars_index[i]);
	//	}

	// Create the Variables.
	printf("ShortestPath::createVariables: Creating the variables.\n");
	for (i = 0; i < no_states; i++){
//		printf("x: %d - x': %d\n", vars_index[i], vars_index[i + (no_states + no_inputs)]);
		x[i]   = mgr_cpp->bddVar(vars_index[i]);
		x_[i]  = mgr_cpp->bddVar(vars_index[i + (no_states + no_inputs)]);
	}

	for (i = 0; i <  no_inputs; i++){
//		printf("u: %d\n", vars_index[i + no_states]);
		u[i]  = mgr_cpp->bddVar(vars_index[no_states + i]);
	}

	status = true; // TODO: check this. maybe void.
	return status;
}

/**/
bool ShortestPath::createADDVariables(BDD *system, std::vector<int> vars_index, int no_states, int no_inputs, ADD *x, ADD *x_){

	int i;
	bool status = false;

	printf("ShortestPath::createVariables: Getting the indexes of the variables (Total: %d).\n", vars_index.size());

		for (std::vector<int>::iterator i = vars_index.begin(); i!= vars_index.end(); ++i){
			printf("%d - ", *i);
		}
		printf("\n");
	//	for(i = 0; i < vars_index.size(); i++){
	//		printf("%d - ", vars_index[i]);
	//	}

	// Create the Variables.
	printf("ShortestPath::createVariables: Creating the variables.\n");
	for (i = 0; i < no_states; i++){
//		printf("x: %d - x': %d\n", vars_index[i], vars_index[i + (no_states + no_inputs)]);
		x[i]   = mgr_cpp->addVar(vars_index[i]);
		x_[i]  = mgr_cpp->addVar(vars_index[i + (no_states + no_inputs)]);
	}


	status = true; // TODO: check this. maybe void.
	return status;
}

/**/
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

	return vars_index;
}

/**/
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
}


/* Important: This functions assumes that there has not been any variable re-oredering. Variables are in the form of x>u>x'.*/
BDD ShortestPath::BDDTransition(BDD *x, BDD *u, BDD *x_, int no_states, int no_inputs, int xi, int ui, int xi_){

	int i;
	int iteration;

	// Define the iteration
	if (no_inputs > no_states){
		iteration = no_inputs;
	}
	else{
		iteration = no_states;
	}

	// Initialize the BDD containing the transition.
	BDD transition = mgr_cpp->bddOne();

	// Create the transition
	for (i = 0; i < iteration; i++){

		if (i < no_states){
			if (i == xi){
				transition *= x[i] *  (!x_[i]);
			}
			else if (i == xi_){
				transition *= (!x[i]) * x_[i];
			}
			else{
				transition *= (!x[i]) * (!x_[i]);
			}
		}

		if (i < no_inputs){
			if (i == ui){
				transition *= u[i];
			}
			else{
				transition *= !u[i];
			}
		}
	}

	return transition;
}


/**/
bool ShortestPath::FloydWarshall(ADD *AG) {

	printf("\nFloyd Warshall Algorithm.\n");

	// C++ to C
	DdNode *AG_ = AG->getNode();

	DdNode *S, *P;
	DdNode *one, *zero, *xminterm, *yminterm, *temp_node;
	DdNode *R, *C;
	DdNode **xx;
	DdNode **yy;

	int i, j;
	int no_vars;
	int matrix_elements;
	int element;

	DdNode *Result[2];
	DdNode *TR, *TR_temp;



//	long start_time = util_cpu_time();

	// Get the index of the variables of the BDD, representing the system.
//	std::vector<int> vars_index = getVarsIndex(AG); // TODO: this is done twice if we call the getCostAdjacencyMatrix function
	std::vector<int> vars_index;
	vars_index.push_back(0);
	vars_index.push_back(2);
	vars_index.push_back(1);
	vars_index.push_back(3);
	no_vars = vars_index.size()/2;
	printf("Number of variables: %d\n", no_vars);

	// Memory allocation
	xx = (DdNode **) malloc(sizeof(DdNode *) * no_vars);
	yy = (DdNode **) malloc(sizeof(DdNode *) * no_vars);





	// "Copy" the AG matrix.
	S = AG_;
	// Initialize the Pointer array.
	TR = Cudd_addConst(mgr_cpp->getManager(), 0);

	// Zero and One (constant) nodes. Used to create the minterms.
	one  = Cudd_ReadOne(mgr_cpp->getManager());
	zero = Cudd_ReadZero(mgr_cpp->getManager());
	Cudd_Ref(one);
	Cudd_Ref(zero);

	// Create the ADD variables.
	printf("Creating the ADD variables\n");
	// It is important here that they have the same index
	// with the ADD they are going to be used with.
	for (i = 0; i < no_vars; i++) {
		printf("x: %d y: %d \n", i, i +no_vars);
		xx[i] = Cudd_addIthVar(mgr_cpp->getManager(), vars_index[i]);
		yy[i] = Cudd_addIthVar(mgr_cpp->getManager(), vars_index[i + no_vars]);
		Cudd_Ref(xx[i]);
		Cudd_Ref(yy[i]);
    }


	/* Iterate over all matrix elements, i.e. nodes of the initial DD. */
	matrix_elements = 1 << no_vars;
	printf("Iterating over all matrix elements (%d).\n", matrix_elements);

	for (i = 0; i < matrix_elements; i++){


//		printf("\nNode (%d).\n", i);

		element = i;
		xminterm = one;
		yminterm = one;

		/* Creating the minterms. */
		// always construct the ADD from the bottom to the top.
		for (j = no_vars - 1; j >=0; j--){
			if (element & 1){
				// row
				temp_node = Cudd_addIte(mgr_cpp->getManager(), xx[j], xminterm, zero);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr_cpp->getManager(),xminterm);
				xminterm = temp_node;
				// column
				temp_node = Cudd_addIte(mgr_cpp->getManager(), yy[j], yminterm, zero);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr_cpp->getManager(),yminterm);
				yminterm = temp_node;
			}
			else{
				// row
				temp_node = Cudd_addIte(mgr_cpp->getManager(), xx[j], zero, xminterm);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr_cpp->getManager(),xminterm);
				xminterm = temp_node;
				// column
				temp_node = Cudd_addIte(mgr_cpp->getManager(), yy[j], zero, yminterm);
				Cudd_Ref(temp_node);
				Cudd_RecursiveDeref(mgr_cpp->getManager(),yminterm);
				yminterm = temp_node;
			}

			element >>= 1;
		}


		/* Co-factor the matrix. */
		// row
		R = Cudd_Cofactor(mgr_cpp->getManager(), S, xminterm);
		Cudd_Ref(R);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), xminterm);
		// column
		C = Cudd_Cofactor(mgr_cpp->getManager(), S, yminterm);
		Cudd_Ref(C);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), yminterm);

		/* Compute the outer sum. */
//		P = Cudd_addOuterSum(mgr_cpp->getManager(),S,R,C);
		AddOuterSumTrace(S,R,C,Result,i+1);
		P = Result[0];
		Cudd_Ref(P);
		TR_temp = Cudd_addApply(mgr_cpp->getManager(),Cudd_addPlus,TR,Result[1]);
		Cudd_Ref(TR_temp);

		// Keep the new matrix P. Delete others.
		Cudd_RecursiveDeref(mgr_cpp->getManager(),R);
		Cudd_RecursiveDeref(mgr_cpp->getManager(),C);
		Cudd_RecursiveDeref(mgr_cpp->getManager(),S);
		Cudd_RecursiveDeref(mgr_cpp->getManager(),TR);

		S = P;
		TR = TR_temp;

//		char buffer[5] = "0";
//		sprintf(buffer, "%d", i);
//		FILE *dotfile1; //output file pointer for .dot file
//		dotfile1 = fopen(buffer, "w");
//		Cudd_DumpDot(mgr_cpp->getManager(), 1, &TR, NULL, NULL, dotfile1);
//		fclose(dotfile1);
	}

	Cudd_Deref(S); // why? check this.

	// Memory De-allocation.
	free(xx);
	free(yy);

	/* Print execution time. */
//	printf("Floyd Warshall execution time: %s\n", util_print_time(util_cpu_time() - start_time));

	FILE *dotfile; //output file pointer for .dot file
	dotfile = fopen("ADD_FW_cpp.dot", "w");
	Cudd_DumpDot(mgr_cpp->getManager(), 1, &S, NULL, NULL, dotfile);
	fclose(dotfile);

	FILE *dotfile1; //output file pointer for .dot file
	dotfile1 = fopen("ADD_FW_cpp_TR.dot", "w");
	Cudd_DumpDot(mgr_cpp->getManager(), 1, &TR, NULL, NULL, dotfile1);
	fclose(dotfile1);

	return 0;
}


/**Function********************************************************************

  Synopsis    [Takes the minimum of a matrix and the outer sum of two vectors.]

  Description [Takes the pointwise minimum of a matrix and the outer
  sum of two vectors.  This procedure is used in the Floyd-Warshall
  all-pair shortest path algorithm.  Returns a pointer to the result if
  successful; NULL otherwise.]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
void
ShortestPath::AddOuterSumTrace(
  DdNode *M,
  DdNode *r,
  DdNode *c,
  DdNode **Result,
  unsigned int node)
{
    //DdNode *res;

    do {
//		printf("\nCudd_addOuterSumTrace...\n");
		mgr_cpp->getManager()->reordered = 0;
		AddOuterSumRecurTrace(M, r, c, Result, node);
	} while (mgr_cpp->getManager()->reordered == 1);
	//return(res);
//	printf("\nCudd_addOuterSumTrace...end! Returning Result\n");
} /* end of Cudd_addOuterSum */



/**Function********************************************************************

  Synopsis    [Performs the recursive step of Cudd_addOuterSum.]

  Description [Performs the recursive step of Cudd_addOuterSum.
  Returns a pointer to the result if successful; NULL otherwise.]

  SideEffects [None]

  SeeAlso     []

******************************************************************************/
void
ShortestPath::AddOuterSumRecurTrace(
  DdNode *M,
  DdNode *r,
  DdNode *c,
  DdNode **Result,
  unsigned int node)
{

//	printf("cuddAddOuterSumRecurTrace\n");

	DdNode *P, *R, *T, *Mt, *Me, *Tt, *Te, *rt, *re, *ct, *ce, *Rt, *Re;
	int topM, topc, topr;
	int v, index;

	statLine(mgr_cpp->getManager());
	/* Check special cases. */
	if (r == DD_PLUS_INFINITY(mgr_cpp->getManager()) || c == DD_PLUS_INFINITY(mgr_cpp->getManager())) {
//		printf("cuddAddOuterSumRecurTrace: return PLUS INFINITY\n");
		Result[0] = M;
		Result[1] = cuddUniqueConst(mgr_cpp->getManager(), 0);
		return;
	}

	if (cuddIsConstant(c) && cuddIsConstant(r)) {
		R = cuddUniqueConst(mgr_cpp->getManager(), Cudd_V(c) + Cudd_V(r) );
		cuddRef(R);
		if (cuddIsConstant(M)) {
			if (cuddV(R) < cuddV(M)) {
//				printf("cuddAddOuterSumRecurTrace: return constant/			less\n");
				cuddDeref(R);
				Result[0] = R;
				Result[1] = cuddUniqueConst(mgr_cpp->getManager(), node);
				return;
			} else {
//				printf("cuddAddOuterSumRecurTrace: return constant/nf\n");
				Cudd_RecursiveDeref(mgr_cpp->getManager(), R);
				Result[0] = M;
				Result[1] = cuddUniqueConst(mgr_cpp->getManager(), 0);
				return;
			}
		} else {
			printf("\nIN\n");
			P = Cudd_addApply(mgr_cpp->getManager(), Cudd_addMinimum, R, M);
			cuddRef(P);
			Cudd_RecursiveDeref(mgr_cpp->getManager(), R);
			cuddDeref(P);
			Result[0] = P;
			return;
		}
	}

	/* Check the cache. */
	R = cuddCacheLookup(mgr_cpp->getManager(), DD_ADD_OUT_SUM_TAG, M, r, c);
	if (R != NULL) {
		printf("cuddAddOuterSumRecurTrace: return cache\n");
		Result[1] = cuddCacheLookup(mgr_cpp->getManager(), DD_ADD_OUT_SUM_TRACE_TAG, M, r, c);
		Result[0] = R;
		return;
	}

	topM = cuddI(mgr_cpp->getManager(),M->index);
	topr = cuddI(mgr_cpp->getManager(),r->index);
	topc = cuddI(mgr_cpp->getManager(),c->index);
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

//	printf("\ncuddAddOuterSumRecurTrace: ELSE\n");

	AddOuterSumRecurTrace(Me, re, ce, Result, node);
	Re = Result[0];
	Te = Result[1];
	if (Re == NULL) {
		Cudd_RecursiveDeref(mgr_cpp->getManager(), Rt);
		Result[0] = NULL;
		if (Te == NULL)
			Result[1] = NULL;
		return;
	}
	cuddRef(Re);
	cuddRef(Te);

	index = mgr_cpp->getManager()->invperm[v];
	R = (Rt == Re) ? Rt : cuddUniqueInter(mgr_cpp->getManager(), index, Rt, Re);

//	printf("cuddAddOuterSumRecurTrace: cuddUniqueInter - index: %d\n", index);
	T = (Tt == Te) ? Tt : cuddUniqueInter(mgr_cpp->getManager(), index, Tt, Te);
	//if (topM == v) {      T = (Tt == Te) ? Tt : Cudd_addIte(mgr_cpp->getManager(),M,Tt,Te); }
	//else if (topr == v) { T = (Tt == Te) ? Tt : Cudd_addIte(mgr_cpp->getManager(),r,Tt,Te); }
	//else if (topc == v) { T = (Tt == Te) ? Tt : Cudd_addIte(mgr_cpp->getManager(),c,Tt,Te);}

	if (R == NULL) {
		Cudd_RecursiveDeref(mgr_cpp->getManager(), Rt);
		Cudd_RecursiveDeref(mgr_cpp->getManager(), Re);
		Result[0] = NULL;

		if (T == NULL) {
			Cudd_RecursiveDeref(mgr_cpp->getManager(), Tt);
			Cudd_RecursiveDeref(mgr_cpp->getManager(), Te);
			Result[1] = NULL;
		}

		return;
	}
	cuddDeref(Rt);
	cuddDeref(Re);
	cuddDeref(Tt);
	cuddDeref(Te);

	/* Store the result in the cache. */
	cuddCacheInsert(mgr_cpp->getManager(), DD_ADD_OUT_SUM_TAG, M, r, c, Result[0]);
	cuddCacheInsert(mgr_cpp->getManager(), DD_ADD_OUT_SUM_TRACE_TAG, M, r, c, T);

//	printf("cuddAddOuterSumRecurTrace: return final\n");
	Result[0] = R;
	Result[1] = T;

} /* end of cuddAddOuterSumRecur */
