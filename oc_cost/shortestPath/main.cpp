/*
 * main.cpp
 *
 *  Created on: Apr 13, 2013
 *      Author: tanasaki
 */

#include "main.hh"


int main(int argc, char* argv[]) {

	printf("Main.\n\n");


	example_FW();

//	example_DSP();

	printf("\n\nExiting Program...\n");

	return 0;
}





void example_DSP(){

	printf("Deterministic Shortest Path Example.\n\n");

	int no_states, no_inputs;

	// Initialize the Manager.
	Cudd mgr(0, 0);

	// Shortest Path Object
	ShortestPath sp(&mgr);

	// BDD of the System
	BDD S;
	// ADD of state costs
	ADD C;
	// ADD representing the Cost Adjacency Matrix
	ADD AG;


	/* Get the system as a BDD that satisfies (x,u,x') -> {0,1}. */
	get_S_xux(&mgr, &S);

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	nodes_bdd.push_back(S);
	FILE *outfile;
	outfile = fopen("System_xux.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

//	mgr.ReduceHeap(CUDD_REORDER_RANDOM, 5);

	/* Get the cost of each transition of the system, described by an ADD. x -> R. */
	get_S_cost_x(&mgr, &C);

	// Create .dot file
	std::vector<ADD> nodes_add;
	nodes_add.push_back(C);
	outfile = fopen("System_Cost_x.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Create the Cost Adjacency Matrix */
	no_states = 4;
	no_inputs = 6;
	AG = sp.getCostAdjacencyMatrix(&S, &C, no_states, no_inputs);

	// Create .dot file
	nodes_add.push_back(AG);
	outfile = fopen("System_CostAdjMatrix.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find All-pair Shortest Path */
	ADD APSP;
	ADD P;
	sp.FloydWarshall(&AG);













//	BDD Z;
//
//	BDD *x;
//	BDD *u;
//	BDD *x_;
//
//	int i;
//
//	x  = new BDD[no_states];
//	u  = new BDD[no_inputs];
//	x_ = new BDD[no_states];
//
//
//	for (i = 0; i < no_states; i++){
//		x[i]  = mgr.bddVar(i);
//		x_[i] = mgr.bddVar(i + 100);
//	}
//
//	for (i = 0; i < no_inputs; i++){
//		u[i]  = mgr.bddVar(i + 10);
//	}
//
//	Z = BDD_transition(&mgr, x, u, x_, no_states, no_inputs, 3,1,3);
//
//	// Create .dot file
//	nodes_bdd.push_back(Z);
//	outfile = fopen("System_input_expl.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//
//
//	Z = S.Restrict(x[0]);
//	Z = Z.Restrict(u[4]);
//
//	Z = S.Restrict(Z);
//	Z = mgr.VectorSupport(nodes_bdd);
//
//	if (Z.IsOne()){
//		printf("Valid input.\n");
//	}
//	if (Z.IsZero()){
//		printf("Invalid input.\n");
//	}
//
//
//
//	Z = BDD_transition(&mgr, x, u, x_, no_states, no_inputs, 3,2,0);
//	Z = S.Restrict(Z);
//
//	if (Z.IsOne()){
//		printf("Valid input.\n");
//	}
//	if (Z.IsZero()){
//		printf("Invalid input.\n");
//	}




	// Create .dot file
//	nodes_bdd.push_back(Z);
//	outfile = fopen("System_output_expl.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();




}


/**/
void get_S_xux(Cudd *mgr, BDD *T){


	BDD *x;
	BDD *u;
	BDD *x_;

	int i;

	int no_states = 4;
	int no_inputs = 6;

	x  = new BDD[no_states];
	u  = new BDD[no_inputs];
	x_ = new BDD[no_states];

//	BDD one  = mgr->bddOne();
//	BDD zero = mgr->bddZero();
//	BDD zero = !one;


	for (i = 0; i < no_states; i++){
		x[i]  = mgr->bddVar(i);
		x_[i] = mgr->bddVar(i + 100);
	}

	for (i = 0; i < no_inputs; i++){
		u[i]  = mgr->bddVar(i + 10);
	}


	/*************** Create the example-system ******************/

	*T =
		  BDD_transition(mgr,x,u,x_,no_states,no_inputs,0,0,1)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,0,1,2)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,0,2,3)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,1,3,2)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,2,4,1)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,3,1,2)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,3,2,0)
		+ BDD_transition(mgr,x,u,x_,no_states,no_inputs,3,5,1)
		;

	/***********************************************************/

	/* Create .dot file */
//	std::vector<BDD> nodes;
//	nodes.push_back(*T);
//	FILE *outfile;
//	outfile = fopen("System_T_xux.dot", "w");
//	mgr->DumpDot(nodes, NULL, NULL, outfile);
//	fclose(outfile);
}

/**/
void get_S_cost_x(Cudd *mgr, ADD *C){

	ADD *x_, *constants;
	int no_states = 4;
	int i;


	x_        = new ADD[no_states];
	constants = new ADD[no_states];


	// set background
	mgr->SetBackground(mgr->plusInfinity());

	for (int i = 0; i < no_states; i++){
		x_[i] = mgr->addVar(i + 100);
	}


	/************** Give the cost of each state  ****************/

	constants[0] = mgr->constant(2);
	constants[1] = mgr->constant(2);
	constants[2] = mgr->constant(4);
	constants[3] = mgr->constant(3);

	/************************************************************/

	x_[no_states - 1] = x_[no_states - 1].Ite(constants[no_states -1], mgr->plusInfinity());

	for (i = (no_states-2); i>= 0; i--){
		x_[i] = x_[i].Ite(constants[i], x_[i+1]);
	}

	*C = x_[0];

	/* Create .dot file */
//	std::vector<ADD> nodes;
//	nodes.push_back(*C);
//	FILE *outfile;
//	outfile = fopen("System_Cost_x.dot", "w");
//	mgr->DumpDot(nodes, NULL, NULL, outfile);
//	fclose(outfile);
}

/**/
BDD BDD_transition(Cudd *mgr, BDD *x, BDD *u, BDD *x_, int no_states, int no_inputs, int xi, int ui, int xi_){

	int i;

	int iteration;

	if (no_inputs > no_states){
		iteration = no_inputs;
	}
	else{
		iteration = no_states;
	}

	BDD transition = mgr->bddOne();

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

	//printf("Transition end.\n");
	return transition;
}


void example_FW(){

		ADD E[1];
		ADD** x;
		ADD** y;
		ADD** xn;
		ADD** yn;
		int * nx;
		int * ny;
		int * m;
		int * n;



		// Memory allocations
		nx = (int *) malloc(sizeof(int));
		ny = (int *) malloc(sizeof(int));
		m  = (int *) malloc(sizeof(int));
		n  = (int *) malloc(sizeof(int));

		// Initialize the Manager.
		Cudd mgr;

		// Set new background. (We need plus-infinity for the all-pair shortest path problem).
		mgr.SetBackground(mgr.plusInfinity());

		// Open the file (read mode).
		FILE* matrixfile = fopen("matrix3.txt", "r");
		if (matrixfile == NULL) {
			printf("Error while opening file.\n");
			exit(1);
		}


		// Import the matrix from the file.

		/*
		 * Caution!!
		 *
		 * The variable bx and by have to have a reasonable spacing between them. They are used to index to row and column variables
		 * respectively. The can be chosen arbitrary but with caution. For example if if bx = 0 and bx =1  (and step 1, i.e. sx=sy=1)
		 * and we have more the one row or column variable, then the row variable is going to "override" the column  variable and they
		 * will be seen as the same and therefore the ADD representation is going to be wrong.
		 *
		 * So for example if we need 2 variables to decode the rows and by using bx = 0 and (step) sx = 1 then
		 * the first row variable is going to be indexed as 0 and the second as 1. This "forbids" us to choose
		 * by = 0 or 1. So use for instance by = 100 with (step) sy = 1 to get column variables indexed as 100 and 101.
		 *
		 * */


		printf("Reading from file... ");

		mgr.Read(matrixfile, E, x, y, xn, yn, nx, ny, m, n, 0, 1, 2, 1);

		std::vector<ADD> nodes;
		nodes.push_back(E[0]);

		// Reorder the ADD.
		int order[4] = { 0, 2, 1, 3 };
		//	printf("%d %d %d %d \n", order[0], order[1], order[2], order[3]);
		mgr.ShuffleHeap(order);

		FILE *outfile; //output file pointer for .dot file
		outfile = fopen("A_matrix_cpp.dot", "w");
		mgr.DumpDot(nodes, NULL, NULL, outfile);
		fclose(outfile);



		/* All-pair shortest path example. */

		printf("\nAll-pair shortest path example.\n");
		FloydWarshallParam *param = new FloydWarshallParam(*nx);


		param->x_index[0] = 0;
		param->x_index[1] = 1;
		param->y_index[0] = 2;
		param->y_index[1] = 3;



		/*
		 *
		 * Floyd Warshall
		 *
		 * */
		ShortestPath sp(&mgr);
		sp.FloydWarshall(&E[0]);
		mgr.FloydWarshall(&E[0], param);

		printf("\nFW-Algorithm end.\n");

	//
	//	/*
	//	 *
	//	 * Repeated Squaring
	//	 *
	//	 *
	//	 */
	//
	//
	//
	//	/* De-allocate Memory. */
	//	printf("Main: De-allocating memory. \n");
	//	delete param;
	//	free(add_matrix);
	//
	//
	//	free(xn[0]);
	//	free(yn[0]);
	//	free(x[0]);
	//	free(y[0]);
	//	free(x);
	//	free(y);
	//	free(xn);
	//	free(yn);
	//	free(nx);
	//	free(ny);
	//	free(m);
	//	free(n);

}
























