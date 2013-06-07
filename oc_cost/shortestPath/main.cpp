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
 * Contains main function used to test the ShortestPath Class.
 *
 * 		No details yet.
 */

#include "main.hh"

#define U_SPACING 100
#define X__SPACING 1000


int main(int argc, char* argv[]) {

	printf("Main.\n\n");

	example_DSP();

	printf("\n\nExiting Program...\n");

	return 0;
}





void example_DSP(){

	printf("***Deterministic Shortest Path Example.***\n\n");

	int no_states, no_inputs;
	int *state_costs;
	std::vector<int> target_set;

	// Initialize the Manager.
	Cudd mgr(0, 0);

	// Shortest Path Object
//	ShortestPath sp(&mgr);

	// BDD of the System
	BDD S;
	// ADD of state costs
	ADD C;
	// ADD representing the Cost Adjacency Matrix
	ADD AG;
	// ADD of the Target Set
	ADD W;




	/************************ Define the System ******************************/

	// Number of states/inputs
	no_states = 5;
	no_inputs = 5;
//	no_states = 23;
//	no_inputs = 36;


	state_costs = new int[no_states];

	// Give the cost of each state
//	for (int i = 0; i < no_states; i++){
//		state_costs[i] = i + 1;
//	}
	state_costs[0] = 1;
	state_costs[1] = 2;
	state_costs[2] = 1;
	state_costs[3] = 3;
	state_costs[4] = 5;

	// Define the System in terms of transitions.
	// Important: Only deterministic transitions!
	#define SYSTEM_TRANSITIONS \
							\
	  TRANSITION(0,0,1)		\
	+ TRANSITION(0,2,2)		\
	+ TRANSITION(1,3,3)     \
	+ TRANSITION(2,2,3)		\
	+ TRANSITION(1,4,4)		\
	+ TRANSITION(2,4,4)

	/*
	+ TRANSITION(3,4,4)		\
	+ TRANSITION(4,4,1)		\
	+ TRANSITION(5,4,4)		\
	+ TRANSITION(6,4,5)		\
	+ TRANSITION(7,4,4)		\
	+ TRANSITION(8,14,4)	\
	+ TRANSITION(9,20,4)	\
	+ TRANSITION(10,4,4)	\
	+ TRANSITION(11,4,6)	\
	+ TRANSITION(12,4,4)	\
	+ TRANSITION(13,4,8)	\
	+ TRANSITION(14,4,4)	\
	+ TRANSITION(15,35,2)	\
	+ TRANSITION(16,4,14)	\
	+ TRANSITION(17,4,10)	\
	+ TRANSITION(18,4,17)	\
	+ TRANSITION(19,29,7)	\
	+ TRANSITION(20,30,20)	\
	+ TRANSITION(21,15,17)	\
	+ TRANSITION(22,26,11)
	*/

	// Give the Target Set W.
	target_set.push_back(3);
	target_set.push_back(4);

	/************************************************************************/


	/* Get the system as a BDD that satisfies (x,u,x') -> {0,1}. */
	get_S_xux(&mgr, &S, no_states, no_inputs);

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	nodes_bdd.push_back(S);
	FILE *outfile;
	outfile = fopen("System_xux.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

	/* Get the cost of each transition of the system, described by an ADD. x -> R. */
	get_S_cost_x(&mgr, &C, no_states, state_costs);

	// Create .dot file
	std::vector<ADD> nodes_add;
	nodes_add.push_back(C);
	outfile = fopen("System_Cost_x.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();


	long long start_time = get_usec();

	/* Create the Shortest Path Object */
	ShortestPath sp(&mgr, &S, no_states, no_inputs); // optimized
//	ShortestPath sp(&mgr);

	/* Create the Cost Adjacency Matrix */
	AG = sp.createCostAdjacencyMatrix(&S, &C, no_states, no_inputs);

	// Create .dot file
	nodes_add.push_back(AG);
	outfile = fopen("System_CostAdjMatrix.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find All-pair Shortest Path */
	ADD APSP;
	ADD PA;
	sp.FloydWarshall(&AG, &APSP, &PA);

	// Create .dot file
	nodes_add.push_back(APSP);
	outfile = fopen("System_APSP.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();
	// Create .dot file
	nodes_add.push_back(PA);
	outfile = fopen("System_APSP_PA.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find the All-pair Shortest Path to  a Set. */
	ADD APSP_W;
	ADD PA_W;
	sp.APtoSetSP(&APSP, &PA, target_set, &APSP_W, &PA_W);

	// Create .dot file
	nodes_add.push_back(APSP_W);
	outfile = fopen("System_APSP_W.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();
	// Create .dot file
	nodes_add.push_back(PA_W);
	outfile = fopen("System_APSP_PA_W.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();


	printf("\n***Deterministic Shortest Path Example END. Execution time: %ds (%dms) (%dus)***\n", (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
}


/**/
void get_S_xux(Cudd *mgr, BDD *T, int no_states, int no_inputs){


	BDD *x;
	BDD *u;
	BDD *x_;

	int i;

	int no_state_vars = no_states/2 + no_states % 2;
	int no_input_vars = no_inputs/2 + no_inputs % 2;

	x  = new BDD[no_state_vars];
	u  = new BDD[no_input_vars];
	x_ = new BDD[no_state_vars];

	// Create the variables
	for (i = 0; i < no_state_vars; i++){
		x[i]  = mgr->bddVar(i);
		x_[i] = mgr->bddVar(i + X__SPACING);
	}

	for (i = 0; i < no_input_vars; i++){
		u[i]  = mgr->bddVar(i + U_SPACING);
	}

	// Create the system
	*T = SYSTEM_TRANSITIONS;
}

/**/
void get_S_cost_x(Cudd *mgr, ADD *C, int no_states, int *costs){

	ADD *x_, *constants;

	int no_state_vars = no_states/2 + no_states % 2;
	int i, j;

	ADD one  = mgr->addOne();
	ADD zero = mgr->addZero();

	x_        = new ADD[no_state_vars];
	constants = new ADD[no_states];


	// Set background (if not set) :P
	mgr->SetBackground(mgr->plusInfinity());


	for (i = 0; i < no_state_vars; i++){
		x_[i]  = mgr->addVar(i + X__SPACING);
		x_[i]  = x_[i].Ite(one, zero);
	}


	// Get the costs
	for(i = 0; i < no_states; i++){
		constants[i] = mgr->constant(costs[i]);
	}


	// Auxiliary Variables
	int state;
	ADD temp;
	ADD minterm;


	// Create the cost matrix
	*C = mgr->background();

	for (i = 0; i < no_states; i++){

		state = i;
		minterm = one;

		for(j = no_state_vars - 1; j >= 0; j--){

			if (state & 0x01){
				temp = minterm * x_[j];
			}
			else{
				temp = minterm * (~x_[j]);
			}
			minterm = temp;
			state >>= 1;
		}

		// Create the constant node.
		temp = minterm.Ite(constants[i], *C);
		*C = temp;
	}
}

/**/
BDD BDD_transition(Cudd *mgr, BDD *x, BDD *u, BDD *x_, int no_state_var, int no_input_var, int xi, int ui, int xi_){

	int i;

	int iteration;

	if (no_input_var > no_state_var){
		iteration = no_input_var;
	}
	else{
		iteration = no_state_var;
	}

	BDD transition = mgr->bddOne();


	for (i = 0; i < iteration; i++){

		// Writing backwards because of the format: MSB--LSB
		if (i < no_state_var){
			if ((xi & 0x01) == 1){
				transition *= x[no_state_var -1 - i];
			}
			else{
				transition *= !x[no_state_var -1 - i];
			}

			if ((xi_ & 0x01) == 1){
				transition *= x_[no_state_var -1 - i];
			}
			else{
				transition *= !x_[no_state_var -1 - i];
			}
		}


		if (i < no_input_var){
			if ((ui & 0x01) == 1){
				transition *= u[no_input_var -1 - i];
			}
			else{
				transition *= !u[no_input_var -1 - i];
			}
		}

		//
		xi  >>= 1;
		ui  >>= 1;
		xi_ >>= 1;
	}

	return transition;
}


long long get_usec(void){
	long long r;
	struct timeval t;


	gettimeofday(&t, NULL);
	r = t.tv_sec * 1000000 + t.tv_usec;
	return r;
}
























