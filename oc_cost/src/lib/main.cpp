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
 * Contains main function used to test the ShortestPath Class.
 *
 * 		No details yet.
 */

#include "main.hh"

// Define if you want some space between the variables.
#define U_SPACING  100
#define X__SPACING 1000

//
#ifdef U_SPACING
	#if (U_SPACING > 0)
		#define HAS_U_SPACING
	#endif
#endif
//
#ifdef X__SPACING
	#if (X__SPACING > 0)
		#define HAS_X__SPACING
	#endif
#endif


///////////////////////////////////////////////////////////
/* MAIN */
int main(int argc, char* argv[]) {

	printf("Main.\n\n");

//	example_DSP();
	example_NDSP();
//	test_actual();

	printf("\n\nExiting Program...\n");

	return 0;
}
///////////////////////////////////////////////////////////


/**/
void test_actual(){


	unsigned int nstates, ninputs;
	int i;

	char *bddsys       = "DCMotor";
	char *bddsys_costs = "DCMotorCosts";
	char *bddsys_tset  = "DCMotorTargetSet";
	char *bddsys_name;
	char *bddsys_costs_name;
	char *bddsys_tset_name;
	char bdd_suffix[5];
	char add_suffix[5];



	nstates = 2583;
	ninputs = 2001;

	// Print some extra-info
	printf("\nSymbolic model size: ");
	printf("%u", nstates);
	printf(" states; \n");
	printf("                     %u", ninputs);
	printf(" inputs. \n");


	/* Check if all necessary files exist! */

	// assign the .bdd/.add file suffix
	strcpy(bdd_suffix,".bdd");
	strcpy(add_suffix,".add");

	/* Check if the System's BDD exists. */
	bddsys_name=(char*)malloc(strlen(bddsys)+5);
	strcpy(bddsys_name, bddsys);
	strcat(bddsys_name, bdd_suffix);
//	mxFree(bddsys);
	//check for file existence
	printf("Checking %s ...\n", bddsys_name);
	FILE * smFile;
	smFile = fopen(bddsys_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the System's-States-Cost ADD file exists. */
	bddsys_costs_name=(char*)malloc(strlen(bddsys_costs)+5);
	strcpy(bddsys_costs_name, bddsys_costs);
	strcat(bddsys_costs_name, add_suffix);
//	mxFree(bddsys);
	//check for file existence
	printf("Checking %s ...\n", bddsys_costs_name);
	smFile = fopen(bddsys_costs_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the Target Set BDD file exists. */
	bddsys_tset_name=(char*)malloc(strlen(bddsys_tset)+5);
	strcpy(bddsys_tset_name, bddsys_tset);
	strcat(bddsys_tset_name, bdd_suffix);
//	mxFree(bddsys);

	//check for file existence
	printf("Checking %s ...\n", bddsys_tset_name);
	smFile = fopen(bddsys_tset_name,"r");
	FILE_EXISTS(smFile)


	/* Files exists... now load them. */

	// Initialize the Manager.
	Cudd mgr(0, 0); //getNoBits(nstates)
	// Set background
	mgr.SetBackground(mgr.plusInfinity());
	//
//	mgr.DisableGarbageCollection();

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	std::vector<ADD> nodes_add;
	FILE *outfile;

//	BDD S  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, "DCMotor.bdd", NULL));

	BDD CS        = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, "DCMotorController.bdd", NULL));
	BDD Contr_dom = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, "DCMotorController_dom.bdd", NULL));

//	// Create .dot file
//	nodes_bdd.push_back(S);
//	outfile = fopen("DCMotor.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();

//	// Create .dot file
//	nodes_bdd.push_back(Contr_dom);
//	outfile = fopen("DCMotorController_dom.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();

	ADD SC = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, "DCMotorCosts.add", NULL));

	// Create .dot file
	nodes_add.push_back(SC);
	outfile = fopen("DCMotorCosts.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	BDD W  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, "DCMotorTargetSet.bdd", NULL));

	// Create .dot file
	nodes_bdd.push_back(W);
	outfile = fopen("DCMotorTargetSet.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();


	/* Create the Shortest Path Object */
	printf("Creating SP Object ...\n");
	ShortestPath sp(&mgr, &CS, nstates, ninputs); // optimized
//	ShortestPath sp(&mgr);

//	/* */
//	printf("Checking Controller domain...\n");
//	sp.checkControllerDom(&CS, &Contr_dom);

	/* Create the Cost Adjacency Matrix */
	printf("Creating Cost Adjacency Matrix ...\n");
	ADD AG;
	AG = sp.createCostAdjacencyMatrix(&CS, &SC, nstates, ninputs);

	printf(" ***Checking createCostAdjacencyMatrix...***\n");
	mgr.DebugCheck();


	// Create .dot file
	nodes_add.push_back(AG);
	outfile = fopen("DCMotor_CostAdjMatrix.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find All-pair Shortest Path */
	printf("Finding All-pair Shortest Path (FW) ...\n");
	ADD APSP;
	ADD PA;
	sp.FloydWarshall(&AG, &APSP, &PA);

	printf(" ***Checking FloydWarshall...***\n");
//	mgr.CheckKeys();
	mgr.DebugCheck();



	// Create .dot file
	nodes_add.push_back(APSP);
	outfile = fopen("DCMotor_APSP.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();
	// Create .dot file
	nodes_add.push_back(PA);
	outfile = fopen("DCMotor_APSP_PA.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find the All-pair Shortest Path to  a Set. */
	printf("Finding All-pair Shortest Path to a Set W ...\n");
	ADD APSP_W;
	ADD PA_W;
	sp.APtoSetSP(&APSP, &PA, &W, &APSP_W, &PA_W);

	printf(" ***Checking APtoSetSP... ***\n");
	mgr.DebugCheck();


	// Create .dot file
	nodes_add.push_back(APSP_W);
	outfile = fopen("DCMotor_APSP_W.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();
	// Create .dot file
	nodes_add.push_back(PA_W);
	outfile = fopen("DCMotor_APSP_PA_W.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();


	/* Create the Controller. */
	BDD controller = sp.createControllerBDD(&CS, &PA, &PA_W);

	printf(" ***Checking createControllerBDD... ***\n");
	mgr.DebugCheck();

	sp.Dddmp_cuddStore(&controller, "DCMotorController_SP.bdd");

	printf("Actual test END\n");
}


/**/
void example_DSP(){

	printf("***Deterministic Shortest Path Example.***\n\n");

	int no_states, no_inputs;
	int *state_costs;
	std::vector<int> target_set;

	// Initialize the Manager.
	Cudd mgr(0, 0);
	// Set background
	mgr.SetBackground(mgr.plusInfinity());

	// Shortest Path Object
//	ShortestPath sp(&mgr);

	// BDD of the System
	BDD S;
	// ADD of state costs
	ADD C;
	// ADD representing the Cost Adjacency Matrix
	ADD AG;
	// ADD of the Target Set
	BDD W;




	/************************ Define the System ******************************/

	// Number of states/inputs
	no_states = 9;
	no_inputs = 9;

	state_costs = new int[no_states];

	// Give the cost of each state
//	state_costs[0] = 1;
//	state_costs[1] = 2;
//	state_costs[2] = 1;
//	state_costs[3] = 3;
//	state_costs[4] = 5;
//	state_costs[5] = 6;

	state_costs[0] = 1;
	state_costs[1] = 1;
	state_costs[2] = 2;
	state_costs[3] = 3;
	state_costs[4] = 1;
	state_costs[5] = 1;
	state_costs[6] = 3;
	state_costs[7] = 7;
	state_costs[8] = 7;

	// Define the System in terms of transitions. (x,u,x')
	// Important: Only deterministic transitions!

	/*

	#define SYSTEM_TRANSITIONS \
		\
	  TRANSITION(0,0,1)		\
	+ TRANSITION(0,2,2)		\
	+ TRANSITION(1,3,3)     \
	+ TRANSITION(1,2,4)		\
	+ TRANSITION(2,2,3)		\
	+ TRANSITION(2,4,4)		\
	+ TRANSITION(3,0,5)     \
	+ TRANSITION(3,1,6)		\
	+ TRANSITION(4,3,5)		\
	+ TRANSITION(4,1,6)		\
	+ TRANSITION(5,0,5)     \
	+ TRANSITION(6,1,5)     \
	+ TRANSITION(6,2,6)		\
	+ TRANSITION(1,1,1)     \
	+ TRANSITION(3,7,4)     \
	+ TRANSITION(4,7,3)     \
	+ TRANSITION(3,8,1)     \
	+ TRANSITION(1,8,0)     \
	+ TRANSITION(0,7,5)		\
	+ TRANSITION(3,6,0)		\
	+ TRANSITION(7,6,8)		\
	+ TRANSITION(8,6,7)		\
	+ TRANSITION(5,7,0)

*/



	/*
	 						\
	  TRANSITION(0,0,1)		\
	+ TRANSITION(0,2,2)		\
	+ TRANSITION(1,3,3)     \
	+ TRANSITION(2,2,3)		\
	+ TRANSITION(1,4,4)		\
	+ TRANSITION(2,4,4)		\
	+ TRANSITION(4,1,5)

	 */

	// Give the Target Set W.
	target_set.push_back(0);
	target_set.push_back(5);
//	target_set.push_back(6);

	/************************************************************************/


	/* Get the system as a BDD that satisfies (x,u,x') -> {0,1}. */
	get_S_xux(&mgr, &S, no_states, no_inputs);

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	FILE *outfile;
	nodes_bdd.push_back(S);
	outfile = fopen("System_xux.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

	/* Get the cost of each transition of the system, described by an ADD. x -> R. */
	get_S_cost_x(&mgr, &C, no_states, no_inputs, state_costs);

	// Create .dot file
	std::vector<ADD> nodes_add;
	nodes_add.push_back(C);
	outfile = fopen("System_Cost_x.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

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


	/* Create the Target set W. */
	W = getTargetSet(&mgr, no_states, target_set);
	BDD W_normal = W;

	// Create .dot file
	nodes_bdd.push_back(W);
	outfile = fopen("System_TargetSet_W.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();


	/* Find the All-pair Shortest Path to  a Set. */
	ADD APSP_W;
	ADD PA_W;
	sp.APtoSetSP(&APSP, &PA, &W, &APSP_W, &PA_W);
//	sp.APtoSetSP(&APSP, &PA, target_set, &APSP_W, &PA_W);

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

	/* Create the Controller. */
	BDD controller = sp.createControllerBDD(&S, &PA, &PA_W);
	nodes_bdd.push_back(controller);
	outfile = fopen("System_Controller.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

//	/* Storing results into .bdd files. */
//	printf("Storing outcome into .bdd and .add files... ");
//	bool stored;
//
//	// BDD's
//	stored = sp.Dddmp_cuddStore(&S, (char *)"System_xux.bdd");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&W_normal, (char *)"System_TargetSet_W.bdd");
//	if(!stored) printf("Error saving ADD into file!\n");
//	// ADD's
//	stored = sp.Dddmp_cuddStore(&C, (char *)"System_Cost_x.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&APSP, (char *)"APSP.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&APSP_W, (char *)"APSP_W.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&PA, (char *)"PA.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&PA_W, (char *)"PA_W.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//
//	stored = writeSysInfo_TBFile(no_states - 1, no_inputs - 1);
//	if(!stored) printf("Error saving System info to .txt!\n");
//
//	printf("done!\n");


	// Memory De-allocation
	delete[] state_costs;

#ifdef ENABLE_TIME_PROFILING
	printf("\n***Deterministic Shortest Path Example END. Execution time: %ds (%dms) (%dus)***\n", (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#else
	printf("\n***Deterministic Shortest Path Example END.***\n");
#endif
}

/**/
void example_NDSP(){

	printf("***Non-Deterministic Shortest Path Example.***\n\n");

	int no_states, no_inputs;
	int *state_costs;
	std::vector<int> target_set;

	// Initialize the Manager.
	Cudd mgr(0, 0);
	// Set background
	mgr.SetBackground(mgr.plusInfinity());

	// BDD of the System
	BDD S;
	// ADD of state costs
	ADD SC;
	// ADD representing the Cost Adjacency Matrix
	ADD AG;
	// ADD of the Target Set
	BDD W;




	/************************ Define the System ******************************/

	// Number of states/inputs
	no_states = 7;
	no_inputs = 4;

	state_costs = new int[no_states];

	// Give the cost of each state
	state_costs[0] = 5;
	state_costs[1] = 2;
	state_costs[2] = 4;
	state_costs[3] = 1;
	state_costs[4] = 1;
	state_costs[5] = 10;
	state_costs[6] = 7;



	// Define the System in terms of transitions.
	#define SYSTEM_TRANSITIONS \
							\
	  TRANSITION(0,2,1)		\
	+ TRANSITION(0,1,6)		\
	  	  	  	  	  	  	\
	+ TRANSITION(1,1,3)		\
	  	  	  	  	  	    \
	+ TRANSITION(2,0,1)		\
	+ TRANSITION(2,0,3)		\
	+ TRANSITION(2,1,5)		\
	+ TRANSITION(2,2,2)		\
	  	  	  	  	  	    \
	+ TRANSITION(3,3,4)		\
	  	  	  	  	  	    \
	+ TRANSITION(4,0,5)		\
	  	  	  	  	  	    \
	+ TRANSITION(5,1,4)		\
	  	  	  	  	  	  	\
	+ TRANSITION(6,0,0)		\
	+ TRANSITION(6,0,2)		\
							\
	+ TRANSITION(2,3,4)     \
	+ TRANSITION(2,0,4)     \
	+ TRANSITION(2,1,4)     \
	+ TRANSITION(4,0,3)     \

	// Give the Target Set W.
	target_set.push_back(4);
	target_set.push_back(5);


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
	get_S_cost_x(&mgr, &SC, no_states, no_inputs, state_costs);

	// Create .dot file
	std::vector<ADD> nodes_add;
	nodes_add.push_back(SC);
	outfile = fopen("System_Cost_x.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Create the Target set W. */
	W = getTargetSet(&mgr, no_states, target_set);
	BDD W_normal = W;

	// Create .dot file
	nodes_bdd.push_back(W);
	outfile = fopen("System_TargetSet_W.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

#ifdef ENABLE_TIME_PROFILING
	long long start_time = get_usec();
#endif

	/* Create the Shortest Path Object */
	ShortestPath sp(&mgr, &S, no_states, no_inputs); // optimized
//	ShortestPath sp(&mgr);


	/* Find the All-pair Shortest Path to  a Set. */
	ADD APSP_W;
	BDD PA_W;
	sp.APtoSetSP(&S, &SC, &W, &APSP_W, &PA_W, no_states, no_inputs);

//	ADD zzz = (~W).Add() * mgr.background();
//
//	// Create .dot file
//	nodes_add.push_back(zzz);
//	outfile = fopen("zzz.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();


//	// Create .dot file
//	nodes_add.push_back(APSP_W);
//	outfile = fopen("System_APSP_W.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//	// Create .dot file
//	nodes_add.push_back(PA_W);
//	outfile = fopen("System_APSP_PA_W.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();




//	/* Storing results into .bdd files. */
//	printf("Storing outcome into .bdd and .add files... ");
//	bool stored;
//
//	// BDD's
//	stored = sp.Dddmp_cuddStore(&S, (char *)"System_xux.bdd");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&W_normal, (char *)"System_TargetSet_W.bdd");
//	if(!stored) printf("Error saving ADD into file!\n");
//	// ADD's
//	stored = sp.Dddmp_cuddStore(&C, (char *)"System_Cost_x.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&APSP, (char *)"APSP.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&APSP_W, (char *)"APSP_W.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&PA, (char *)"PA.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//	stored = sp.Dddmp_cuddStore(&PA_W, (char *)"PA_W.add");
//	if(!stored) printf("Error saving ADD into file!\n");
//
//	stored = writeSysInfo_TBFile(no_states - 1, no_inputs - 1);
//	if(!stored) printf("Error saving System info to .txt!\n");
//	printf("done!\n");


	// Memory De-allocation
	delete[] state_costs;

#ifdef ENABLE_TIME_PROFILING
	printf("\n***Non-Deterministic Shortest Path Example END. Execution time: %ds (%dms) (%dus)***\n", (int)(get_usec() - start_time)/1000000, (int)(get_usec() - start_time)/1000, (int)(get_usec() - start_time));
#else
	printf("\n***Non-Deterministic Shortest Path Example END.***\n");
#endif
}


/**/
void get_S_xux(Cudd *mgr, BDD *T, int no_states, int no_inputs){


	BDD *x;
	BDD *u;
	BDD *x_;

	int i;

	int no_state_vars = getNoBits(no_states - 1);
	int no_input_vars = getNoBits(no_inputs - 1);

	x  = new BDD[no_state_vars];
	u  = new BDD[no_input_vars];
	x_ = new BDD[no_state_vars];

	// Create the variables
#ifdef HAS_X__SPACING
	for (i = 0; i < no_state_vars; i++){
		x[i]  = mgr->bddVar(i);
		x_[i] = mgr->bddVar(i + X__SPACING);
	}
#else
	for (i = 0; i < no_state_vars; i++){
		x[i]  = mgr->bddVar(i);
		x_[i] = mgr->bddVar(i + (no_state_vars + no_input_vars));
	}
#endif

#ifdef HAS_U_SPACING
	for (i = 0; i < no_input_vars; i++){
		u[i]  = mgr->bddVar(i + U_SPACING);
	}
#else
	for (i = 0; i < no_input_vars; i++){
		u[i]  = mgr->bddVar(i + no_state_vars);
	}
#endif

	// Create the system
	*T = SYSTEM_TRANSITIONS;

	// Memory De-allocation
	delete[] x;
	delete[] u;
	delete[] x_;
}

/**/
void get_S_cost_x(Cudd *mgr, ADD *C, int no_states, int no_inputs, int *costs){

	ADD *x_, *constants;

	int no_state_vars = getNoBits(no_states - 1);
	int i, j;

	ADD one  = mgr->addOne();
	ADD zero = mgr->addZero();

	x_        = new ADD[no_state_vars];
	constants = new ADD[no_states];

	for (i = 0; i < no_state_vars; i++){
		x_[i]  = mgr->addVar(i);
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

		for(j = 0; j < no_state_vars; j++){

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

	// Memory De-allocation
	delete[] x_;
	delete[] constants;
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

		// Writing in the format: LSB -- MSB.
		if (i < no_state_var){
			if ((xi & 0x01) == 1){
				transition *= x[i];
			}
			else{
				transition *= !x[i];
			}

			if ((xi_ & 0x01) == 1){
				transition *= x_[i];
			}
			else{
				transition *= !x_[i];
			}
		}


		if (i < no_input_var){
			if ((ui & 0x01) == 1){
				transition *= u[i];
			}
			else{
				transition *= !u[i];
			}
		}

		//
		xi  >>= 1;
		ui  >>= 1;
		xi_ >>= 1;
	}

	return transition;
}


/* Important: The target set is given as a function of x. */
BDD getTargetSet(Cudd *mgr, int no_states, std::vector<int> target_set){

	//
	BDD minterm;
	BDD one  = mgr->bddOne();
	BDD zero = mgr->bddZero();
	BDD cofactor;
	BDD temp;

	unsigned int i;
	int j;
	int node;

	// Get the number of variables
	int no_state_vars = getNoBits(no_states - 1);

	BDD x[no_state_vars];

	// Create the variables
	for (i = 0; i < (unsigned int)no_state_vars; i++){
		x[i]  = mgr->bddVar(i);
	}


	/* Create the target set */
	cofactor = zero;

	for (i = 0; i < target_set.size(); i++){

		minterm = one;
		node    = target_set[i];

		for (j = 0; j < no_state_vars; j++){
			if (node & 0x01){
				temp = x[j].Ite(minterm, zero);
			}
			else{
				// row
				temp = x[j].Ite(zero, minterm);
			}
			minterm = temp;
			node >>= 1;
		}

		// Create the constant node.
		temp     = minterm.Ite(one, cofactor);
		cofactor = temp;
	}

	return cofactor;
} /* createTargetSet */


unsigned int getNoBits(unsigned int number){

	unsigned int no_bits = 0;

	for(;;){
		number >>= 1;
		no_bits++;
		if (number == 0)
			break;
	}
	return no_bits;
}

bool writeSysInfo_TBFile(int no_states, int no_inputs){
	std::ofstream myfile;

	myfile.open ("tb_sys_info.txt");
	myfile << no_states << std::endl;
	myfile << no_inputs;
	myfile.close();
	return true;
}

#ifdef ENABLE_TIME_PROFILING
long long get_usec(void){
	long long r;
	struct timeval t;


	gettimeofday(&t, NULL);
	r = t.tv_sec * 1000000 + t.tv_usec;
	return r;
}
#endif
























