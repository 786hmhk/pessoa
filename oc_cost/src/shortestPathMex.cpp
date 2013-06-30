/*
 * shortestPathMex.cpp
 *
 *  Created on: Jun 15, 2013
 *      Author: tanasaki
 */

#include "shortestPathMex.hh"

shortestPathMex::shortestPathMex() {
	// TODO Auto-generated constructor stub

}

shortestPathMex::~shortestPathMex() {
	// TODO Auto-generated destructor stub
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){


	unsigned int nstates, ninputs;
	int i;

	s_vector params_symb;
	mxArray *psv;

	char *bddsys       = "System_xux";
	char *bddsys_costs = "System_Cost_x";
	char *bddsys_tset  = "System_TargetSet_W";
	char *bddsys_name;
	char *bddsys_costs_name;
	char *bddsys_tset_name;
	char bdd_suffix[5];
	char add_suffix[5];



	mexPrintf("\n---------------- OC_SP: TEST BEGIN --------------- \n\n");

	/* Get the characteristics of the system. */

	// Copy data from MATLAB workspace.
	psv = mexGetVariable("caller","params_symb");

	params_symb.n    = (int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m    = (int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume = (double *)mxGetPr(mxGetField(psv,0,"nume"));


	// Get the number of states and inputs.
    nstates = 1;
	for (i = 0; i < params_symb.n; i++)
		nstates *= ((unsigned int)params_symb.nume[i]+1);

	ninputs = 1;
	for (i = params_symb.n; i < params_symb.n + params_symb.m; i++)
		ninputs *= ((unsigned int)params_symb.nume[i]+1);

	// Print some extra-info
	mexPrintf("\nSymbolic model size: ");
	mexPrintf("%u", nstates);
	mexPrintf(" states; \n");
	mexPrintf("                     %u", ninputs);
	mexPrintf(" inputs. \n");


	/* Check if all necessary files exist! */

	// assign the .bdd/.add file suffix
	strcpy(bdd_suffix,".bdd");
	strcpy(add_suffix,".add");

	/* Check if the System's BDD exists. */
	bddsys_name=(char*)mxMalloc(strlen(bddsys)+5);
	strcpy(bddsys_name, bddsys);
	strcat(bddsys_name, bdd_suffix);
//	mxFree(bddsys);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_name);
	FILE * smFile;
	smFile = fopen(bddsys_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the System's-States-Cost ADD file exists. */
	bddsys_costs_name=(char*)mxMalloc(strlen(bddsys_costs)+5);
	strcpy(bddsys_costs_name, bddsys_costs);
	strcat(bddsys_costs_name, add_suffix);
//	mxFree(bddsys);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_costs_name);
	smFile = fopen(bddsys_costs_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the Target Set BDD file exists. */
	bddsys_tset_name=(char*)mxMalloc(strlen(bddsys_tset)+5);
	strcpy(bddsys_tset_name, bddsys_tset);
	strcat(bddsys_tset_name, bdd_suffix);
//	mxFree(bddsys);

	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_tset_name);
	smFile = fopen(bddsys_tset_name,"r");
	FILE_EXISTS(smFile)


	/* Files exists... now load them. */

	// Initialize the Manager.
	Cudd mgr(0, 0);
	// Set background
	mgr.SetBackground(mgr.plusInfinity());
	//

	BDD S  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_name, NULL));

	// Create .dot file
	std::vector<BDD> nodes_bdd;
	nodes_bdd.push_back(S);
	FILE *outfile;
	outfile = fopen("System_xux.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();

	ADD SC = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_costs_name, NULL));

	// Create .dot file
	std::vector<ADD> nodes_add;
	nodes_add.push_back(SC);
	outfile = fopen("System_Cost_x.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	BDD W  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_tset_name, NULL));

	// Create .dot file
	nodes_bdd.push_back(W);
	outfile = fopen("System_TargetSet_W.dot", "w");
	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
	fclose(outfile);
	nodes_bdd.clear();


	/* Create the Shortest Path Object */
	mexPrintf("Creating SP Object ...\n");
	ShortestPath sp(&mgr, &S, nstates, ninputs); // optimized
//	ShortestPath sp(&mgr);

//	/* Create the Cost Adjacency Matrix */
	mexPrintf("Creating Cost Adjacency Matrix ...\n");
	ADD AG;
	AG = sp.createCostAdjacencyMatrix(&S, &SC, nstates, ninputs);

	// Create .dot file
	nodes_add.push_back(AG);
	outfile = fopen("System_CostAdjMatrix.dot", "w");
	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();

	/* Find All-pair Shortest Path */
	mexPrintf("Finding All-pair Shortest Path (FW) ...\n");
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
	mexPrintf("Finding All-pair Shortest Path to a Set W ...\n");
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

	mexPrintf("\n\n---------------- OC_SP: TEST END   --------------- \n");
}
