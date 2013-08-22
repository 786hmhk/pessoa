/*
 * pessoandsp.cpp
 *
 *  Created on: Aug 22, 2013
 *      Author: tanasaki
 */

#include "pessoa_ndsp.hh"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){


	unsigned int nstates, ninputs;
	int i;

	s_vector params_symb;
	mxArray *psv;

	char *sys_name;
	char *bddsys_name;
	char *bddsys_costs_name;
	char *bddsys_tset_name;
	char bdd_suffix[5];
	char add_suffix[5];

	int verbose;



	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>2) {
		mexErrMsgTxt("Two inputs required: system filename. Optional second input: verbose flag.");
		return;
	}
	if(nrhs==2)
		verbose = (int)mxGetScalar(prhs[1]);
	else
		verbose = 0;

	// get the system's name.
	sys_name = mxArrayToString(prhs[0]);

	if (verbose > 0){
	}


	/* Check if all necessary files exist! */

	// TODO: change all these strcpy, strcat... make it nicer plz :). thx.
	// assign the .bdd/.add file suffix
	strcpy(bdd_suffix,".bdd");
	strcpy(add_suffix,".add");

//	/* Check if the System's BDD exists. */
//	bddsys_name=(char*)mxMalloc(strlen(sys_name)+5);
//	strcpy(bddsys_name, sys_name);
//	strcat(bddsys_name, bdd_suffix);
//	//check for file existence
//	mexPrintf("Checking %s ...\n", bddsys_name);
//	FILE * smFile;
//	smFile = fopen(bddsys_name,"r");
//	FILE_EXISTS(smFile)

	/* Check if the System's BDD exists. */
	bddsys_name=(char*)mxMalloc(strlen(sys_name)+5+10);
	strcpy(bddsys_name, sys_name);
	strcat(bddsys_name, "Controller");
	strcat(bddsys_name, bdd_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_name);
	FILE * smFile;
	smFile = fopen(bddsys_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the System's-States-Cost ADD file exists. */
	bddsys_costs_name=(char*)mxMalloc(strlen(sys_name)+5+6);
	strcpy(bddsys_costs_name, sys_name);
	strcat(bddsys_costs_name, "Costs");
	strcat(bddsys_costs_name, add_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_costs_name);
	smFile = fopen(bddsys_costs_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the Target Set BDD file exists. */
	bddsys_tset_name=(char*)mxMalloc(strlen(sys_name)+5+10);
	strcpy(bddsys_tset_name, sys_name);
	strcat(bddsys_tset_name, "TargetSet");
	strcat(bddsys_tset_name, bdd_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_tset_name);
	smFile = fopen(bddsys_tset_name,"r");
	FILE_EXISTS(smFile)



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


	/* Files exist... now load them. */

	// Initialize the Manager.
	Cudd mgr(0, 0);
	// Set background
	mgr.SetBackground(mgr.plusInfinity());

	// Loading .bdd and .add files.
	BDD S  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_name, NULL));
	ADD SC = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_costs_name, NULL));
	BDD W  = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_tset_name, NULL));

	/* Create the Shortest Path Object */
	ShortestPath sp(&mgr, &S, nstates, ninputs, false); // optimized
//	ShortestPath sp(&mgr);


	/* Find the All-pair Shortest Path to  a Set. */
	mexPrintf("Finding All-pair Shortest Path to a Set W ...\n");
	ADD APSP_W;
	BDD PA_W;
	sp.APtoSetSP(&S, &SC, &W, &APSP_W, &PA_W, nstates, ninputs);

	mexPrintf("Calculating ND-SP done!\n");


	/* Dump results into files. */
	if (verbose == 3){
		mexPrintf("Dumping Controller into (.bdd) file...\n");
	}


	// Create the file name.
	char *name = (char*)mxMalloc(strlen(sys_name)+5+9);
	strcpy(name, sys_name);
	strcat(name, "Controller");
	strcat(name, ".bdd");

	DdNode *fn = PA_W.getNode();
	/* Dump the PA ADD to an .add file. */
	if (!Dddmp_cuddBddStore(mgr.getManager(), NULL, fn, NULL, NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS, name, NULL))
		mexErrMsgTxt("Error dumping the PA_W (Controller) BDD into an .bdd file!");

	mxFree(name);







	/* Create .dot files. */
//	FILE *outfile;
//	std::vector<BDD> nodes_bdd;
//	std::vector<ADD> nodes_add;
//
//
//	// Create .dot file
//	nodes_bdd.push_back(S);
//	outfile = fopen("DCMotor.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();
//
//	// Create .dot file
//	nodes_bdd.push_back(W);
//	outfile = fopen("DCMotorTargetSet.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();

//
//	// Create .dot file
//	nodes_add.push_back(APSP);
//	outfile = fopen("System_APSP.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//


//	// Create .dot file
//	nodes_add.push_back(APSP_W);
//	outfile = fopen("DCMotorController_APSP_W.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
//	// Create .dot file
//	nodes_bdd.push_back(PA_W);
//	outfile = fopen("DCMotorController.dot", "w");
//	mgr.DumpDot(nodes_bdd, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_bdd.clear();


	/* De-allocate Memory */
	mxFree(bddsys_name);
	mxFree(bddsys_costs_name);
	mxFree(bddsys_tset_name);
	mxFree(sys_name);
}
