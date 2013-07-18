/*
 * pessoacontrollerc.cpp
 *
 *  Created on: Jul 14, 2013
 *      Author: tanasaki
 */

#include "pessoa_controller.hh"


//! Mex Function.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	unsigned int nstates, ninputs;
	int i;

	s_vector params_symb;
	mxArray *psv;

	char *sys_name;
	char *bddcntrl_ref_name;
	char *bddsys_APSP_PA_name;
	char *bddsys_APSP_PA_W_name;
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


	/* Check if all necessary files exist! */

	// TODO: change all these strcpy, strcat... make it nicer plz :). thx.
	// assign the .bdd/.add file suffix
	strcpy(bdd_suffix,".bdd");
	strcpy(add_suffix,".add");

//	/* Check if the System's BDD exists. */
//	bddcntrl_ref_name = (char*)mxMalloc(strlen(sys_name)+5+15);
//	strcpy(bddcntrl_ref_name, sys_name);
//	strcat(bddcntrl_ref_name, "Controller_dom");
//	strcat(bddcntrl_ref_name, bdd_suffix);
//	//check for file existence
//	mexPrintf("Checking %s ...\n", bddcntrl_ref_name);
//	FILE * smFile;
//	smFile = fopen(bddcntrl_ref_name,"r");
//	FILE_EXISTS(smFile)


	/* Check if the System's BDD exists. */
	bddcntrl_ref_name = (char*)mxMalloc(strlen(sys_name)+5);
	strcpy(bddcntrl_ref_name, sys_name);
	strcat(bddcntrl_ref_name, bdd_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddcntrl_ref_name);
	FILE * smFile;
	smFile = fopen(bddcntrl_ref_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the System's-States-Cost ADD file exists. */
	bddsys_APSP_PA_name = (char*)mxMalloc(strlen(sys_name)+5+9);
	strcpy(bddsys_APSP_PA_name, sys_name);
	strcat(bddsys_APSP_PA_name, "_APSP_PA");
	strcat(bddsys_APSP_PA_name, add_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_APSP_PA_name);
	smFile = fopen(bddsys_APSP_PA_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the Target Set BDD file exists. */
	bddsys_APSP_PA_W_name = (char*)mxMalloc(strlen(sys_name)+5+11);
	strcpy(bddsys_APSP_PA_W_name, sys_name);
	strcat(bddsys_APSP_PA_W_name, "_APSP_PA_W");
	strcat(bddsys_APSP_PA_W_name, add_suffix);
	//check for file existence
	mexPrintf("Checking %s ...\n", bddsys_APSP_PA_W_name);
	smFile = fopen(bddsys_APSP_PA_W_name,"r");
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

	if (verbose == 3){
		// Print some extra-info
		mexPrintf("\nSymbolic model size: ");
		mexPrintf("%u", nstates);
		mexPrintf(" states; \n");
		mexPrintf("                     %u", ninputs);
		mexPrintf(" inputs. \n");
	}

	/* Files exist... now load them. */

	// Initialize the Manager.
	Cudd mgr(0, 0);
	// Set background
	mgr.SetBackground(mgr.plusInfinity());

	// Loading .bdd and .add files.
	BDD CR   = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddcntrl_ref_name, NULL));
	ADD PA   = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_APSP_PA_name, NULL));
	ADD PA_W = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_APSP_PA_W_name, NULL));

	/* Create the Shortest Path Object */
	mexPrintf("Creating SP Object ...\n");
	ShortestPath sp(&mgr, &CR, nstates, ninputs); // optimized


	/* Create the Controller. */
	BDD controller = sp.createControllerBDD(&CR, &PA, &PA_W);


	// Save as .bdd.
	sp.Dddmp_cuddStore(&controller, bddcntrl_ref_name);


	/* De-allocate Memory */
	mxFree(bddcntrl_ref_name);
	mxFree(bddsys_APSP_PA_name);
	mxFree(bddsys_APSP_PA_W_name);
}
