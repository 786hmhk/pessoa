/**
 * @file
 * @author Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>
 * @version 0.92
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
 * Contains...
 *
 * 		No details yet.
 */

#include "pessoa_dspController.hh"


//! Mex Function.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	printf("pessoa_dspController.cpp\n");

	unsigned int nstates, ninputs;

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

	/* Check if the System's BDD exists. */
	bddcntrl_ref_name = (char*)mxMalloc(strlen(sys_name)+5+11);
	strcpy(bddcntrl_ref_name, sys_name);
	strcat(bddcntrl_ref_name, "Controller");
	strcat(bddcntrl_ref_name, bdd_suffix);
	//check for file existence
//	mexPrintf("Checking %s ...\n", bddcntrl_ref_name);
	FILE * smFile;
	smFile = fopen(bddcntrl_ref_name,"r");
	FILE_EXISTS(smFile)

	/* Check if the System's-States-Cost ADD file exists. */
	bddsys_APSP_PA_name = (char*)mxMalloc(strlen(sys_name)+5+9);
	strcpy(bddsys_APSP_PA_name, sys_name);
	strcat(bddsys_APSP_PA_name, "_APSP_PA");
	strcat(bddsys_APSP_PA_name, add_suffix);
	//check for file existence
//	mexPrintf("Checking %s ...\n", bddsys_APSP_PA_name);
	smFile = fopen(bddsys_APSP_PA_name,"r");
	FILE_EXISTS(smFile)


	/* Check if the Target Set BDD file exists. */
	bddsys_APSP_PA_W_name = (char*)mxMalloc(strlen(sys_name)+5+11);
	strcpy(bddsys_APSP_PA_W_name, sys_name);
	strcat(bddsys_APSP_PA_W_name, "_APSP_PA_W");
	strcat(bddsys_APSP_PA_W_name, add_suffix);
	//check for file existence
//	mexPrintf("Checking %s ...\n", bddsys_APSP_PA_W_name);
	smFile = fopen(bddsys_APSP_PA_W_name,"r");
	FILE_EXISTS(smFile)



	/* Get the characteristics of the system. */

	// Copy data from MATLAB workspace.
	psv = mexGetVariable("caller","params_symb");

	params_symb.n    = (int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m    = (int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume = (double *)mxGetPr(mxGetField(psv,0,"nume"));


	//TODO: FIX THIS! This is not correct! Given wrong number of states/inputs.
//	// Get the number of states and inputs.
//    nstates = 1;
//	for (i = 0; i < params_symb.n; i++)
//		nstates *= ((unsigned int)params_symb.nume[i]+1);
//
//	ninputs = 1;
//	for (i = params_symb.n; i < params_symb.n + params_symb.m; i++)
//		ninputs *= ((unsigned int)params_symb.nume[i]+1);

	// Temp fix.
	double *nxbits = (double *)mxGetPr(mxGetField(psv,0,"nxbits"));
	double *nubits = (double *)mxGetPr(mxGetField(psv,0,"nubits"));

	// Get the number of states and inputs.
	unsigned int nxbits_total = 0;
	for (int i = 0; i < params_symb.n; i++)
		nxbits_total += (unsigned int)nxbits[i];


	unsigned int nubits_total = 0;
	for (int i = 0; i < params_symb.m; i++)
		nubits_total += (unsigned int)nubits[i];

	nstates = 2<<(nxbits_total - 1);
	ninputs = 2<<(nubits_total - 1);

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
	BDD S    = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, sys_name, NULL));
	BDD C    = BDD(mgr, Dddmp_cuddBddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddcntrl_ref_name, NULL));
	ADD PA   = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_APSP_PA_name, NULL));
	ADD PA_W = ADD(mgr, Dddmp_cuddAddLoad(mgr.getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_APSP_PA_W_name, NULL));

	/* Create the Shortest Path Object */
	ShortestPath sp(&mgr, &C, nstates, ninputs, false); // optimized

	/* Create the Controller. */
	BDD controller = sp.createControllerBDD(&S, &C, &PA, &PA_W);

	// Save as .bdd.
	if (sp.Dddmp_cuddStore(&controller, bddcntrl_ref_name)){
		mexPrintf("Symbolic controller successfully saved to '%s' file.\n", bddcntrl_ref_name);
	}


	/* De-allocate Memory */
	mxFree(bddcntrl_ref_name);
	mxFree(bddsys_APSP_PA_name);
	mxFree(bddsys_APSP_PA_W_name);


	mexPrintf("\n------------- Pessoa: Designing the DSP Controller Terminated ----------- \n");
}
