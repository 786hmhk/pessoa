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

#include "pessoa_dsp.hh"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){


	printf("pessoa_dsp.cpp\n");

	unsigned int nstates, ninputs;

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


	//TODO: FIX THIS! This is not correct! Given wrong number of states/inputs.
//	// Get the number of states and inputs.
//    nstates = 1;
//	for (int i = 0; i < params_symb.n; i++)
//		nstates *= ((unsigned int)params_symb.nume[i]+1);
//
//	ninputs = 1;
//	for (int i = params_symb.n; i < params_symb.n + params_symb.m; i++)
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

	// Print some extra-info
	mexPrintf("\nSymbolic model size: ");
	mexPrintf("%u", nstates);
	mexPrintf(" states; \n");
	mexPrintf("                     %u", ninputs);
	mexPrintf(" inputs. \n");


	/* Files exist... now load them. */

	// Initialize the Manager.
	unsigned int numVars = (2 * (nxbits_total - 1)) + (nubits_total - 1);
	Cudd *mgr = new Cudd(numVars, 0);

	// Set background
	mgr->SetBackground(mgr->plusInfinity());

	// Loading .bdd and .add files.
	BDD S  = BDD(*mgr, Dddmp_cuddBddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_name, NULL));
	ADD SC = ADD(*mgr, Dddmp_cuddAddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_costs_name, NULL));

	/* Create the Shortest Path Object */
	ShortestPath *sp = new ShortestPath(mgr, &S, nstates, ninputs, false); // optimized

	/* Create the Cost Adjacency Matrix */
	mexPrintf("Creating Cost Adjacency Matrix ...\n");
	ADD AG = sp->createCostAdjacencyMatrix(&S, &SC, nstates, ninputs);

	/* Find All-pair Shortest Path */
	mexPrintf("Finding All-pair Shortest Path (FW) ...\n");
	ADD APSP;
	ADD PA;
	sp->FloydWarshall(&AG, &APSP, &PA);

	// Create the file name.
	char *apsp_tmp_name = (char*)mxMalloc(strlen(sys_name)+5+10);
	strcpy(apsp_tmp_name, sys_name);
	strcat(apsp_tmp_name, "_APSP_tmp");
	strcat(apsp_tmp_name, ".add");
	sp->Dddmp_cuddStore(&APSP, apsp_tmp_name);

	char *apsp_pa_tmp_name = (char*)mxMalloc(strlen(sys_name)+5+13);
	strcpy(apsp_pa_tmp_name, sys_name);
	strcat(apsp_pa_tmp_name, "_APSP_PA_tmp");
	strcat(apsp_pa_tmp_name, ".add");
	sp->Dddmp_cuddStore(&PA, apsp_pa_tmp_name);

	/*
	 * Necessary evil/hack! I need to destroy and re-create the manager because of the
	 * weird behaviour of the cudd lib, when it comes for the carbage collector.
	 */

	delete sp;
	delete mgr;

	mgr = new Cudd(numVars, 0);
	BDD S1 = BDD(*mgr, Dddmp_cuddBddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_name, NULL));
	BDD W  = BDD(*mgr, Dddmp_cuddBddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, bddsys_tset_name, NULL));

	ADD APSP1 = ADD(*mgr, Dddmp_cuddAddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, apsp_tmp_name, NULL));
	ADD PA1   = ADD(*mgr, Dddmp_cuddAddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, apsp_pa_tmp_name, NULL));

	delete apsp_tmp_name;
	delete apsp_pa_tmp_name;

	sp = new ShortestPath(mgr, &S1, nstates, ninputs, false); // optimized

	/* Find the All-pair Shortest Path to  a Set. */
	mexPrintf("Finding All-pair Shortest Path to a Set W ...\n");
	ADD APSP_W;
	ADD PA_W;
	sp->APtoSetSP(&APSP1, &PA1, &W, &APSP_W, &PA_W);

	mexPrintf("Calculating D-SP done!\n");


	/* Dump results into files. */
	if (verbose == 3){
		mexPrintf("Dumping PA and PA_W into files...\n");
	}


	// Create the file name.
	char *name = (char*)mxMalloc(strlen(sys_name)+5+9);
	strcpy(name, sys_name);
	strcat(name, "_APSP_PA");
	strcat(name, ".add");

	DdNode *fn = PA1.getNode();
	/* Dump the PA ADD to an .add file. */
	if (!Dddmp_cuddAddStore(mgr->getManager(), NULL, fn, NULL, NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS, name, NULL))
		mexErrMsgTxt("Error dumping the APSP PA ADD into an .add file!");

	delete name;

	name = (char*)mxMalloc(strlen(sys_name)+5+11);
	strcpy(name, sys_name);
	strcat(name, "_APSP_PA_W");
	strcat(name, ".add");
	fn = PA_W.getNode();
	/* Dump the PA_W ADD to an .add file. */
	if (!Dddmp_cuddAddStore(mgr->getManager(), NULL, fn, NULL, NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS, name, NULL))
		mexErrMsgTxt("Error dumping the APSP PA_W ADD into an .add file!");

	delete name;





//	/* Create .dot files. */
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
//	nodes_add.push_back(AG);
//	outfile = fopen("System_CostAdjMatrix.dot", "w");
//	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
//	fclose(outfile);
//	nodes_add.clear();
////
////	// Create .dot file
////	nodes_add.push_back(APSP);
////	outfile = fopen("System_APSP.dot", "w");
////	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
////	fclose(outfile);
////	nodes_add.clear();
////
////	// Create .dot file
////	nodes_add.push_back(PA);
////	outfile = fopen("System_APSP_PA.dot", "w");
////	mgr.DumpDot(nodes_add, NULL, NULL, outfile);
////	fclose(outfile);
////	nodes_add.clear();
//
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


	/* De-allocate Memory */
	delete bddsys_name;
	delete bddsys_costs_name;
	delete bddsys_tset_name;
	delete sys_name;

	delete sp;
	delete mgr;
}

