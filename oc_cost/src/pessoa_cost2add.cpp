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
 * Contains...
 *
 * 		No details yet.
 */

#include "pessoa_cost2add.hh"

pessoa_cost2add::pessoa_cost2add(mxArray **plhs, const mxArray **prhs, Cudd *mgr, int verbose) {

	this->verbose = verbose;

	// Mex function inputs.
	wmin = (double *)mxGetPr(prhs[1]);
	wmax = (double *)mxGetPr(prhs[2]);

	// Copy data from MATLAB workspace.
	nbatch   = (long)mxGetScalar(mexGetVariable("caller","nbatch"));
	psv      = mexGetVariable("caller","params_symb");
#ifdef WAITBAR_ENABLE
	totloops = mxGetScalar(mexGetVariable("caller","totloops"));
#endif
	params_symb.n         = (int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m         = (int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume      = (double *)mxGetPr(mxGetField(psv,0,"nume"));
	params_symb.totbits   = (int)mxGetScalar(mxGetField(psv,0,"totbits"));
	params_symb.nbitsloop = (int)mxGetScalar(mxGetField(psv,0,"nbitsloop"));
	params_symb.nbits     = (double *)mxGetPr(mxGetField(psv,0,"nbits"));
	params_symb.deter     = (int)mxGetScalar(mxGetField(psv,0,"deter"));
	params_symb.nbitsx    = (int)mxGetScalar(mxGetField(psv,0,"nbitsx"));

	// Get & create the filename of the ADD holding the state costs.
	SysStateCostADD_name = mxArrayToString(prhs[0]);

	// CUDD Manager.
	this->mgr = mgr;
	// Set background
	mgr->SetBackground(mgr->plusInfinity());
	// Initiate the Cost ADD
	sys_state_cost = mgr->background();

}

pessoa_cost2add::~pessoa_cost2add() {
	mxFree(SysStateCostADD_name);
}

//! Creates the ADD holding the System's States Cost.
void pessoa_cost2add::createSysStatesCost(){

	int end = 0;
	long i, j = 0;
	//
	double *tempv, *state_array, *state_cost_array, *set_state_array;
#ifdef WAITBAR_ENABLE
	double *wfsteps;
#endif
//	double *jnbatch; // TODO: DELETE This... (?)
	//
	mxArray *pstate_array, *pset_state_array, *pjnbatch, *pstate_cost_array;
	//
	long nstates;
	long total_states=0;

	tempv             = (double *)mxMalloc(params_symb.n*sizeof(double));
	state_array       = (double *)mxMalloc(nbatch*params_symb.n*sizeof(double));

	memcpy(tempv, wmin, params_symb.n*sizeof(double));

	pstate_array      = mxCreateNumericMatrix(params_symb.n,nbatch,mxDOUBLE_CLASS,mxREAL);
	pjnbatch          = mxCreateDoubleScalar(nbatch);

//	jnbatch           = (double *) mxGetData(pjnbatch);
	state_array       = (double *) mxGetData(pstate_array);


#ifdef WAITBAR_ENABLE
	// Wait bar.
	if (mexEvalString("h = waitbar(0,'Assigning Costs to the System's States','Name','Pessoa','CreateCancelBtn','global stopbuild');"))
		mexPrintf("ERRORRRRRRRRRRRRRRRRRRRRRRRRR\n\n");
	mxArray *pwfsteps;
	pwfsteps    = mxCreateDoubleScalar(j);
	wfsteps     = (double *) mxGetData(pwfsteps);
	wfsteps[0]  = 0;
#endif



	//Initialize
	mexPutVariable("caller", "jnbatch", pjnbatch);
	while(!end)
	{
		for(j=0; j<nbatch; j++)
		{
			memcpy(&state_array[j*params_symb.n], tempv, params_symb.n*sizeof(double));
			i=0;
			while(i<params_symb.n && (++tempv[i])>wmax[i])
			{
				tempv[i++]=wmin[i];
			}
			if(i==params_symb.n)
			{
				end=1;
				j++;
				break;
			}

		}

		mexPutVariable("caller", "state_array", pstate_array);

		// TOFIX: In the case of an Error below, the BDD manager will not be freed, causing a memory leak!!
		if(mexEvalString("[set_state_array state_cost_array nstates] = pss_build_charfcost(params_symb,state_array,jnbatch);")){
			mexEvalString("warndlg('pss_build_charfcost.m failed or missing!','Pessoa Error')");
#ifdef WAITBAR_ENABLE
			mexEvalString("delete(h); clear global stopbuild;");
#endif
			mexErrMsgTxt("ERROR: pss_build_charfcost.m failed or missing!");
		}

		pset_state_array  = mexGetVariable("caller","set_state_array");
		pstate_cost_array = mexGetVariable("caller","state_cost_array");
		set_state_array   = mxGetPr(pset_state_array);
		state_cost_array  = mxGetPr(pstate_cost_array);

		nstates = mxGetScalar(mexGetVariable("caller","nstates"));


		addBatchToADD(set_state_array, state_cost_array, &params_symb, nstates, params_symb.n);

		total_states = total_states + nstates;

#ifdef WAITBAR_ENABLE
		// Used for waitbar
		wfsteps[0]++;
		mexPutVariable("caller", "wfsteps", pwfsteps);
		mexEvalString("waitbar(wfsteps/totloops,h)");
		if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
		{
			mexEvalString("delete(h); clear global stopbuild;");
		}
#endif
	}

#ifdef WAITBAR_ENABLE
	mexEvalString("delete(h);");
#endif
}


//! Adds the next missing batch to the ADD.
bool pessoa_cost2add::addBatchToADD(double *array, double *state_cost_array, s_vector* params_symb, long nbatch, int nnm){
	long i, j;
	double val_bool = 0;
	int k;

	ADD temp;
	ADD minterm;

	for (j = 0; j < nbatch; j++) {

		k = 0;
		minterm = mgr->addOne();

		for (i = 0; i < nnm; i++) {
			int ind = 1 << ((int) params_symb->nbits[i] - 1);
			while (ind > 0) {
				val_bool = (int) array[j * nnm + i] & ind;
				//
				if (val_bool) {
					temp = minterm * mgr->addVar(k);

				}
				else{
					temp = minterm * (~mgr->addVar(k));
				}

				minterm = temp;
				k++;
				ind >>= 1;

			}
		}

		// Create the constant node.
		temp = minterm.Ite(mgr->constant(state_cost_array[j]), sys_state_cost);
		sys_state_cost = temp;

	}

	return true;
}

//! Returns the Ssstem's Cost ADD.
ADD pessoa_cost2add::getSysStateCost(){
	return sys_state_cost;
}

//! Dumps the System's State Cost ADD to a .add file.
bool pessoa_cost2add::dumpSysStateCost(){
	bool ok;

	if (verbose == 3){
		mexPrintf("Creating .add file...\n");
	}

	// Create the file name.
	char *name = (char*)mxMalloc(strlen(SysStateCostADD_name)+5);
	strcpy(name, SysStateCostADD_name);
	strcat(name, ".add");

	DdNode *fn = sys_state_cost.getNode();
	ok = Dddmp_cuddAddStore(mgr->getManager(), NULL, fn, NULL, NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS, name, NULL);

	mxFree(name);

	if (ok) return true;
	else return false;
}

//! Dumps the System's State Cost ADD to a .dot file.
void pessoa_cost2add::dumpSysStateCostDot(){

	if (verbose == 3){
		mexPrintf("Creating .dot file...\n");
	}
	// Create .dot file
	std::vector<ADD> nodes_add;
	FILE *outfile;
	nodes_add.push_back(sys_state_cost);
	// Create the file name.
	char *name = (char*)mxMalloc(strlen(SysStateCostADD_name)+5);
	strcpy(name, SysStateCostADD_name);
	strcat(name, ".dot");
	outfile = fopen(name, "w");
	mgr->DumpDot(nodes_add, NULL, NULL, outfile);
	fclose(outfile);
	nodes_add.clear();
	mxFree(name);
}


//! Mex Function.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	int verbose;

	/* Initialize the Manager. */
	Cudd mgr(0, 0);


	mexPrintf("\n------------------ pessoa_cost2add: TEST BEGIN ----------------- \n\n");

	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>4) {
		mexErrMsgTxt("Four inputs required: filename, min and max for the System's state space. Optional fourth input: verbose flag.");
		return;
	}
	if(nrhs==4)
		verbose = (int)mxGetScalar(prhs[3]);
	else
		verbose = 0;


	/* Create the pessoa_cost2add object. */
	pessoa_cost2add cost2add(plhs, prhs, &mgr, verbose);
	/* Create the ADD. */
	cost2add.createSysStatesCost();
	/* Dump the ADD to a .add file. */
	if (!cost2add.dumpSysStateCost())
		mexErrMsgTxt("Error dumping the Systems' Cost ADD into an .add file!");

	// Create .dot file
	cost2add.dumpSysStateCostDot();


	mexPrintf("\n\n------------------ pessoa_cost2add: TEST END   ----------------- \n");
	mexPrintf("\n------------------ Pessoa: Creating Cost ADD Terminated ------------------- \n");
}