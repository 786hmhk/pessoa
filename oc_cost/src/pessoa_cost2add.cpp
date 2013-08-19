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


	// Get the number of states and inputs.
    nstates = 1;
	for (int i = 0; i < params_symb.n; i++)
		nstates *= ((unsigned int)params_symb.nume[i]+1);

	ninputs = 1;
	for (int i = params_symb.n; i < params_symb.n + params_symb.m; i++)
		ninputs *= ((unsigned int)params_symb.nume[i]+1);

	no_figure = 1;

	// Get & create the filename of the ADD holding the state costs.
	SysBDD_name = mxArrayToString(prhs[0]);
	SysStateCostADD_name = (char*)mxMalloc(strlen(SysBDD_name)+6);
	strcpy(SysStateCostADD_name, SysBDD_name);
	strcat(SysStateCostADD_name, "Costs");

	// CUDD Manager.
	this->mgr = mgr;
	// Set background
	mgr->SetBackground(mgr->plusInfinity());
	// Initiate the Cost ADD
	sys_state_cost = mgr->background();

	total_add_states = 0;
}

pessoa_cost2add::~pessoa_cost2add() {
	mxFree(SysStateCostADD_name);
	mxFree(SysBDD_name);
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
				tempv[i]=wmin[i];
				i++;
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

	total_add_states = (unsigned int) total_states;
//	mexPrintf("Total ADD states: %d\n", total_add_states);

#ifdef WAITBAR_ENABLE
	mexEvalString("delete(h);");
#endif

	// TODO: FREE MEMORY.
	mxDestroyArray(pstate_array);
	mxDestroyArray(pset_state_array);
	mxDestroyArray(pjnbatch);
	mxDestroyArray(pstate_cost_array);

	mxFree(tempv);
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

//! Plots the System's State Cost. Currently supporting only 3D.
void pessoa_cost2add::plotSysStateCost(){

	mexPrintf("Plotting ADD... \n");

	double *tempv, cost = 0;

	long i;


	if (params_symb.n != 2) {
		mexEvalString(
				"warndlg(Can't plot systems other than 2-D!','Pessoa Error')");
		return;
	}

	//
	mxArray *states_array_mx = mxCreateNumericMatrix(total_add_states, params_symb.n, mxDOUBLE_CLASS, mxREAL);
	mxArray *cost_array_mx   = mxCreateNumericMatrix(total_add_states, 1, mxDOUBLE_CLASS, mxREAL);
	double *states_array    = (double *) mxGetData(states_array_mx);
	double *cost_array      = (double *) mxGetData(cost_array_mx);

	tempv = (double *) malloc(params_symb.n * sizeof(double));
	for (i = 0; i < params_symb.n; i++)
		tempv[i] = 0;



#ifdef WAITBAR_ENABLE // TODO: Have not tested it! Also take care of Memory leakage...
	int totbits = params_symb.totbits;
	double *ps_steps, *stotstates;
	long j;
	double totstates = 1;

	for (j = 0; j < n; j++)
		totstates *= params_symb.nume[j] + 1;

	mxArray *pps_steps, *pstotstates;
	pps_steps     = mxCreateDoubleScalar(0);
	pstotstates   = mxCreateDoubleScalar(totstates);
	ps_steps      = (double *) mxGetData(pps_steps);
	stotstates    = (double *) mxGetData(pstotstates);
	ps_steps[0]   = 0;
	stotstates[0] = totstates;

	mexPutVariable("caller", "pstotstates", pstotstates);
	mexEvalString(
			"h=waitbar(0,'Plotting set in progress','Name','Pessoa','CreateCancelBtn','global stopbuild;');");
#endif

	unsigned int ii = 0;
	while (1) {

		if (isInADD(tempv, &cost)) {

			// Get the state.
			for (int m = 0; m < params_symb.n; m++){
				states_array[ii + m*total_add_states] = tempv[m];
			}
			// Get the State Cost.
			cost_array[ii] = cost;
			ii++;
		}
#ifdef WAITBAR_ENABLE
		ps_steps[0]++;
#endif

		i = 0;
		while (i < params_symb.n && (++tempv[i]) > params_symb.nume[i]) {
			tempv[i++] = 0;
		}
		if (i == params_symb.n) {
			break;
		}

#ifdef WAITBAR_ENABLE
		// Used for waitbar
		mexPutVariable("caller", "ps_steps", pps_steps);
		mexEvalString("waitbar(ps_steps/pstotstates,h)");
		if (mexGetVariable("base", "stopbuild") != NULL) // This is a dirty Hack, but I don't know any other way.
		{
			mexEvalString("delete(h);");
			return;
		}
#endif

	}

	for (unsigned int kk = ii; kk < total_add_states; kk++){
		cost_array[kk] = mxGetInf();
	}

	// Plot the Results.
	mexPutVariable("caller", "states_array", states_array_mx);
	mexPutVariable("caller", "cost_array",   cost_array_mx);
	mexEvalString(	"for ii = 1:length(states_array(:,1)) states_array(ii,:) = params_symb.eta * (states_array(ii,:)' + params_symb.min(params_symb.xoind)); end");
//	mexEvalString("plot3(states_array(:,1), states_array(:,2), cost_array, 'r+'); grid on;");

//	xlabel('length'); ylabel('width'); zlabel('height');
//	mexEvalString("states_costs_ = [states_array cost_array]");

	switch (no_figure){

	case 1:
		mexEvalString(
				"h = figure('name','State Costs of the System'); scatter3(states_array(:,1), states_array(:,2), cost_array, 5, cost_array); colormap(jet); xlabel('State Variable x(1)'); ylabel('State Variable x(2)'); zlabel('State Cost'); "
						"saveas(h,'costs_orig','fig')");
		break;
	case 2:
		mexEvalString(
				"h = figure('name','State Costs of the System'); scatter3(states_array(:,1), states_array(:,2), cost_array, 5, cost_array); colormap(jet); xlabel('State Variable x(1)'); ylabel('State Variable x(2)'); zlabel('State Cost'); "
						"saveas(h,'costs_vs','fig')");
		break;
	case 3:
		mexEvalString(
				"h = figure('name','State Costs of the System'); scatter3(states_array(:,1), states_array(:,2), cost_array, 5, cost_array); colormap(jet); xlabel('State Variable x(1)'); ylabel('State Variable x(2)'); zlabel('State Cost'); "
						"saveas(h,'costs_vt','fig')");
		break;
	}

	no_figure++;




#ifdef WAITBAR_ENABLE
	mexEvalString("delete(h);");
#endif

	// TODO: FREE MEMORY.
	mxDestroyArray(states_array_mx);
	mxDestroyArray(cost_array_mx);
	free(tempv);
}

//! Checks whether a current state has some cost assigned to it.
bool pessoa_cost2add::isInADD(double *array, double *cost){
	long i;
	double val_bool = 0;
	int k;

	ADD temp;
	ADD minterm, restrct;

	k = 0;
	minterm = mgr->addOne();

	for (i = 0; i < params_symb.n; i++) {
		int ind = 1 << ((int) params_symb.nbits[i] - 1);
		while (ind > 0) {
			val_bool = (int) array[i] & ind;
			//
			if (val_bool) {
				temp = minterm * mgr->addVar(k);

			} else {
				temp = minterm * (~mgr->addVar(k));
			}

			minterm = temp;
			k++;
			ind >>= 1;
		}
	}


	restrct = sys_state_cost.Restrict(minterm);

	if (restrct == mgr->background()){
		return false;
	}
	else{
		*cost = (double)restrct.getNode()->type.value;
		return true;
	}
	return false;
}

//! Returns the System's Cost ADD.
ADD pessoa_cost2add::getSysStateCost(){
	return sys_state_cost;
}

//! Sets the System's Cost ADD.
void pessoa_cost2add::setSysStateCost(ADD sysCosts){
	sys_state_cost = sysCosts;
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

//! Experimental. Filters out (assigns infinity) to states that are invalid.
ADD pessoa_cost2add::filterCosts(int mode){

	char *SysBDD_filename = (char*)mxMalloc(strlen(SysBDD_name)+5);
	strcpy(SysBDD_filename, SysBDD_name);
	strcat(SysBDD_filename, ".bdd");

	// check for file existence
	mexPrintf("Checking %s ...\n", SysBDD_filename);
	FILE * smFile;
	smFile = fopen(SysBDD_filename,"r");
	FILE_EXISTS(smFile)

	// Loading .bdd and .add files.
	BDD S = BDD(*mgr, Dddmp_cuddBddLoad(mgr->getManager(), DDDMP_VAR_MATCHIDS, NULL, NULL, NULL, DDDMP_MODE_DEFAULT, SysBDD_filename, NULL));

	/* Create the Shortest Path Object */
	ShortestPath sp(mgr, &S, nstates, ninputs);

	return sp.filterCosts(&S, &sys_state_cost, mode);
}


//! Mex Function.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	int verbose;

	/* Initialize the Manager. */
	Cudd mgr(0, 0);

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

	if (verbose == 3){
		// Create .dot file
		cost2add.dumpSysStateCostDot();

		// Plot ADD.
		cost2add.plotSysStateCost();
	}

	if (verbose == 4){

		ADD backup = cost2add.getSysStateCost();

		// Plot Initial ADD.
		cost2add.plotSysStateCost();

		// KEEP_VALID_STATES
		cost2add.setSysStateCost(cost2add.filterCosts(KEEP_VALID_STATES));
		// Plot ADD.
		cost2add.plotSysStateCost();

		// restore ADD.
		cost2add.setSysStateCost(backup);

		// KEEP_VALID_TRANSITIONS
		cost2add.setSysStateCost(cost2add.filterCosts(KEEP_VALID_TRANSITIONS));
		// Plot ADD.
		cost2add.plotSysStateCost();

		// restore ADD.
		cost2add.setSysStateCost(backup);

		// Filter/fix the figures.
		mexEvalString("analyzeCostFig");
	}

	/* Dump the ADD to an .add file. */
	if (!cost2add.dumpSysStateCost())
		mexErrMsgTxt("Error dumping the Systems' Cost ADD into an .add file!");


	mexPrintf("\n------------------ Pessoa: Creating Cost ADD Terminated ------------------- \n");
}
