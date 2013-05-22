/***********************************************************************

	PESSOA Version 1.0  
	------------------

Cyber-Physical Systems Laboratory
http://www.cyphylab.ee.ucla.edu
Author(s):	Manuel Mazo Jr. - mmazo@ee.ucla.edu
		Anna Davitian	- davitian@ee.ucla.edu

Dependencies:	CUDD library, MEX library
University of California, Los Angeles.
September 2009. 

************************************************************************/
 
#include "mex.h"
#include "pessoa.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "cudd.h"
#include "dddmp.h"


using namespace std;

static DdManager* ddman;     // global ddmanager pointer variable
static DdNode** ddnodearray; // global array of pointers to nodes 

const int BDD_W =0;

int numroots_design=1;    // number of roots/BDDs/functions

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), calls "buildBDDW", and plots the set W.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

	s_vector params_symb;
	int ret, ok;
	int type, verbose;
	double *wmin, *wmax;
	FILE *fp;
   	mxArray *psv, *psv_w;
	long nbatch,r,s;
	double totloops;
	double total_time;
	long memuse1, memuse2, mem_res;
	int plot_ok;

	char *bddW;
	char suffix[5];
	char *bddW_name;
	
	//BDD Variables		
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>4) {
		mexErrMsgTxt("Four inputs required: filename, min and max for W. Optional fourth input: verbose flag.");
		return;
	} 
	if(nrhs==4)
		verbose=(int)mxGetScalar(prhs[3]);
	else
		verbose=0;

	wmin=(double *)mxGetPr(prhs[1]);
	wmax=(double *)mxGetPr(prhs[2]);

	//Copy data from Matlab workspace 
	psv=mexGetVariable("caller","params_symb");

	//Copy data to variables
	nbatch=(long)mxGetScalar(mexGetVariable("caller","nbatch"));
	totloops=mxGetScalar(mexGetVariable("caller","totloops"));
	params_symb.n=(int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m=(int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume=(double *)mxGetPr(mxGetField(psv,0,"nume"));   
	params_symb.totbits=(int)mxGetScalar(mxGetField(psv,0,"totbits"));
	params_symb.nbitsloop=(int)mxGetScalar(mxGetField(psv,0,"nbitsloop"));
	params_symb.nbits=(double *)mxGetPr(mxGetField(psv,0,"nbits"));
	params_symb.deter=(int)mxGetScalar(mxGetField(psv,0,"deter"));
	params_symb.nbitsx=(int)mxGetScalar(mxGetField(psv,0,"nbitsx"));

	bddW=mxArrayToString(prhs[0]);

	strcpy(suffix,".bdd");

	bddW_name=(char*)mxMalloc(strlen(bddW)+5);
	strcpy(bddW_name, bddW);
	strcat(bddW_name,suffix);
	mxFree(bddW);

	//BDD INITIALIZATIONS
	numVars=params_symb.totbits;  // num of vars; if unknown set to 0
	numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
	numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
	cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
	//maxCacheSize=10485760*2;   // default for CUDD package

	ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, 0); //maxCacheSize);

	// BDD manager size
	memuse1=Cudd_ReadMemoryInUse(ddman);

 	ddnodearray = (DdNode**)malloc(numroots_design * sizeof(DdNode*));
 	for (r=0; r<numroots_design; r++) 
		ddnodearray[r] = Cudd_ReadLogicZero(ddman);

	// Plot first some info about the simulations:
	if(verbose>0)
	{
		mexEvalString("t_start = tic;");
		if(verbose==3)
			mexEvalString("profile on");
	}

	//buildBDDW returns number of states
	ret=buildBDDW(ddman,ddnodearray,BDD_W,wmin,wmax,&params_symb, nbatch);

	memuse2=Cudd_ReadMemoryInUse(ddman);

	// If target set is empty, give a warning message, beep, and stop
	if(ret==0){ 
		mexEvalString("beep");
		mexPrintf("\n----------------------- Pessoa: Target Set Aborted ---------------------  \n");
		for (s=0; s<numroots_design; s++) 
			Cudd_RecursiveDeref(ddman, ddnodearray[s]);
		Cudd_Quit(ddman);
		mexEvalString("warndlg('Target set is empty!','Pessoa Warning')");
		mexErrMsgTxt("\nERROR: Target set is empty! Process Stopped.\n");
	}

	if(verbose>0){
		// Print some extra-info
		mexPrintf("\nTarget set size: ");
		mexPrintf("%d", ret);
		mexPrintf(" states. \n");
		mexPrintf("Memory in use: %ld",memuse1);
		mexPrintf(" bytes (BDD manager); \n");
		memuse2=Cudd_ReadMemoryInUse(ddman);
		mem_res = memuse2-memuse1;
		if(mem_res<0)
			mem_res = 0;
		mexPrintf("               %ld",mem_res);
		mexPrintf(" bytes (symbolic set). \n");
	}

	//Dump BDD to file for viewing in DOT format
	if(verbose==3){
		fp=fopen("dumpfile_set2bdd.dot", "w");
		Cudd_DumpDot(ddman, numroots_design, ddnodearray,NULL,NULL,fp); 
		fclose(fp);
	}

	// Plotting
	if(verbose==2)
	{
		mexPrintf("\nPlotting target set...\n");
		plot_ok = PlotSet(ddman,ddnodearray,BDD_W, &params_symb);
		if(plot_ok == -1)
			mexPrintf("Plotting target set stopped... \n");
	}

	// Store the W set
	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[0], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddW_name, NULL);
	if(ok)
		mexPrintf("Target set successfully saved to '.bdd' file.\n");
	else
		mexPrintf("Target set failed to save.\n");
    

    for (s=0; s<numroots_design; s++) 
      Cudd_RecursiveDeref(ddman, ddnodearray[s]);

	Cudd_Quit(ddman);

	if(verbose>0)
	{
		mexEvalString("t_end = toc(t_start);");
		total_time = mxGetScalar(mexGetVariable("caller","t_end"));
		if(total_time < 0.001)
			mexPrintf("Elapsed time: 0.001 seconds. \n");
		else		
		{	
			mexPrintf("Elapsed time: %.3f", total_time);
			mexPrintf(" seconds. \n", total_time);
		}
	}

	if(verbose==3){	
		mexEvalString("profile off");
		mexEvalString("profsave(profile('info'),'profile_results')");}

	mexPrintf("\n--------------------- Pessoa: Target Set Terminated -------------------- \n");
}
