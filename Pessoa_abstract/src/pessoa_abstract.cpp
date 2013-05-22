/***********************************************************************

	PESSOA Version 1.4  
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

// BDD nodearray number for transitions
const int BDD_T = 0;

static DdManager* ddman;     // global ddmanager pointer variable
static DdNode** ddnodearray; // global array of pointers to nodes 

int numroots_abstract=1;    // number of roots/BDDs/functions

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), calls "buildBDDT", and plots the transitions.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

	s_vector params_symb;
	int ret,ok;
	FILE *fp;
	int verbose;
        mxArray *psv;
	double total_time;
	double ninputs,nstates;
	long nbatch,i,j,r,s;
	double totloops;
	long memuse1, memuse2, mem_res;
	int plot_ok;

	char *bddT;
	char suffix[5];
	char *bddT_name;

	//BDD Variables	
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package
	
	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>2) {
		mexErrMsgTxt("At least one string input required. Optional: a verbose level (default 0).");
		return;
	} 
	if(nrhs==2)
		verbose=(int)mxGetScalar(prhs[1]);
	else
		verbose=0;
	
	// Plot first some info about the simulations:
	if(verbose>0)
	{
		mexEvalString("t_start = tic;");
		if(verbose==3)
			mexEvalString("profile on");
	}

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

	//BDD INITIALIZATIONS		
	numVars=params_symb.totbits;  // num of vars; if unknown set to 0
	numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
	numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
	cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
	//maxCacheSize=10485760*2;   // default for CUDD package
		
	ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, 0); //maxCacheSize);

	// BDD manager size
	memuse1=Cudd_ReadMemoryInUse(ddman);

	//Each DdNode* needs to be properly init. Watch it for multi BDDs.
 	ddnodearray = (DdNode**)malloc(numroots_abstract * sizeof(DdNode*));
 	for (r=0; r<numroots_abstract; r++) 
		ddnodearray[r] = Cudd_ReadLogicZero(ddman);

       nstates=1;
       for (i=0;i<params_symb.n;i++)
	    	nstates*=(params_symb.nume[i]+1);
	
       ninputs=1;
       for (i=params_symb.n;i<params_symb.n+params_symb.m;i++)
	    	ninputs*=(params_symb.nume[i]+1);
	
	if(verbose>0){
		// Print some extra-info
		mexPrintf("\nSymbolic model size: ");
		mexPrintf("%.0f", nstates);
		mexPrintf(" states; \n");
		mexPrintf("                     %.0f", ninputs);
		mexPrintf(" inputs. \n");
	}

	// Build model abstraction
	ret=buildBDDT(ddman,ddnodearray,BDD_T,&params_symb, nbatch, totloops);

	if(ret==-1){
		mexEvalString("beep");
		mexPrintf("\n------------------------- Pessoa: Abstraction Aborted ----------------  \n");
		for (s=0; s<numroots_abstract; s++) 
			Cudd_RecursiveDeref(ddman, ddnodearray[s]);
		Cudd_Quit(ddman);
		mexErrMsgTxt("Pessoa: Abstraction Aborted");
	}

	if(verbose>0){
		mexPrintf("Memory in use: %ld",memuse1);
		mexPrintf(" bytes (BDD manager); \n");
		memuse2=Cudd_ReadMemoryInUse(ddman);
		mem_res = memuse2-memuse1;
		if(mem_res<0)
			mem_res = 0;
		mexPrintf("               %ld",mem_res);
		mexPrintf(" bytes (symbolic model). \n");
	}

	// Plotting
	if(verbose==2)
	{
		mexPrintf("\nPlotting Transitions...\n");
		plot_ok = PlotFSM(ddman,ddnodearray,BDD_T, &params_symb);
		if(plot_ok == -1)
			mexPrintf("Plotting Transitions stopped... \n");
	}

	if(verbose==3)
	{
		//Dump BDD to file for viewing in DOT format
		fp=fopen("dumpfile_abstr.dot", "w");
		Cudd_DumpDot(ddman, numroots_abstract, ddnodearray,NULL,NULL,fp); 
		fclose(fp);
	}


	bddT=mxArrayToString(prhs[0]);

	strcpy(suffix,".bdd");

	bddT_name=(char*)mxMalloc(strlen(bddT)+5);
	strcpy(bddT_name, bddT);
	strcat(bddT_name,suffix);
	mxFree(bddT);

	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[BDD_T], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddT_name, NULL);
	if(ok)
		mexPrintf("Symbolic model successfully saved to '.bdd' file.\n");
	else
		mexPrintf("Symbolic model failed to be saved.\n");
			

      for (s=0; s<numroots_abstract; s++) 
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
		
	mexPrintf("\n---------------------- Pessoa: Abstraction Terminated ------------------ \n");
}	
