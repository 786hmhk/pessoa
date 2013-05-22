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

const int BDD_T1=0;
const int BDD_T2=1;
const int BDD_W=2;
const int BDD_TR=3;
const int BDD_WR=4;

int numroots_design=5;    // number of roots/BDDs/functions

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), computes combination controller, and plots it.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	s_vector params_symb;
	int ok;
	int type, type_flag,verbose;
	double *wmin, *wmax;
	FILE *fp;
   	mxArray *psv, *psv_w;
	long nbatch,l,r,s;
	double totloops;
	double num_trans;
	int conjdis;
	int *existential;
	int plot_ok, plotd_ok;

	// Return variable: 1 if controller is anything but empty, -1 if empty
	int ctrl=1;

	char *bddcont1;
	char *bddcont2;
	char *bddrescont;
	char *bddresW;
	char *bddcont1_name;
	char *bddcont2_name;
	char *bddrescont_name;
	char *bddresW_name;
	char suffix[5];

	//BDD Variables		
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	/* Check for proper number of arguments. */
	if(nrhs<6 || nrhs>7) {
		mexErrMsgTxt("Six inputs required:\n controller filename, controller 2 filename, controller result filename, control set result filename, 1 (AND) or 0 (OR), and type flag.\n Optional: Verbose Level (defaulty 0).");
		return;
	} 
	if(nrhs==6)
		verbose=(int)mxGetScalar(prhs[6]);
	else
		verbose=0;

	conjdis=(int)mxGetScalar(prhs[4]);
	type_flag=(int)mxGetScalar(prhs[5]);

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

	mexEvalString("ntransits=prod([params_symb.num]+ones(params_symb.n+params_symb.m,1));");
	num_trans=mxGetScalar(mexGetVariable("caller","ntransits"));

	bddcont1=mxArrayToString(prhs[0]);
	bddcont2=mxArrayToString(prhs[1]);
	bddrescont=mxArrayToString(prhs[2]);
	bddresW=mxArrayToString(prhs[3]);

	bddcont1_name=(char*)mxMalloc(strlen(bddcont1)+5);

	bddcont2_name=(char*)mxMalloc(strlen(bddcont2)+5);

	bddrescont_name=(char*)mxMalloc(strlen(bddrescont)+5);

	bddresW_name=(char*)mxMalloc(strlen(bddresW)+5);

	strcpy(suffix,".bdd");

	strcpy(bddcont1_name, bddcont1);
	strcat(bddcont1_name,suffix);
	mxFree(bddcont1);

	strcpy(bddcont2_name, bddcont2);
	strcat(bddcont2_name,suffix);
	mxFree(bddcont2);

	strcpy(bddrescont_name, bddrescont);
	strcat(bddrescont_name,suffix);
	mxFree(bddrescont);

	strcpy(bddresW_name, bddresW);
	strcat(bddresW_name,suffix);
	mxFree(bddresW);

	//check for file existence
	FILE * T1File;
	T1File = fopen(bddcont1_name,"r");

	if (T1File==NULL)
	{
		mexEvalString("warndlg('Controller1 file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Controller1 file not found!.\n");
	}
	else
	{
		fclose(T1File);
	}

	//check for file existence
	FILE * T2File;
	T2File = fopen(bddcont2_name,"r");

	if (T2File==NULL)
	{
		mexEvalString("warndlg('Controller2 file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Controller2 file not found!.\n");
	}
	else
	{
		fclose(T2File);
	}

	//BDD INITIALIZATIONS		

	// We want all the BDD's to be the same number of variables
	numVars=params_symb.totbits;  // num of vars; if unknown set to 0
	numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
	numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
	cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
	//maxCacheSize=10485760*2;   // default for CUDD package

	ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, 0); //maxCacheSize);

 	ddnodearray = (DdNode**)malloc(numroots_design * sizeof(DdNode*));
 	for (r=0; r<numroots_design; r++) 
		ddnodearray[r] = Cudd_ReadLogicZero(ddman);
	
	ddnodearray[BDD_T1]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddcont1_name, NULL);
	ddnodearray[BDD_T2]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddcont2_name, NULL);

	existential = new int[params_symb.totbits];

	//initialize everything to 2
	for (s=0; s<params_symb.totbits; s++)
		existential[s] = 2;

	//final x and u set to 1 -- rest is ignored
	for (l=params_symb.totbits-params_symb.nbitsloop; l<params_symb.totbits; l++)
		existential[l] = 1;

	DdNode* cube_array;

	cube_array = Cudd_CubeArrayToBdd(ddman, existential);
	Cudd_Ref(cube_array);
	
	if(conjdis)
		ddnodearray[BDD_TR]=Cudd_bddAnd(ddman,ddnodearray[BDD_T1],ddnodearray[BDD_T2]);
	else
		ddnodearray[BDD_TR]=Cudd_bddOr(ddman,ddnodearray[BDD_T1],ddnodearray[BDD_T2]);

	Cudd_Ref(ddnodearray[BDD_TR]);
	ddnodearray[BDD_WR] = Cudd_bddExistAbstract(ddman, ddnodearray[BDD_TR], cube_array);
	Cudd_Ref(ddnodearray[BDD_WR]);
	Cudd_RecursiveDeref(ddman,cube_array);

	if(Cudd_IsConstant(ddnodearray[BDD_TR]) && Cudd_IsComplement(ddnodearray[BDD_TR]))
	{
		if(type_flag == 1)
		{
			//mexPrintf("WARNING: The target set controller has an empty domain of definition.\n");
			//mexEvalString("warndlg('The target set controller has an empty domain of definition.','Pessoa Warning')");
			ctrl = -1;
		}
		else if(type_flag == 2)	
		{
			//mexPrintf("WARNING: The constraint set controller has an empty domain of definition.\n");
			//mexEvalString("warndlg('The constraint set controller has an empty domain of definition.','Pessoa Warning')");
			ctrl = -1;
		}
		else
		{
			//mexPrintf("WARNING: The controller has an empty domain of definition.\n");
			//mexEvalString("warndlg('The controller has an empty domain of definition.','Pessoa Warning')");
			ctrl = -1;		
		}
	}

	if(Cudd_IsConstant(ddnodearray[BDD_WR]) && Cudd_IsComplement(ddnodearray[BDD_WR]))
	{
		if(type_flag == 1)
		{
			//mexPrintf("WARNING: The target set controller domain is empty.\n");
			//mexEvalString("warndlg('The target set controller domain is empty.','Pessoa Warning')");
			ctrl = -1;
		}
		else if(type_flag == 2)	
		{
			//mexPrintf("WARNING: The constraint set controller domain controller is empty.\n");
			//mexEvalString("warndlg('The constraint set controller domain is empty.','Pessoa Warning')");
			ctrl = -1;
		}
		else
		{
			//mexPrintf("WARNING: The controller domain is empty.\n");
			//mexEvalString("warndlg('The controller domain is empty.','Pessoa Warning')");
			ctrl = -1;		
		}
	}

	//Dump BDD to file for viewing in DOT format
	if(verbose==3){
		fp=fopen("dumpfile_combine.dot", "w");
		Cudd_DumpDot(ddman, numroots_design, ddnodearray,NULL,NULL,fp); 
		fclose(fp);
	}

	// Plotting
	if(verbose==2)
	{
		if(ctrl==1)
		{
			mexPrintf("\nPlotting controller domain...\n");
			plotd_ok = PlotSet(ddman,ddnodearray,BDD_WR, &params_symb);
			if(plotd_ok == -1)
				mexPrintf("Plotting controller domain stopped... \n");

			mexPrintf("\nPlotting controller...\n");
			plot_ok = PlotFSM(ddman,ddnodearray,BDD_TR, &params_symb);
			if(plot_ok == -1)
				mexPrintf("Plotting controller stopped... \n");
		}
	}

	// Store the resulting controller
	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[BDD_TR], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddrescont_name, NULL);

	// Store the W set
	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[BDD_WR], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddresW_name, NULL);

    for (s=0; s<numroots_design; s++) 
      Cudd_RecursiveDeref(ddman, ddnodearray[s]);


	Cudd_Quit(ddman);

	if(verbose==3){	
		mexEvalString("profile off");
		mexEvalString("profsave(profile('info'),'profile_results')");}
}
