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
// BDD nodearray number for W(x) 
const int BDD_W = 1;
// BDD with the fixed point algorithm solution
const int BDD_FP = 2;

static DdManager* ddman;     // global ddmanager pointer variable
static DdNode** ddnodearray; // global array of pointers to nodes 

int numroots_design=3;    // number of roots/BDDs/functions

int bdd_main(int type, s_vector* params_symb, int type_flag, long nbatch)
// Function called from the mex function. Computes the controller based on the type (Safety, Reachability, etc).
{ 
	int deter = params_symb->deter;
	int *existential;
	long i;
	
	// Return variable: 1 if controller is anything but empty, -1 if empty
	int ctrl=1;

	struct timeval result, start, finish;

	// CHANGED!: Line below
	DdNode *bdd_safety_set, *bdd_reach_set;
	
	if(type>4 && deter==2)
	{
		type=type-4;
		deter=1;
	}
	
	// Safety
	if(type == 1)
	{
		//get a copy of the safety target set bdd
		bdd_safety_set=ddnodearray[BDD_W];
		Cudd_Ref(bdd_safety_set);
		
		//get init time
		gettimeofday(&start, NULL);

		// This function calls to calculate the fixed point algorithm for safety
		if (deter == 1)
			FPBDD_safety(ddman,ddnodearray,BDD_T,BDD_W,BDD_FP,params_symb);

		if (deter == 2)
			FPBDD_safety_nondeter(ddman,ddnodearray,BDD_T,BDD_W,BDD_FP,params_symb);

		//get final time
		gettimeofday(&finish, NULL);
		timeval_subtract(&result, &finish, &start);

		DdNode* bdd_result;

		// AND the defined set BDD with the BDD defined in the name
		bdd_result=Cudd_bddAnd(ddman,bdd_safety_set,ddnodearray[BDD_W]);
		
		Cudd_Ref(bdd_result);

		if(Cudd_IsConstant(bdd_result) && Cudd_IsComplement(bdd_result))
		{
			if(type_flag == 1)
			{
				mexPrintf("WARNING: The (constraint set) safety controller has an empty domain.\n");
				mexEvalString("warndlg('The (constraint set) safety controller has an empty domain.','Pessoa Warning')");
				ctrl = -1;
			}
			else if(type_flag == 2)	
			{
				mexPrintf("WARNING: The (target set) safety controller has an empty domain.\n");
				mexEvalString("warndlg('The (target set) safety controller has an empty domain.','Pessoa Warning')");
				ctrl = -1;
			}
			else
			{
				mexPrintf("WARNING: The safety controller has an empty domain.\n");
				mexEvalString("warndlg('The safety controller has an empty domain.','Pessoa Warning')");
				ctrl = -1;		
			}
		}
		else if(Cudd_bddLeq(ddman, bdd_result, bdd_safety_set) && Cudd_bddLeq(ddman, bdd_safety_set, bdd_result))
		{
			if(type_flag == 1)
			{
				mexPrintf("The (constraint set) safety controller domain is the whole constraint set .\n");
				//mexEvalString("warndlg('The controller domain is the whole target set .','Pessoa Warning')");
				ctrl = 1;
			}
			else if(type_flag == 2)	
			{
				mexPrintf("The (target set) safety controller domain is the whole target set.\n");
				//mexEvalString("warndlg('The controller domain is the whole safety set.','Pessoa Warning')");
				ctrl = 1;
			}
			else
			{
				mexPrintf("The safety controller domain is the whole target set.\n");
				//mexEvalString("warndlg('The safety controller domain is the whole target set.','Pessoa Warning')");
				ctrl = 1;		
			}
		}
		else
		{
			if(type_flag == 1)
			{
				mexPrintf("WARNING: The (constraint set) safety controller domain is NOT\n");
				mexPrintf("         the whole constraint set.\n");
				mexEvalString("warndlg('The (constraint set) safety controller domain is NOT the whole constraint set.','Pessoa Warning')");
				ctrl = 1;
			}
			else if(type_flag == 2)	
			{
				mexPrintf("WARNING: The (target set) safety controller domain is NOT\n");
				mexPrintf("         the whole target set.\n");
				mexEvalString("warndlg('The (target set) safety controller domain is NOT the whole target set.','Pessoa Warning')");
				ctrl = 1;
			}
			else
			{
				mexPrintf("WARNING: The safety controller domain is NOT the whole target set.\n");
				mexEvalString("warndlg('The safety controller domain is NOT the whole target set.','Pessoa Warning')");
				ctrl = 1;		
			}
		}

		Cudd_RecursiveDeref(ddman, bdd_result);
		Cudd_RecursiveDeref(ddman, bdd_safety_set);
	}

	// Reachability
	if(type == 2)
	{
		//get init time
		gettimeofday(&start, NULL);

		//CHANGE!: Added these lines
		bdd_reach_set=ddnodearray[BDD_W];
		Cudd_Ref(bdd_reach_set);
		// Til here

		// This function calls to calculate the fixed point algorithm for reachability
		if (deter == 1)
			FPBDD_reachability(ddman,ddnodearray,BDD_T,BDD_W,BDD_FP,params_symb);

		if (deter == 2)
			FPBDD_reachability_nondeter(ddman,ddnodearray,BDD_T,BDD_W,BDD_FP,params_symb);

		//get final time
		gettimeofday(&finish, NULL);
		timeval_subtract(&result, &finish, &start);

		// existential for u and x' only
		existential = new int[params_symb->totbits];
	
		//u and x' to 1 -- x marked to 2
		for (i=0; i<params_symb->nbitsx; i++)
			existential[i] = 2;

		for (i=params_symb->nbitsx; i<params_symb->totbits; i++)
			existential[i] = 1;

		DdNode* cube_array_exist, *bdd_xset, *bdd_not_Wset, *bdd_result;
		cube_array_exist = Cudd_CubeArrayToBdd(ddman, existential);
		Cudd_Ref(cube_array_exist);

		//Compute full domain of transition system
	
		bdd_xset = Cudd_bddExistAbstract(ddman,ddnodearray[BDD_T] , cube_array_exist);
		Cudd_Ref(bdd_xset);

		// CHANGE!: Next block of code
		// NOT target set
		bdd_not_Wset=Cudd_bddNand(ddman, bdd_reach_set, Cudd_ReadOne(ddman));
		Cudd_Ref(bdd_not_Wset);
		Cudd_RecursiveDeref(ddman, bdd_reach_set);
		// Til here

		// AND the x_set BDD with NOT of target set = Potential full domain
		bdd_result=Cudd_bddAnd(ddman,bdd_xset,bdd_not_Wset);
		Cudd_Ref(bdd_result);
		
		//CHANGE!: Two recursive derefs Moved here:
		Cudd_RecursiveDeref(ddman, bdd_xset);
		Cudd_RecursiveDeref(ddman, bdd_not_Wset);


		//CHANGE!: The check for emptyness is on the DOMAIN of the controller
		if(Cudd_IsConstant(ddnodearray[BDD_W]) && Cudd_IsComplement(ddnodearray[BDD_W]))
		{
				mexPrintf("WARNING: The reachability controller has an empty domain.\n");
				mexEvalString("warndlg('The reachability controller has an empty domain.','Pessoa Warning')");
				ctrl = -1;		
		//CHANGE!: the check of equality is between the domain of the controller and (x_set AND (NOT W))
		}else if(Cudd_bddLeq(ddman, bdd_result, ddnodearray[BDD_W]) && Cudd_bddLeq(ddman, ddnodearray[BDD_W], bdd_result))
		{
				mexPrintf("The reachability controller domain is the whole state set.\n");
				//mexEvalString("warndlg('The reachability controller domain is the whole state set.','Pessoa Warning')");
				ctrl = 1;		
		}
		else
		{
				mexPrintf("WARNING: The reachability controller domain is NOT the whole state set.\n");
				mexEvalString("warndlg('The reachability controller domain is NOT the whole state set.','Pessoa Warning')");
				ctrl = 1;		
		}

		Cudd_RecursiveDeref(ddman, bdd_result);
		//CHANGE!: I moved Two recursive_derefs UP
	}
	
  return ctrl; 
}

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), calls "bdd_main()", decides which control problem is asked for and solves it.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	s_vector params_symb;
	int mxout = 0;
	int ret, ok;
	int type, type_flag, verbose;
	FILE *fp;
   	mxArray *psv, *psv_w;
	long nbatch,r,s;
	double totloops;
	double total_time;
	double num_trans;
	long memuse1, memuse2, memuse3, mem_res;
	int plot_ok, plotd_ok;

	char *bddsys;
	char *bddW;
	char *bddcont;
	char *bddset;
	char suffix[5];
	char *bddsys_name;
	char *bddW_name;
	char *bddcont_name;
	char *bddset_name;

	//BDD Variables	
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	/* Check for proper number of arguments. */
	if(nrhs<6 || nrhs>7) {
		mexErrMsgTxt("Six inputs required:\n system filename, target set filename, controller result filename, target set result filename, type (1-safety or 2-reachabililty), and type flag. \n Optional: Verbose Level (defaulty 0).");
		return;
	} 

	if(nrhs==6)
		verbose=(int)mxGetScalar(prhs[6]);
	else
		verbose=0;

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

	type=(int)mxGetScalar(prhs[4]);
	type_flag=(int)mxGetScalar(prhs[5]);

	bddsys=mxArrayToString(prhs[0]);
	bddW=mxArrayToString(prhs[1]);
	bddcont=mxArrayToString(prhs[2]);
	bddset=mxArrayToString(prhs[3]);

	strcpy(suffix,".bdd");

	bddsys_name=(char*)mxMalloc(strlen(bddsys)+5);
	strcpy(bddsys_name, bddsys);
	strcat(bddsys_name,suffix);
	mxFree(bddsys);

	//check for file existence
	FILE * smFile;
	smFile = fopen(bddsys_name,"r");

	if (smFile==NULL)
	{
		mexEvalString("warndlg('Symbolic model file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Symbolic model file not found!.\n");
	}
	else
	{
		fclose(smFile);
	}

	bddW_name=(char*)mxMalloc(strlen(bddW)+5);
	strcpy(bddW_name, bddW);
	strcat(bddW_name,suffix);
	mxFree(bddW);
	
	//check for file existence
	FILE * tsFile;
	tsFile = fopen(bddW_name,"r");

	if (tsFile==NULL)
	{
		mexEvalString("warndlg('Target set file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Target set file not found!.\n");
	}
	else
	{
		fclose(tsFile);
	}

	bddcont_name=(char*)mxMalloc(strlen(bddcont)+5);
	strcpy(bddcont_name, bddcont);
	strcat(bddcont_name,suffix);
	mxFree(bddcont);

	bddset_name=(char*)mxMalloc(strlen(bddset)+5);
	strcpy(bddset_name, bddset);
	strcat(bddset_name,suffix);
	mxFree(bddset);

	//BDD Static INITIALIZATIONS 	

	// We want all the BDD's to be the same number of variables
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

	ddnodearray[BDD_T]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddsys_name, NULL);

	// Plot first some info about the simulations:
	if(verbose>0)
	{
		mexEvalString("t_start = tic;");
		if(verbose==3)
			mexEvalString("profile on");
	}

	ddnodearray[BDD_W]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddW_name, NULL);

	memuse2=Cudd_ReadMemoryInUse(ddman);

	//Call the main function to manipulate the BDD's (Reachability, Safety)
	ret=bdd_main(type,&params_symb,type_flag,nbatch);

	memuse3=Cudd_ReadMemoryInUse(ddman);

	//Dump BDD to file for viewing in DOT format
	if(verbose==3){
		fp=fopen("dumpfile_desgn.dot", "w");
		Cudd_DumpDot(ddman, numroots_design, ddnodearray,NULL,NULL,fp); 
		fclose(fp);
	}

	// Plotting
	if(verbose==2)
	{
		if(ret == 1)
		{

			mexPrintf("\nPlotting controller domain...\n");
			plotd_ok = PlotSet(ddman,ddnodearray,BDD_W, &params_symb);
			if(plotd_ok == -1)
				mexPrintf("Plotting controller domain stopped... \n");

			if(type == 1)
			{
				mexPrintf("\nPlotting safety controller...\n");
			}
			if(type == 2)	
			{
				mexPrintf("\nPlotting reachability controller...\n");
			}

			plot_ok = PlotFSM(ddman,ddnodearray,BDD_FP,&params_symb);
			if(plot_ok == -1)
					mexPrintf("Plotting controller stopped... \n");


		}
	}

	// Store the controller BDD
	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[BDD_FP], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddcont_name, NULL);

	if(ok)
		mexPrintf("Symbolic controller sucessfully saved to '.bdd' file.\n");
	else
		mexErrMsgTxt("\nERROR: Symbolic controller failed to save.\n");


	// Store the W set
	ok = Dddmp_cuddBddStore(ddman, NULL, ddnodearray[BDD_W], NULL,
					NULL, DDDMP_MODE_TEXT, DDDMP_VARIDS,
					bddset_name, NULL);
	if(ok)
		mexPrintf("Symbolic controller domain sucessfully saved to '.bdd' file.\n");
	else
		mexErrMsgTxt("ERROR: Symbolic controller domain failed to save.\n");


	if(verbose>0){
		// Print some extra-info
		mexPrintf("Memory in use: %ld",memuse1);
		mexPrintf(" bytes (BDD manager); \n");
		mem_res = memuse3-memuse2;
		if(mem_res<0)
			mem_res = 0;
		mexPrintf("               %ld",mem_res);
		mexPrintf(" bytes (symbolic controller). \n");
	}

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

	mexPrintf("\n---------------- Pessoa: Controller Synthesis Terminated --------------- \n");
}
