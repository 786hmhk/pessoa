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

// BDD nodearray number for user-defined rectangle
const int BDD_set = 0;
// BDD nodearray number for controller domain
const int BDD_ctrldomain = 1;
// BDD nodearray number for AND of the two above
const int BDD_result = 2;

int numroots_design=3;  // number of roots/BDDs/functions

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), compares two sets (returns 0 if there is no intersection between sets, returns 1  if the second set is contained in the first set, and returns 2 if there is partial intersection.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

	mexEvalString("t_start = tic;");

	s_vector params_symb;
	int ret, ok;
	int type;
	double *wmin, *wmax, *output;
	double total_time;
	long r,s;
	FILE *fp;
   	mxArray *psv;

	char *bddctrl;
	char suffix[5];
	char *bddctrl_name;
	char *bddset;
	char *bddset_name;

	//BDD Variables		
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	// output
	plhs[0] = mxCreateDoubleScalar(2);
	output = mxGetPr(plhs[0]);
	
	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>2) {
		mexErrMsgTxt("Two inputs required: names of two compatible sets (domain of a controller and an arbitrarty set).");
		return;
	} 

	//Copy data from Matlab workspace 
	psv=mexGetVariable("caller","params_symb");
	
	//Copy data to variables
	params_symb.n=(int)mxGetScalar(mxGetField(psv,0,"n"));
	params_symb.m=(int)mxGetScalar(mxGetField(psv,0,"m"));
	params_symb.nume=(double *)mxGetPr(mxGetField(psv,0,"nume"));   
	params_symb.totbits=(int)mxGetScalar(mxGetField(psv,0,"totbits"));
	params_symb.nbitsloop=(int)mxGetScalar(mxGetField(psv,0,"nbitsloop"));
	params_symb.nbits=(double *)mxGetPr(mxGetField(psv,0,"nbits"));
	params_symb.deter=(int)mxGetScalar(mxGetField(psv,0,"deter"));
	params_symb.nbitsx=(int)mxGetScalar(mxGetField(psv,0,"nbitsx"));
	
	bddctrl=mxArrayToString(prhs[0]);

	strcpy(suffix,".bdd");

	bddctrl_name=(char*)mxMalloc(strlen(bddctrl)+5);
	strcpy(bddctrl_name, bddctrl);
	strcat(bddctrl_name,suffix);
	mxFree(bddctrl);

	// BDD name of set
	bddset=mxArrayToString(prhs[1]);

	bddset_name=(char*)mxMalloc(strlen(bddset)+5);
	strcpy(bddset_name, bddset);
	strcat(bddset_name,suffix);
	mxFree(bddset);

	//check for file existence
	FILE * ctrlFile;
	ctrlFile = fopen(bddctrl_name,"r");

	if (ctrlFile==NULL)
	{
		mexEvalString("warndlg('Pessoa Compare: Controller domain file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Pessoa Compare: Controller domain file not found!.\n");
	}
	else
	{
		fclose(ctrlFile);
	}

	//check for file existence
	FILE * csFile;
	csFile = fopen(bddset_name,"r");

	if (csFile==NULL)
	{
		mexEvalString("warndlg('Pessoa Compare: Set file not found!','Pessoa Error')");
		mexErrMsgTxt("\nERROR: Pessoa Compare: Set file set file not found!.\n");
	}
	else
	{
		fclose(csFile);
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

	// Load the domain set
	ddnodearray[BDD_ctrldomain]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddctrl_name, NULL);

	// Load the set to compare with
	ddnodearray[BDD_set]=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
				DDDMP_MODE_DEFAULT, bddset_name, NULL);

	// AND the defined set BDD with the BDD defined in the name
	ddnodearray[BDD_result]=Cudd_bddAnd(ddman,ddnodearray[BDD_set],ddnodearray[BDD_ctrldomain]);
	Cudd_Ref(ddnodearray[BDD_result]);

	if(Cudd_IsConstant(ddnodearray[BDD_result]) && Cudd_IsComplement(ddnodearray[BDD_result]))
	{
			mexPrintf("\nThere is no intersection between the sets.\n");
			output[0] = 0;
	}else if(Cudd_bddLeq(ddman, ddnodearray[BDD_result], ddnodearray[BDD_set]) && Cudd_bddLeq(ddman, ddnodearray[BDD_set], ddnodearray[BDD_result]))
	{
		mexPrintf("\nThe second set is contained in the first set (domain).\n");
		output[0] = 1;
	}
	else
	{
		mexPrintf("\nThere is a partial intersection between the first  \n");
		mexPrintf(" set (domain) and the second one.\n");
		output[0] = 2;	
	}

	for (s=0; s<numroots_design; s++) 
	      Cudd_RecursiveDeref(ddman, ddnodearray[s]);

	Cudd_Quit(ddman);


	mexEvalString("t_end = toc(t_start);");
	total_time = mxGetScalar(mexGetVariable("caller","t_end"));
	if(total_time < 0.001)
		mexPrintf("Elapsed time: 0.001 seconds. \n");
	else		
	{	
		mexPrintf("Elapsed time: %.3f", total_time);
		mexPrintf(" seconds. \n", total_time);
	}
	mexPrintf("\n----------------- Pessoa: Comparing Sets Terminated -------------------- \n");
}
