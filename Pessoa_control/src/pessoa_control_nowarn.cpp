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
#include <limits>

using namespace std;

static DdManager* ddman;     // global ddmanager pointer variable
static DdNode** ddnodearray; // global array of pointers to nodes 

static DdNode* controldd;
static DdNode* existbdd;

static int firstrun=1;
static char contname[30];

static void CloseBDDman(void)
// Call back function for the event of "clear mex". Forces the unload from memory of the BDD's in use as well as the BDD manager.
{
			Cudd_RecursiveDeref(ddman,existbdd);
			Cudd_RecursiveDeref(ddman,controldd);	
			Cudd_Quit(ddman);
			firstrun=1;
}

/************************************

 MEX-MAIN function "mexFunction"
 Initializes data (params_symb structure and pointer to params_symb), in charge of controller execution.

*************************************/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

	s_vector params_symb;
	double *arrayu, *arrayx, *arrayuold;
 	mxArray *psv;
	int k,ind,ret;
	double *nbits;
	long i,j;
	int error=0;
	
	char *contfile;
	char suffix[5];
	char *contfile_name;


	//BDD Variables		
	short numVars;  // num of vars; if unknown set to 0
	short numVarsZ; // num of vars for ZBDDs; if unknown set to 0
	int numSlots; // default for CUDD package
	int cacheSize; // default for CUDD package
	//int maxCacheSize;   // default for CUDD package

	/* Check for proper number of arguments. */
	if(nrhs<3 || nrhs>3) {
		mexErrMsgTxt("Three inputs required (controller name, last input, state).");
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

	// Existvect for existential elimination of x and x' in [x u x'].
	int existvect[params_symb.totbits];
	
	//To open or check if it changed the file with the controller
	contfile=mxArrayToString(prhs[0]);

	strcpy(suffix,".bdd");

	contfile_name=(char*)mxMalloc(strlen(contfile)+5);
	strcpy(contfile_name, contfile);
	strcat(contfile_name,suffix);
	mxFree(contfile);


	if(firstrun || strcmp(contname,contfile_name))
	{
	
		firstrun--;
		mexAtExit(CloseBDDman);

		//Initialize BDD manager
		if (strlen(contfile_name)>30)
		{
			strncpy(contname,contfile_name,29);
			contname[29]='\0';
		}else
		{
			strcpy(contname,contfile_name);
		}	
		if(existbdd!=NULL)
			Cudd_RecursiveDeref(ddman,existbdd);
		
		if(controldd!=NULL)
			Cudd_RecursiveDeref(ddman,controldd);	
		
		// We want all the BDD's to be the same number of variables
		numVars=params_symb.totbits;  // num of vars; if unknown set to 0
		numVarsZ=0; // num of vars for ZBDDs; if unknown set to 0
		numSlots=CUDD_UNIQUE_SLOTS; // default for CUDD package
		cacheSize=CUDD_CACHE_SLOTS; // default for CUDD package
		//maxCacheSize=10485760*2;   // default for CUDD package

		ddman = Cudd_Init(numVars, numVarsZ, numSlots, cacheSize, 0); //maxCacheSize);

		//check for file existence
		FILE * pFile;
		pFile = fopen (contfile_name,"r");

		if (pFile==NULL)
		{
			mexEvalString("warndlg('Controller file not found!','Pessoa Error')");
			mexErrMsgTxt("\nERROR: Controller file not found!.\n");
		}
		else
		{
			fclose(pFile);
		}
		
		controldd=Dddmp_cuddBddLoad(ddman, DDDMP_VAR_MATCHIDS, NULL, NULL, NULL,
					DDDMP_MODE_DEFAULT, contfile_name, NULL);

		mxFree(contfile_name);
		
		k=params_symb.totbits-params_symb.nbitsloop;
		for(i=0;i<params_symb.totbits;i++)
			existvect[i]=1;

		nbits = params_symb.nbits;
		for (i = 0; i < params_symb.m; i++) {
			ind=1<<((int)nbits[params_symb.n+i]-1);
			while (ind>0) { 
				existvect[k]=2;
				k++;
				ind >>= 1;
			}
		}
		
		existbdd=Cudd_CubeArrayToBdd(ddman,existvect);
		Cudd_Ref(existbdd);
		
	}

	// Other initializations

	plhs[0] = mxCreateDoubleMatrix(1,params_symb.m,mxREAL);
	arrayx=mxGetPr(prhs[2]);
	arrayuold=mxGetPr(prhs[1]);
	arrayu=mxGetPr(plhs[0]);
	/* Check the state is in XSET */
	
	for(i=0;i<params_symb.n;i++){
		if(arrayx[i]>params_symb.nume[i] || arrayx[i]<0){
			for(j=0;j<params_symb.m;j++){
				arrayu[j] = numeric_limits<float>::quiet_NaN();
//				mexEvalString("beep");
//				mexEvalString("warndlg('The controller has no input defined for current state (outside the domain of the abstraction).','Pessoa Warning')");
				error=1;
			}
		}
	}
	

	/*	Controller execution	*/
	if(!error){
		ret=ValidUInBDD(ddman,ddnodearray,controldd,existbdd,arrayx,arrayuold,arrayu,&params_symb);

		if(!ret){
			for(j=0;j<params_symb.m;j++){
				arrayu[j] = numeric_limits<float>::quiet_NaN();
//				mexEvalString("beep");
//				mexEvalString("warndlg('The controller has no input defined for current state.','Pessoa Warning')");
			}
		}
	}
}
