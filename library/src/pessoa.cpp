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
#include "pessoa.h"

using namespace std;

/***********************************
 
	This function computes a time difference (between start and end time)
	
	INPUTS: result - time difference
		x - initial time
		y - ending time
	OUTPUTS: returns the time difference
	
************************************/
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
	int nsec;

       /* Perform the carry for the later subtraction by updating y. */
       if (x->tv_usec < y->tv_usec) {
         nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
         y->tv_usec -= 1000000 * nsec;
         y->tv_sec += nsec;
       }
       if (x->tv_usec - y->tv_usec > 1000000) {
         nsec = (x->tv_usec - y->tv_usec) / 1000000;
         y->tv_usec += 1000000 * nsec;
         y->tv_sec -= nsec;
       }
     
       /* Compute the time remaining to wait.
          tv_usec is certainly positive. */
       result->tv_sec = x->tv_sec - y->tv_sec;
       result->tv_usec = x->tv_usec - y->tv_usec;
     
       /* Return 1 if result is negative. */
       return x->tv_sec < y->tv_sec;
}

/***********************************
 
	This function takes a transition and checks if it belongs to the BDD.
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		bddnum - the index number for which BDD to check
		array - input string to check against BDD
 		params_symb   - data structure containing discretization
                      information
		nnm - length of the array
	OUTPUTS: 1 if the transition is in the BDD; 0 otherwise
	
************************************/
bool IsInBDD(DdManager* ddmanl,DdNode** ddnodearrayl, int bddnum, double *array, s_vector* params_symb, int nnm)
{
    int totbits = params_symb->totbits;
    long i;
    double val_bool=0;
    int bddarray[totbits];
    DdNode *basecube, *var, *tmp;
    double k=0;

   for(i=0;i<totbits;i++)
	bddarray[i]=1;
	
   basecube = Cudd_ReadOne(ddmanl);
   Cudd_Ref(basecube);

   for (i = 0; i < nnm; i++) 
   {
		int ind=1<<((int)params_symb->nbits[i]-1);
		while (ind>0) 
		{ 
			val_bool=(int)array[i] & ind;
			var = Cudd_bddIthVar(ddmanl,k);
			k++;
			tmp = Cudd_bddAnd(ddmanl,basecube,Cudd_NotCond(var,val_bool==0));
			if (tmp == NULL) {
				Cudd_RecursiveDeref(ddmanl,basecube);
				return(NULL);
			}
			Cudd_Ref(tmp);
			Cudd_RecursiveDeref(ddmanl,basecube);
			basecube = tmp;
			ind >>= 1;
		}
   }

	Cudd_BddToCubeArray(ddmanl, basecube, bddarray);

	Cudd_RecursiveDeref(ddmanl,basecube);

	return(!Cudd_IsComplement(Cudd_Eval(ddmanl,ddnodearrayl[bddnum],bddarray)));

} 


/***********************************
 
	Plots the FSM given by fsm number (Only works for 2D systems)
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array 
		fsm - Finite State Machine to plot
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: -1 if the plotting fails.
	
************************************/
int PlotFSM(DdManager* ddmanl,DdNode** ddnodearrayl, int fsm, s_vector* params_symb)
{
	int totbits = params_symb->totbits;
	int n = params_symb->n;
	int m = params_symb->m;	

	int mxout = 0;
	double *array, *tempv, *p1, *p2;
	long i, j, countin, countout;
        double totstates=1;

	mxArray *pp1, *pp2;
	mxArray *ppf_steps, *pftotstates;
	double *pf_steps, *ftotstates; 

	if(n!=2)
	{
		mexEvalString("warndlg(Can't plot systems other than 2-D!','Pessoa Error')");
		return -1;	
	}

	pp1=mxCreateNumericMatrix(1,n,mxDOUBLE_CLASS,mxREAL);
	p1=(double *) mxGetData(pp1);
	pp2=mxCreateNumericMatrix(1,n,mxDOUBLE_CLASS,mxREAL);
	p2=(double *) mxGetData(pp2);
	int oobs=(1<<((int)params_symb->nbits[0]-1));

	array=(double *)mxMalloc((2*n+m)*sizeof(double));
	tempv=(double *)mxMalloc((2*n+m)*sizeof(double));
	for(i=0;i<2*n+m;i++)
		tempv[i]=0;

	countin=0;
	countout=0;	

       for (j=0;j<(2*n+m);j++)
	    totstates*=params_symb->nume[j]+1;

	ppf_steps = mxCreateDoubleScalar(0);
	pftotstates = mxCreateDoubleScalar(totstates);
	pf_steps=(double *) mxGetData(ppf_steps);
	ftotstates=(double *) mxGetData(pftotstates);
	pf_steps[0]=0;
	ftotstates[0]=totstates;
	
	mexPutVariable("caller", "pftotstates", pftotstates);

	mexEvalString("h=waitbar(0,'Plotting FSM in progress','Name','Pessoa','CreateCancelBtn','global stopbuild;');");

	mexEvalString("figure; hold on");		
	while(1)
	{

			if(IsInBDD(ddmanl,ddnodearrayl,fsm, tempv, params_symb, 2*n+m))
			{
					memcpy(p1, &tempv[0], n*sizeof(double));
					memcpy(p2, &tempv[n+m], n*sizeof(double));

					mexPutVariable("caller", "p1", pp1);
					mexPutVariable("caller", "p2", pp2);

	// TOFIX: In the case of an Error below, the BDD manager will not be freed, causing a memory leak!!
					if (mexEvalString("hold on; pss_vectarrow(params_symb.eta*(p1'+params_symb.min(params_symb.xoind)),params_symb.eta*(p2'+params_symb.min(params_symb.xoind)));"))
					{
						mexEvalString("warndlg('Plot FSM failed!','Pessoa Error')");
						mexEvalString("delete(h);");
						return(-1);
					}
			}

			pf_steps[0]++;

			i=0;			
			while(i<(2*n+m) && (++tempv[i])>params_symb->nume[i])
			{
				tempv[i++]=0;
			}
			if(i==(2*n+m))
				break;

			// Used for waitbar
			mexPutVariable("caller", "pf_steps", ppf_steps);
			mexEvalString("waitbar(pf_steps/pftotstates,h)");
			if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
			{		
				mexEvalString("delete(h);");
				return(-1);
			}

	}

	mexEvalString("delete(h);");

	mxDestroyArray(pp1);
	mxDestroyArray(pp2);
	mxFree(array);
	mxFree(tempv);
}

/***********************************
 
	This function takes a list of transition relations and generates the 
	binary strings codifying the transitions (xini)x(u)x(xend_i). Then 
	appends it to the BDD.
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		bddnum - the index number for which BDD to check
		array - list of strings to check against BDD
 		params_symb   - data structure containing discretization
                      information
		nbatch	- size of the batch (number of transitions to add)
		nnm - length of the array
	OUTPUTS: return 1 if batch is added to BDD successfully
	
************************************/
int AddToBddBatch(DdManager* ddmanl,DdNode** ddnodearrayl, int bddnum, double *array, s_vector* params_symb, long nbatch, int nnm)
{
   DdNode *basecube, *var, *tmp, *tmp2;
   long i, j;
   double val_bool=0;
   double k;	
	
   for (j=0; j<nbatch; j++) 
   {
	
    basecube = Cudd_ReadOne(ddmanl);
    Cudd_Ref(basecube);

    k=0;
    for (i = 0; i < nnm; i++) {
		int ind=1<<((int)params_symb->nbits[i]-1);
		while (ind>0) {
			val_bool=(int)array[j*nnm+i] & ind;
			var = Cudd_bddIthVar(ddmanl,k);
			k++;
			tmp = Cudd_bddAnd(ddmanl,basecube,Cudd_NotCond(var,val_bool==0));
			if (tmp == NULL) {
				Cudd_RecursiveDeref(ddmanl,basecube);
				return(NULL);
			}
			Cudd_Ref(tmp);
			Cudd_RecursiveDeref(ddmanl,basecube);
			basecube = tmp;
			ind >>= 1;
		}
	}

	tmp2 = Cudd_bddOr(ddmanl, ddnodearrayl[bddnum], basecube);
	Cudd_Ref(tmp2);
	Cudd_RecursiveDeref(ddmanl,basecube);
	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[bddnum]);
	ddnodearrayl[bddnum] = tmp2;

	}

	return(1);

}  

/***********************************

	This function constructs a BDD representation for
	the FSM approximatelly bisimilar to the control system
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		BDD_Tl - index for the Transition BDD in ddnodearray
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: returns 1 if the Transition BDD is built successfully; 0 otherwise
	
************************************/
int buildBDDT(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_Tl, s_vector* params_symb, long nbatch, double totloops)
{
	int n = params_symb->n;
	int m = params_symb->m;

	double *array;
	long i=0;
	mxArray *ptrans;
	int ok;
	struct timeval result, start, finish;
	const mwSize *sizevect;

	double ndone;	
	int ndisp = 0;

	mxArray *psteps;
	psteps = mxCreateDoubleScalar(i);
	double *steps=(double *) mxGetData(psteps);

	//get init time
	gettimeofday(&start, NULL);	
	
	mexEvalString("xnext=zeros(params_symb.n,1);");	mexEvalString("unext=zeros(params_symb.m,1);");
	
	mexEvalString("h=waitbar(0,'Abstraction in progress','Name','Pessoa','CreateCancelBtn','global stopbuild;');");

	for(i=0;i<totloops;i++)
	{
		
		if(mexEvalString("[trans xnext unext]=pss_sim_transit(xnext,unext,nbatch,params_symb,params_cont);"))
		{
			mexEvalString("warndlg('pss_sim_transit.m failed!','Pessoa Error')");
			mexEvalString("delete(h);");
			return(-1);
		}

		ptrans=mexGetVariable("caller","trans");
		array=mxGetPr(ptrans);
		int ndim=(int)mxGetNumberOfDimensions(ptrans);
		if(ndim==2){
			sizevect=mxGetDimensions(ptrans);
			nbatch=(long)sizevect[1];}
		else if(ndim==1)
			nbatch=1;
		else
			nbatch=0;
			
		if(nbatch)
		{

			ok=AddToBddBatch(ddmanl,ddnodearrayl,BDD_Tl, array, params_symb, nbatch, 2*n+m);
		}
		// Used for waitbar
		steps[0] = i+1;
		mexPutVariable("caller", "steps", psteps);
		mexEvalString("waitbar(steps/totloops,h)");
		if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
		{		
			mexEvalString("delete(h); clear global stopbuild;");
			return(-1);
		}

	}

	mexEvalString("delete(h);");

	//get final time
	gettimeofday(&finish, NULL);

	timeval_subtract(&result, &finish, &start);

	return(1);
}

/***********************************
 
	This function computes the control input "u" at some given 
	state "x" (in lazy fashion, i.e. if last u is still ok, then use that u).
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		controldd - CUDD node for the controller
		existdd - CUDD node for existential check
		arrayx - state x
		arrayuold - input at first
		arrayu - the new input
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: returns 1 if the input u is in BDD; 0 otherwise
	
************************************/
int ValidUInBDD(DdManager* ddmanl,DdNode** ddnodearrayl,DdNode* controldd, DdNode* existbdd, double *arrayx, double *arrayuold, double *arrayu, s_vector* params_symb)
{
	int n = params_symb->n;
	int m = params_symb->m;
	int totbits = params_symb->totbits;
	int nbitsloop = params_symb->nbitsloop;

	int indu, lenlit;
	double val_bool=0;
	int bddarray[totbits], utemp[totbits];
	DdNode *absbdd, *ubdd;
	long i,j, k;
	
	for(i=0;i<totbits;i++)
		bddarray[i]=2;

	k=0;
    	for (i = 0; i < n; i++) {
		int ind=1<<((int)params_symb->nbits[i]-1);
		while (ind>0) { 
			bddarray[k]=((int)arrayx[i] & ind)>0;
			k++;
			ind >>= 1;
		}
	}
	indu=k;
	
	absbdd=Cudd_bddAndAbstract(ddmanl, controldd, Cudd_CubeArrayToBdd(ddmanl,bddarray), existbdd);
	Cudd_Ref(absbdd);

	// First check for the current input ("lazy control")
		k=indu;
		for (i = 0; i < m; i++) {
			int ind=1<<((int)params_symb->nbits[n+i]-1);
			while (ind>0) { 
				bddarray[k]=((int)arrayuold[i] & ind)>0;
				k++;
				ind >>= 1;
			}
		}

		if(!Cudd_IsComplement(Cudd_Eval(ddmanl,absbdd,bddarray)))
		{
			memcpy(arrayu, arrayuold, m*sizeof(double));
			Cudd_RecursiveDeref(ddmanl,absbdd);
			return(1);
		}
		//ubdd=Cudd_LargestCube(ddmanl,absbdd,&lenlit);
		ubdd=Cudd_ShortestPath(ddmanl,absbdd,NULL,NULL,&lenlit);
		Cudd_Ref(ubdd);
		if(Cudd_BddToCubeArray(ddmanl,ubdd,utemp)>0)
		{
			k=totbits-nbitsloop; // k starts at the 1st (leftmost) bit of u
			for(i=0;i<m;i++){
				arrayu[i]=0;
				for(j=params_symb->nbits[n+i];j>0;j--) { 
					if(utemp[k]==1)
						arrayu[i]+=utemp[k]<<(j-1);
					k++;
				}
			}

			Cudd_RecursiveDeref(ddmanl,absbdd);
			Cudd_RecursiveDeref(ddmanl,ubdd);
			return(1);
		}
		else
		{
			Cudd_RecursiveDeref(ddmanl,absbdd);
			Cudd_RecursiveDeref(ddmanl,ubdd);			
			return(0);
		}
} 

/***********************************
 
	Plots the Set given by fsm number (Only works for 2D sets)
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		fsm - set to plot
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: -1 if the plotting fails
	
************************************/
int PlotSet(DdManager* ddmanl,DdNode** ddnodearrayl, int fsm, s_vector* params_symb)
{

	int totbits = params_symb->totbits;
	int n = params_symb->n;
	int m = params_symb->m;

	double *array, *tempv, *p1;
	double *ps_steps, *stotstates;
	long i,j;
	mxArray *pp1;
	double totstates=1;
	
	if(n!=2)
	{
		mexEvalString("warndlg(Can't plot systems other than 2-D!','Pessoa Error')");
		return -1;
	}


	pp1=mxCreateNumericMatrix(1,n,mxDOUBLE_CLASS,mxREAL);
	p1=(double *) mxGetData(pp1);
	
	array=(double *)malloc(n*sizeof(double));
	tempv=(double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		tempv[i]=0;
	
	for (j=0;j<n;j++)
		totstates*=params_symb->nume[j]+1;

	mxArray *pps_steps, *pstotstates;
	pps_steps = mxCreateDoubleScalar(0);
	pstotstates = mxCreateDoubleScalar(totstates);
	ps_steps=(double *) mxGetData(pps_steps);
	stotstates=(double *) mxGetData(pstotstates);
	ps_steps[0]=0;
	stotstates[0]=totstates;

	mexPutVariable("caller", "pstotstates", pstotstates);

	mexEvalString("h=waitbar(0,'Plotting set in progress','Name','Pessoa','CreateCancelBtn','global stopbuild;');");

	mexEvalString("figure; hold on");		
	while(1)
	{

			if(IsInBDD(ddmanl,ddnodearrayl,fsm, tempv, params_symb, n))
			{
					memcpy(p1, &tempv[0], n*sizeof(double));

					mexPutVariable("caller", "p1", pp1);
					mexEvalString("hold on; pp1=params_symb.eta*(p1'+params_symb.min(params_symb.xoind)); plot(pp1(1),pp1(2),'*');");

			}

			ps_steps[0]++;

			i=0;
			while(i<n && (++tempv[i])>params_symb->nume[i])
			{
				tempv[i++]=0;
			}
			if(i==n)
			{
				break;
			}

			// Used for waitbar
			mexPutVariable("caller", "ps_steps", pps_steps);
			mexEvalString("waitbar(ps_steps/pstotstates,h)");
			if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
			{		
				mexEvalString("delete(h);");
				return(-1);
			}


	}

	mexEvalString("delete(h);");

	mxDestroyArray(pp1);
}

/***********************************

	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		BDD_W - index for the W BDD in ddnodearray
		wmin - the minimum point of W
		wmax - the maximum point of W
 		params_symb   - data structure containing discretization
                      information
		nbatch	- size of the batch (number of transitions to add)
	OUTPUTS: returns the total number of states in target set (W BDD); -1 if failed
	
************************************/
long buildBDDW(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_W, double *wmin, double *wmax, s_vector* params_symb, long nbatch)
{
	int n = params_symb->n;

	double *array, *tempv, *wsteps;
	int ok, end=0;
	long i, j=0;
	struct timeval result, start, finish;
	long total_states = 0;
	mxArray *pwsteps;

	//get init time
	gettimeofday(&start, NULL);
	
	array=(double *)mxMalloc(nbatch*n*sizeof(double));
	tempv=(double *)mxMalloc(n*sizeof(double));
	memcpy(tempv, wmin, n*sizeof(double));

	mexEvalString("h=waitbar(0,'Symbolic set construction in progress','Name','Pessoa','CreateCancelBtn','global stopbuild');");

	pwsteps = mxCreateDoubleScalar(j);
	wsteps=(double *) mxGetData(pwsteps);
	wsteps[0]=0;
	
	while(end==0)
	{
		for(j=0; j<nbatch; j++)
		{
			memcpy(&array[j*n], tempv, n*sizeof(double));
			i=0;
			while(i<n && (++tempv[i])>wmax[i])
			{
				tempv[i++]=wmin[i];
			}
			if(i==n)
			{
				end=1;
				j++;
				break;
			}
		}

		ok=AddToBddBatch(ddmanl,ddnodearrayl,BDD_W, array, params_symb, j, n);

		total_states = total_states + j;

		// Used for waitbar
		wsteps[0]++;
		mexPutVariable("caller", "wsteps", pwsteps);
		mexEvalString("waitbar(wsteps/totloops,h)");
		if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
		{		
			mexEvalString("delete(h); clear global stopbuild;");
			return(-1);
		}
	}

	mexEvalString("delete(h);");

	//get final time
	gettimeofday(&finish, NULL);

	timeval_subtract(&result, &finish, &start);

	mxFree(array);
	mxFree(tempv);
	return total_states;
}

/***********************************

	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array
		BDD_W - index for the W BDD in ddnodearray
		wmin - the minimum point of W
		wmax - the maximum point of W
 		params_symb   - data structure containing discretization
                      information
		nbatch - the number of transitions per batch
	OUTPUTS: returns the total number of states in target set (W BDD); -1 if failed
	
************************************/
long build_charfBDDW(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_W, double *wmin, double *wmax, s_vector* params_symb, long nbatch)
{
	int n = params_symb->n;
	int m = params_symb->m;

	double *tempv, *state_array, *set_state_array;
	double *jnbatch, *wfsteps;
	int ok, end=0;
	long i, j=0;
	struct timeval result, start, finish;
	mxArray *pstate_array, *pset_state_array, *pjnbatch;
	long nstates;
	long total_states=0;

	//get init time
	gettimeofday(&start, NULL);
	
	tempv=(double *)mxMalloc(n*sizeof(double));
	state_array=(double *)mxMalloc(nbatch*n*sizeof(double));

	memcpy(tempv, wmin, n*sizeof(double));

	pstate_array=mxCreateNumericMatrix(n,nbatch,mxDOUBLE_CLASS,mxREAL);
	pjnbatch = mxCreateDoubleScalar(nbatch);

	jnbatch=(double *) mxGetData(pjnbatch);

	state_array=(double *) mxGetData(pstate_array);	

	//Initialize
	mexPutVariable("caller", "jnbatch", pjnbatch);

	mexEvalString("h=waitbar(0,'Symbolic set construction in progress','Name','Pessoa','CreateCancelBtn','global stopbuild');");
	
	mxArray *pwfsteps;
	pwfsteps = mxCreateDoubleScalar(j);
	wfsteps=(double *) mxGetData(pwfsteps);
	wfsteps[0]=0;
	
	while(end==0)
	{
		for(j=0; j<nbatch; j++)
		{
			memcpy(&state_array[j*n], tempv, n*sizeof(double));
			i=0;
			while(i<n && (++tempv[i])>wmax[i])
			{
				tempv[i++]=wmin[i];
			}
			if(i==n)
			{
				end=1;
				j++;
				break;
			}

		}

		mexPutVariable("caller", "state_array", pstate_array);

		// TOFIX: In the case of an Error below, the BDD manager will not be freed, causing a memory leak!!
		if(mexEvalString("[set_state_array nstates] = pss_build_charfset(params_symb,state_array,jnbatch);")){
			mexEvalString("warndlg('pss_build_charfset.m failed or missing!','Pessoa Error')");
			mexEvalString("delete(h); clear global stopbuild;");
			mexErrMsgTxt("ERROR: pss_build_charfset.m failed or missing!");
			return(-1);
		}
		
		pset_state_array=mexGetVariable("caller","set_state_array");
		set_state_array=mxGetPr(pset_state_array);

		nstates=mxGetScalar(mexGetVariable("caller","nstates"));

		ok=AddToBddBatch(ddmanl,ddnodearrayl,BDD_W, set_state_array, params_symb, nstates, n);

		total_states = total_states + nstates;

		// Used for waitbar
		wfsteps[0]++;
		mexPutVariable("caller", "wfsteps", pwfsteps);
		mexEvalString("waitbar(wfsteps/totloops,h)");
		if(mexGetVariable("base","stopbuild")!=NULL) // This is a dirty Hack, but I don't know any other way.
		{		
			mexEvalString("delete(h); clear global stopbuild;");
			return(-1);
		}
	}

	mexEvalString("delete(h);");

	//get final time
	gettimeofday(&finish, NULL);

	timeval_subtract(&result, &finish, &start);

	// Matlab will free these variables when exiting
	// mxFree(state_array);
	// mxFree(tempv);
	// mxDestroyArray(pstate_array);
	// mxDestroyArray(pjnbatch);

	return total_states;
}


/////////////////////////////////////////////////////DESIGN CONTROLLERS///////////////////////////////////////////////////////////

/***********************************
 
	This function solves the safety problem (deterministic)...
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array 
		BDD_TS - index for the Transition BDD for safety
		BDD_WS - index for the W BDD for safety
		BDD_FPS - index for the result of the fixed-point algorithm for safety
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: none
	
************************************/
void FPBDD_safety(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_TS, int BDD_WS, int BDD_FPS, s_vector* params_symb)
{
	int FP_count = 0;
	int totbits = params_symb->totbits;
	int nbitsloop = params_symb->nbitsloop;
	long i,j,t,s,l;
	long k=1;
	long r=0;
	double FP;

	// Temporary nodes
	// Zp (Z prime) holds the final answer
	DdNode *Z, *Zp;
	DdNode* cube_array, *tmpz;
	DdNode* tmp, *tmp2, *tmpperm;

	// Count iterations through fixed point algorithm
	double *nFP;
	mxArray *pnFP;
	pnFP = mxCreateDoubleScalar(0);
	nFP=(double *) mxGetData(pnFP);

	mexPrintf("\n                           Safety Controller \n");

////////////////////////////////////////////////////////////////////////////////////////

	//Build for permutation of x,u,x' to x',u,x
	int *permutation =  new int[totbits];
	
	//totbits=2*bits(x)+bits(u), nbitsloop=bits(x)+bits(u)

	//initialize u
	for (t=totbits-nbitsloop; t<nbitsloop; t++)
		permutation[t] = t;

	//reorder initial x (u remains the same)
	for (j=0; j<totbits-nbitsloop; j++)
		permutation[j] = nbitsloop+j;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	int *existential = new int[totbits];

	//initialize everything to 2
	for (s=0; s<totbits-nbitsloop; s++)
		existential[s] = 2;

	//final x and u set to 1 -- rest is ignored
	for (l=totbits-nbitsloop; l<totbits; l++)
		existential[l] = 1;	

	cube_array = Cudd_CubeArrayToBdd(ddmanl, existential);
	Cudd_Ref(cube_array);

////////////////////////////////////////////////////////////////////////////////////////

	// Here we initialize Z to 1, then copy W (target) onto Zp
 	Z=Cudd_ReadOne(ddmanl);
	Cudd_Ref(Z);
	Zp=ddnodearrayl[BDD_WS];
	Cudd_Ref(Zp);

	while (Cudd_bddLeq(ddmanl,Z,Zp) != 1) // while Z>Zp
	{	
	
		// Pop-up messages for # of iterations through the fixed point algorithm
		if(FP_count >= k*25)
		{
			if(r>0)
				mexEvalString("delete(g)");
		
			nFP[0] = FP_count;
			mexPutVariable("caller", "nFP", pnFP);
			mexEvalString("nFP=num2str(nFP);");
			mexEvalString("FP_string=strcat('Fixed-point iteration: ',nFP);");
			mexEvalString("g=helpdlg(FP_string,'Pessoa');");
			k++;
			r++;
		}

		// Counting the number of iterations
		FP_count++;

		// Eliminate an iteration if Z and Zp are already equal
		if(Cudd_bddLeq(ddmanl, Z, Zp) && Cudd_bddLeq(ddmanl, Zp, Z))
			FP_count--;
	
		// Decrease reference count of Z, set equal to Zp
		Cudd_RecursiveDeref(ddmanl, Z);
		Z = Zp;
		
		// Zp=exists u,x': T AND Z AND Z(x')
		// Permute Z (x,u,x' to x',u,x)
		tmpz = Cudd_bddAnd(ddmanl,Z,Cudd_bddPermute(ddmanl, Z, permutation));
		Cudd_Ref(tmpz);

		// Existential abstraction of the AND of permuted Z and transitions BDD
		Zp = Cudd_bddAndAbstract(ddmanl, ddnodearrayl[BDD_TS],tmpz, cube_array);
		Cudd_Ref(Zp);
		Cudd_RecursiveDeref(ddmanl, tmpz);
	}
		
	if(FP_count > 25)		
		mexEvalString("delete(g);");

	mexPrintf("\nFixed point algorithm finished after ");
	mexPrintf("%i", FP_count);
	mexPrintf(" iteration(s). \n");
	
	Cudd_RecursiveDeref(ddmanl, cube_array);
	Cudd_RecursiveDeref(ddmanl, Zp);
	Cudd_RecursiveDeref(ddmanl,ddnodearrayl[BDD_WS]);

	// Final controller domain is copied from Z to BDD_WS
	ddnodearrayl[BDD_WS]=Z;
	Cudd_Ref(ddnodearrayl[BDD_WS]);

	// Permute Z (x,u,x' to x',u,x)
	tmpperm=Cudd_bddPermute(ddmanl, Z, permutation);
	Cudd_Ref(tmpperm);

	// Result TFP= T AND W AND W(x')
	// Take the AND of original Z (target) and permuted Z
	tmp=Cudd_bddAnd(ddmanl, Z, tmpperm);
	Cudd_Ref(tmp);
	Cudd_RecursiveDeref(ddmanl, tmpperm);
	Cudd_RecursiveDeref(ddmanl, Z);
	
	// AND the result with transitions BDD
	tmp2=Cudd_bddAnd(ddmanl,tmp, ddnodearrayl[BDD_TS]);
	Cudd_Ref(tmp2);
	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_FPS]);
	// Final controller answer is copied from tmp2 to BDD_FPS
	ddnodearrayl[BDD_FPS] = tmp2;
	Cudd_RecursiveDeref(ddmanl, tmp);
}

/***********************************
 
	This function solves the reachability problem...
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array 
		BDD_TR - index for the Transition BDD for reachability
		BDD_WR - index for the W BDD for reachability
		BDD_FPR - index for the result of the fixed-point algorithm for reachability
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: none
	
************************************/
void FPBDD_reachability(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_TR, int BDD_WR, int BDD_FPR, s_vector* params_symb)
{
	int FP_count = 0;
	int totbits = params_symb->totbits;
	int nbitsloop = params_symb->nbitsloop;
	int n = params_symb->n;
	int m = params_symb->m;
	long i,j,t,l,s;
	long r=0;
	long k=1;
	double FP;

	DdNode *Z, *Zp;
	DdNode* cube_array, *tmpperm;
 	DdNode *tmpTc, *tmpTcp, *tmpZ1, *tmpZ2, *tmpZ3;

	//Count iterations through fixed point algorithm
	double *nFP;
	mxArray *pnFP;
	pnFP = mxCreateDoubleScalar(0);
	nFP=(double *) mxGetData(pnFP);

	mexPrintf("\n                        Reachability Controller \n");

////////////////////////////////////////////////////////////////////////////////////////

	//Build for permutation of x,u,x' to x',u,x
	int *permutation =  new int[totbits];
	
	//totbits=2*bits(x)+bits(u), nbitsloop=bits(x)+bits(u)

	//initialize u
	for (t=totbits-nbitsloop; t<nbitsloop; t++)
		permutation[t] = t;

	//reorder initial x (u remains the same)
	for (j=0; j<totbits-nbitsloop; j++)
		permutation[j] = nbitsloop+j;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	int *existential = new int[totbits];

	//initialize everything to 2
	for (s=0; s<totbits-nbitsloop; s++)
		existential[s] = 2;

	//final x and u set to 1 -- rest is ignored
	for (l=totbits-nbitsloop; l<totbits; l++)
		existential[l] = 1;	

	cube_array = Cudd_CubeArrayToBdd(ddmanl, existential);
	Cudd_Ref(cube_array);

////////////////////////////////////////////////////////////////////////////////////////

	// Here we initialize Z to 0, then copy W(target) onto Zp
	Z=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(Z);
	Zp=ddnodearrayl[BDD_WR];
	Cudd_Ref(Zp);

	// Initialize Tc to one
	tmpTc=Cudd_ReadOne(ddmanl);
	Cudd_Ref(tmpTc);

	// Initiazlie Tc' to zero
	tmpTcp=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(tmpTcp);

	do
	{

		// Pop-up messages for # of iterations through the fixed point algorithm
		if(FP_count >= k*25)
		{
			if(r>0)
				mexEvalString("delete(g)");
		
			nFP[0] = FP_count;
			mexPutVariable("caller", "nFP", pnFP);
			mexEvalString("nFP=num2str(nFP);");
			mexEvalString("FP_string=strcat('Fixed-point iteration: ',nFP);");
			mexEvalString("g=helpdlg(FP_string,'Pessoa');");

			k++;
			r++;
		}

		// Counting the number of iterations
		FP_count++;	

		// Eliminate an iteration if Z and Zp are already equal
		if(Cudd_bddLeq(ddmanl, tmpTcp, tmpTc) && Cudd_bddLeq(ddmanl, tmpTc, tmpTcp))
			FP_count--;

		// Decrease reference count of Z, set equal to Zp
		Cudd_RecursiveDeref(ddmanl, Z);

		// Take the OR of Z' and original W (target)
		Z = Cudd_bddOr(ddmanl, Zp, ddnodearrayl[BDD_WR]);
		Cudd_Ref(Z);

		Cudd_RecursiveDeref(ddmanl, tmpTc);
		tmpTc = tmpTcp;

		// Permute Z (x,u,x' to x',u,x)
		tmpperm=Cudd_bddPermute(ddmanl, Z, permutation);
		Cudd_Ref(tmpperm);

		///////////////////////////////////////////////////////////////////////////////
		//CONTROLLER
		///////////////////////////////////////////////////////////////////////////////

		// Zp=exists u,x': Tc'(result) AND T AND (NOT Z AND Z_permuted)
		// (NOT is done by NAND-ing BDD and 1)
		tmpZ2=Cudd_bddAnd(ddmanl, Cudd_bddNand(ddmanl, Z, Cudd_ReadOne(ddmanl)), tmpperm);
		Cudd_Ref(tmpZ2);
		Cudd_RecursiveDeref(ddmanl, tmpperm);

		tmpZ3=Cudd_bddAnd(ddmanl, ddnodearrayl[BDD_TR], tmpZ2);
		Cudd_Ref(tmpZ3);
		Cudd_RecursiveDeref(ddmanl, tmpZ2);

		tmpTcp=Cudd_bddOr(ddmanl, tmpTc, tmpZ3);
		Cudd_Ref(tmpTcp);
		Cudd_RecursiveDeref(ddmanl, tmpZ3);

		// Existential abstraction 
		Zp = Cudd_bddExistAbstract(ddmanl, tmpTcp, cube_array);
		Cudd_Ref(Zp);

	}while (Cudd_bddLeq(ddmanl, tmpTcp, tmpTc) != 1);

	if(FP_count > 25)
		mexEvalString("delete(g);");

	mexPrintf("\nFixed point algorithm finished after ");
	mexPrintf("%i", FP_count);
	mexPrintf(" iteration(s). \n");

	Cudd_RecursiveDeref(ddmanl, cube_array);

	// Final controller answer is copied from Tc' to BDD_FPR
	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_FPR]);
	ddnodearrayl[BDD_FPR] = tmpTcp;
	Cudd_Ref(ddnodearrayl[BDD_FPR]);
	Cudd_RecursiveDeref(ddmanl, tmpTc);
	Cudd_RecursiveDeref(ddmanl, tmpTcp);

	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_WR]);
	// Final controller domain is copied from Z' to BDD_WR
	ddnodearrayl[BDD_WR] = Zp;
	Cudd_Ref(ddnodearrayl[BDD_WR]);
	Cudd_RecursiveDeref(ddmanl, Z);
	Cudd_RecursiveDeref(ddmanl, Zp);

}

/***********************************
 
	This function solves the safety problem (non-deterministic)...
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array 
		BDD_TS - index for the Transition BDD for safety
		BDD_WS - index for the W BDD for safety
		BDD_FPS - index for the result of the fixed-point algorithm for safety
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: none

************************************/
void FPBDD_safety_nondeter(DdManager* ddmanl,DdNode** ddnodearrayl,int BDD_TS, int BDD_WS, int BDD_FPS, s_vector* params_symb)
{
	int FP_count = 0;
	int totbits = params_symb->totbits;
	int nbitsloop = params_symb->nbitsloop;
	int nbitsx = params_symb->nbitsx;
	long i,j,t,s,l,s2,l2;
	long r=0;
	long k=1;
	double FP;

	DdNode *Z, *Zp;
	DdNode *cube_array, *cube_array2;
 	DdNode *tmpTc, *tmpTcp;
	DdNode *tmpperm, *tmpZ1, *tmpZ2, *tmpZ3, *tmpZ4;

	//Count iterations through fixed point algorithm
	double *nFP;
	mxArray *pnFP;
	pnFP = mxCreateDoubleScalar(0);
	nFP=(double *) mxGetData(pnFP);

	mexPrintf("\n                           Safety Controller \n");

////////////////////////////////////////////////////////////////////////////////////////

	//Build for permutation of x,u,x' to x',u,x
	int *permutation =  new int[totbits];
	
	// nbitsloop=bits(x)+bits(u)

	//initialize u
	for (t=totbits-nbitsloop; t<nbitsloop; t++)
		permutation[t] = t;

	//reorder initial x (u remains the same)
	for (j=0; j<totbits-nbitsloop; j++)
		permutation[j] = nbitsloop+j;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	int *existential = new int[totbits];

	//initialize everything to 2
	for (s=0; s<totbits-nbitsloop; s++)
		existential[s] = 2;

	//final x and u set to 1 -- rest is ignored
	for (l=totbits-nbitsloop; l<totbits; l++)
		existential[l] = 1;	

	// another existential for x' only
	int *existential2 = new int[totbits];

	//initialize everything to 2 
	for (s2=0; s2<totbits-nbitsx; s2++)
		existential2[s2] = 2;

	//final x to 1 -- rest is ignored
	for (l2=totbits-nbitsx; l2<totbits; l2++)
		existential2[l2] = 1;

	cube_array = Cudd_CubeArrayToBdd(ddmanl, existential);
	Cudd_Ref(cube_array);

	cube_array2 = Cudd_CubeArrayToBdd(ddmanl, existential2);
	Cudd_Ref(cube_array2);

////////////////////////////////////////////////////////////////////////////////////////

	// Here we initialize Z to 1, then copy W onto Zp
 	Z=Cudd_ReadOne(ddmanl);
	Cudd_Ref(Z);
	Zp=ddnodearrayl[BDD_WS];
	Cudd_Ref(Zp);

	// Initialize Tc to zero
	tmpTc=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(tmpTc);

	// Initiazlie Tc' to transitions BDD
	tmpTcp=ddnodearrayl[BDD_TS];
	Cudd_Ref(tmpTcp);

	do
	{	

		// Pop-up messages for # of iterations through the fixed point algorithm
		if(FP_count >= k*25)
		{
			if(r>0)
				mexEvalString("delete(g)");
		
			nFP[0] = FP_count;
			mexPutVariable("caller", "nFP", pnFP);
			mexEvalString("nFP=num2str(nFP);");
			mexEvalString("FP_string=strcat('Fixed-point iteration: ',nFP);");
			mexEvalString("g=helpdlg(FP_string,'Pessoa');");
			k++;
			r++;
		}

		// Counting the number of iterations
		FP_count++;

		// Eliminate an iteration if Z and Zp are already equal
		if(Cudd_bddLeq(ddmanl, tmpTcp, tmpTc) && Cudd_bddLeq(ddmanl, tmpTc, tmpTcp))
			FP_count--;

		// Decrease reference count of Z, set equal to Zp
		Cudd_RecursiveDeref(ddmanl, Z);
		Z = Zp;

		// Decrease reference count of Tc, set equal to Tc'
		Cudd_RecursiveDeref(ddmanl, tmpTc);
		tmpTc = tmpTcp;

		// Permute Z (x,u,x' to x',u,x)
		tmpperm=Cudd_bddPermute(ddmanl, Z, permutation);
		Cudd_Ref(tmpperm);

		// Z1 is NOT tmpperm (NOT is done by NAND-ing BDD and 1)
		tmpZ1=Cudd_bddNand(ddmanl, tmpperm, Cudd_ReadOne(ddmanl));
		Cudd_Ref(tmpZ1);
		Cudd_RecursiveDeref(ddmanl, tmpperm);

		//tmpZ1=NOT(Z(x'))
		// Existential abstraction x' of the AND of Tc' and tmpZ1
		tmpZ2 = Cudd_bddAndAbstract(ddmanl, tmpTcp,tmpZ1, cube_array2);
		Cudd_Ref(tmpZ2);
		Cudd_RecursiveDeref(ddmanl, tmpZ1);
		
		//tmpZ3=NOT(exists x': Tc AND NOT(Z(x')))
		tmpZ3=Cudd_bddNand(ddmanl, tmpZ2, Cudd_ReadOne(ddmanl));
		Cudd_Ref(tmpZ3);
		Cudd_RecursiveDeref(ddmanl, tmpZ2);

		//tmpZ4=Z AND NOT(exists x': Tc AND NOT(Z(x')))
		tmpZ4=Cudd_bddAnd(ddmanl, tmpZ3, Z);
		Cudd_Ref(tmpZ4);
		Cudd_RecursiveDeref(ddmanl, tmpZ3);

		// tmpTcp= Tc AND NOT(tmpZ2) AND Z(x)
		tmpTcp=Cudd_bddAnd(ddmanl, tmpTc, tmpZ4);
		Cudd_Ref(tmpTcp);
		Cudd_RecursiveDeref(ddmanl, tmpZ4);

		// Zp=exists u,x': Tc AND NOT(tmpZ2) AND Z(x)
		Zp = Cudd_bddExistAbstract(ddmanl, tmpTcp, cube_array);
		Cudd_Ref(Zp);

	}while (Cudd_bddLeq(ddmanl,tmpTc,tmpTcp) != 1); // Tc>Tp

	if(FP_count > 25)
		mexEvalString("delete(g);");		
	
	mexPrintf("\nFixed point algorithm finished after ");
	mexPrintf("%i", FP_count);
	mexPrintf(" iteration(s). \n");

	Cudd_RecursiveDeref(ddmanl, cube_array);
	Cudd_RecursiveDeref(ddmanl, cube_array2);

	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_FPS]);
	// Final controller answer is copied from Tc' to BDD_FPS
	ddnodearrayl[BDD_FPS] = tmpTcp;
	Cudd_Ref(ddnodearrayl[BDD_FPS]);
	Cudd_RecursiveDeref(ddmanl, tmpTc);
	Cudd_RecursiveDeref(ddmanl, tmpTcp);

	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_WS]);
	// Final controller domain is copied from Z' to BDD_WS
	ddnodearrayl[BDD_WS] = Zp;
	Cudd_Ref(ddnodearrayl[BDD_WS]);
	Cudd_RecursiveDeref(ddmanl, Z);
	Cudd_RecursiveDeref(ddmanl, Zp);
}


/***********************************
 
	This function solves the reachability problem (non-deterministic)...
	
	INPUTS: ddmanl - CUDD manager 
		ddnodearrayl - CUDD node array 
		BDD_TR - index for the Transition BDD for reachability
		BDD_WR - index for the W BDD for reachability
		BDD_FPR - index for the result of the fixed-point algorithm for reachability
 		params_symb   - data structure containing discretization
                      information
	OUTPUTS: none
	
************************************/
void FPBDD_reachability_nondeter(DdManager* ddmanl,DdNode** ddnodearrayl, int BDD_TR, int BDD_WR, int BDD_FPR, s_vector* params_symb)
{
	int FP_count = 0;
	int totbits = params_symb->totbits;
	int nbitsloop = params_symb->nbitsloop;
	int n = params_symb->n;
	int m = params_symb->m;
	int nbitsx = params_symb->nbitsx;
	long i,j,t,s,l,s2,l2;
	long r=0;
	long k=1;
	double FP;

	// Temporary BDD's
	DdNode *Z, *Zp;
	DdNode *cube_array, *cube_array2;
 	DdNode *tmpTc, *tmpTcp;
	DdNode *tmpperm, *tmpZ1, *tmpZ2, *tmpZ3, *tmpZ4, *tmpZ5, *tmpZ6;

	//Count iterations through fixed point algorithm
	double *nFP;
	mxArray *pnFP;
	pnFP = mxCreateDoubleScalar(0);
	nFP=(double *) mxGetData(pnFP);

	mexPrintf("\n                        Reachability Controller \n");

////////////////////////////////////////////////////////////////////////////////////////

	//Build for permutation of x,u,x' to x',u,x
	int *permutation =  new int[totbits];
	
	//initialize u
	for (t=totbits-nbitsloop; t<nbitsloop; t++)
		permutation[t] = t;

	//reorder initial x (u remains the same)
	for (j=0; j<totbits-nbitsloop; j++)
		permutation[j] = nbitsloop+j;

	//reorder final x
	for (i=nbitsloop; i<totbits; i++)
		permutation[i] = i-nbitsloop;	

////////////////////////////////////////////////////////////////////////////////////////

	// Build cube for existential check of u and x'
	int *existential = new int[totbits];

	//initialize everything to 2
	for (s=0; s<totbits-nbitsloop; s++)
		existential[s] = 2;

	//final x and u set to 1 -- rest is ignored
	for (l=totbits-nbitsloop; l<totbits; l++)
		existential[l] = 1;	

	// another existential for x' only
	int *existential2 = new int[totbits];

	//initialize everything to 2
	for (s2=0; s2<totbits-nbitsx; s2++)
		existential2[s2] = 2;

	//final x to 1 -- rest is ignored
	for (l2=totbits-nbitsx; l2<totbits; l2++)
		existential2[l2] = 1;

	cube_array = Cudd_CubeArrayToBdd(ddmanl, existential);
	Cudd_Ref(cube_array);

	cube_array2 = Cudd_CubeArrayToBdd(ddmanl, existential2);
	Cudd_Ref(cube_array2);

////////////////////////////////////////////////////////////////////////////////////////

	// Here we initialize Z to 0, then copy W (target) onto Zp 
 	Zp=ddnodearrayl[BDD_WR];
	Cudd_Ref(Zp);
	Z=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(Z);

	// Initialize Tc to zero
	tmpTc=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(tmpTc);

	// Initialize Tc' to zero
	tmpTcp=Cudd_ReadLogicZero(ddmanl);
	Cudd_Ref(tmpTcp);

	do
	{	
		// Pop-up messages for # of iterations through the fixed point algorithm
		if(FP_count >= k*25)
		{
			if(r>0)
				mexEvalString("delete(g)");
		
			nFP[0] = FP_count;
			mexPutVariable("caller", "nFP", pnFP);
			mexEvalString("nFP=num2str(nFP);");
			mexEvalString("FP_string=strcat('Fixed-point iteration: ',nFP);");
			mexEvalString("g=helpdlg(FP_string,'Pessoa');");

			k++;
			r++;
		}

		// Counting the number of iterations
		FP_count++;

		// Eliminate an iteration if Z and Zp are already equal
		if(Cudd_bddLeq(ddmanl, tmpTcp, tmpTc) && Cudd_bddLeq(ddmanl, tmpTc, tmpTcp))
			FP_count--;

		// Decrease reference count of Z
		Cudd_RecursiveDeref(ddmanl, Z);
		// Z = Z' OR BDD_WR (target)
		Z = Cudd_bddOr(ddmanl, Zp, ddnodearrayl[BDD_WR]);
		Cudd_Ref(Z);

		// Decrease reference count of Tc, set equal to Tc'
		Cudd_RecursiveDeref(ddmanl, tmpTc);
		tmpTc = tmpTcp;

		// Permute Z (x,u,x' to x',u,x)
		tmpperm=Cudd_bddPermute(ddmanl, Z, permutation);
		Cudd_Ref(tmpperm);

		// Z1 is NOT tmpperm (NOT is done by NAND-ing BDD and 1)
		tmpZ1=Cudd_bddNand(ddmanl, tmpperm, Cudd_ReadOne(ddmanl));
		Cudd_Ref(tmpZ1);
		Cudd_RecursiveDeref(ddmanl, tmpperm);

		//tmpZ1=NOT(Z(x'))
		// Existential abstraction x' of the AND of T and tmpZ1
		// tmpZ2= exists x': T AND NOT(Z(x'))
		tmpZ2 = Cudd_bddAndAbstract(ddmanl, ddnodearrayl[BDD_TR],tmpZ1, cube_array2);
		Cudd_Ref(tmpZ2);
		Cudd_RecursiveDeref(ddmanl, tmpZ1);

		//tmpZ3=NOT(exists x': Tc AND NOT(Z(x')))	
		tmpZ3=Cudd_bddNand(ddmanl, tmpZ2, Cudd_ReadOne(ddmanl));
		Cudd_Ref(tmpZ3);
		Cudd_RecursiveDeref(ddmanl, tmpZ2);

		//tmpZ4=Z AND NOT(exists x': Tc AND NOT(Z(x')))
		//tmpZ4 = T AND NOT(tmpZ2)
		tmpZ4=Cudd_bddAnd(ddmanl, tmpZ3, ddnodearrayl[BDD_TR]);
		Cudd_Ref(tmpZ4);
		Cudd_RecursiveDeref(ddmanl, tmpZ3);
		
		//tmpZ5 = NOT(Z) 
		tmpZ5=Cudd_bddNand(ddmanl, Z, Cudd_ReadOne(ddmanl));
		Cudd_Ref(tmpZ5);

		//tmpZ6= tmpZ4 AND NOT(Z)		
		tmpZ6=Cudd_bddAnd(ddmanl, tmpZ4, tmpZ5);
		Cudd_Ref(tmpZ6);
		Cudd_RecursiveDeref(ddmanl, tmpZ4);
		Cudd_RecursiveDeref(ddmanl, tmpZ5);
		
		//Tp=Tc OR tmpZ6= Tc OR (T AND NOT(exists x':T AND NOT(Z(x'))) AND NOT(Z))
		tmpTcp=Cudd_bddOr(ddmanl, tmpTc, tmpZ6);
		Cudd_Ref(tmpTcp);
		Cudd_RecursiveDeref(ddmanl, tmpZ6);

		// Zp=exists u,x': Tc AND NOT(tmpZ2) AND Z(x)
		Zp = Cudd_bddExistAbstract(ddmanl, tmpTcp, cube_array);
		Cudd_Ref(Zp);

	}while (Cudd_bddLeq(ddmanl, tmpTcp, tmpTc) != 1); // Tp>Tc

	if(FP_count > 25)
		mexEvalString("delete(g);");		
	
	mexPrintf("\nFixed point algorithm finished after ");
	mexPrintf("%i", FP_count);
	mexPrintf(" iteration(s). \n");

	Cudd_RecursiveDeref(ddmanl, cube_array);
	Cudd_RecursiveDeref(ddmanl, cube_array2);

	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_FPR]);
	// Final controller answer is copied from Tc' to BDD_FPR
	ddnodearrayl[BDD_FPR] = tmpTcp;
	Cudd_Ref(ddnodearrayl[BDD_FPR]);
	Cudd_RecursiveDeref(ddmanl, tmpTc);
	Cudd_RecursiveDeref(ddmanl, tmpTcp);

	Cudd_RecursiveDeref(ddmanl, ddnodearrayl[BDD_WR]);
	// Final controller domain is copied from Z' to BDD_WR
	ddnodearrayl[BDD_WR] = Zp;
	Cudd_Ref(ddnodearrayl[BDD_WR]);
	Cudd_RecursiveDeref(ddmanl, Z);
	Cudd_RecursiveDeref(ddmanl, Zp);

}

