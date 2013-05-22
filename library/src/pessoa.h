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

#ifndef _PESSOA_H_
#define _PESSOA_H_

#include "mex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "../../cudd-2.4.2/cudd/cudd.h"
#include "../../cudd-2.4.2/dddmp/dddmp.h"

typedef struct t_sv{
	int n;
	int m;
	double* nume;
	double* nbits; 
	int totbits;
	int nbitsloop;
	int deter;
	int nbitsx;
	double* nxbits;
}s_vector;

int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);
bool IsInBDD(DdManager* ddman,DdNode** ddnodearray, int bddnum, double *array, s_vector* params_symb, int nnm);
int PlotFSM(DdManager* ddman,DdNode** ddnodearray, int fsm, s_vector* params_symb);
int AddToBddBatch(DdManager* ddman,DdNode** ddnodearray, int bddnum, double *array, s_vector* params_symb, long nbatch, int nnm);
int buildBDDT(DdManager* ddman,DdNode** ddnodearray, int BDD_T, s_vector* params_symb, long nbatch, double totloops);
int ValidUInBDD(DdManager* ddman,DdNode** ddnodearray,DdNode* controldd, DdNode* existbdd, double *arrayx, double *arrayuold, double *arrayu, s_vector* params_symb);
int PlotSet(DdManager* ddman,DdNode** ddnodearray, int fsm, s_vector* params_symb);
long buildBDDW(DdManager* ddman,DdNode** ddnodearray, int BDD_W, double *wmin, double *wmax, s_vector* params_symb, long nbatch);
long build_charfBDDW(DdManager* ddman,DdNode** ddnodearray, int BDD_W, double *wmin, double *wmax, s_vector* params_symb, long nbatch);
void FPBDD_safety(DdManager* ddman,DdNode** ddnodearray, int BDD_T, int BDD_W, int BDD_FPS, s_vector* params_symb);
void FPBDD_reachability(DdManager* ddman,DdNode** ddnodearray, int BDD_T, int BDD_W, int BDD_FPR, s_vector* params_symb);
void FPBDD_safety_nondeter(DdManager* ddman,DdNode** ddnodearray,int BDD_T, int BDD_W, int BDD_FPS, s_vector* params_symb);
void FPBDD_reachability_nondeter(DdManager* ddman,DdNode** ddnodearray, int BDD_T, int BDD_W, int BDD_FPR, s_vector* params_symb);

#endif
