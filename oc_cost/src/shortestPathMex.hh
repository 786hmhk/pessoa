/*
 * shortestPathMex.hh
 *
 *  Created on: Jun 15, 2013
 *      Author: tanasaki
 */

#ifndef SHORTESTPATHMEX_HH_
#define SHORTESTPATHMEX_HH_


#include "mex.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

// CUDD Library
#include "cuddObj.hh"
// Pessoa Header
//#include "pessoa.h"
//
#include "ShortestPath.hh"

// TODO: delete this when used together with the Pessoa project.
// Add "pessoa.h" instead.
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


#define FILE_EXISTS(file) 	if ((file)==NULL) \
								{ \
									mexEvalString("warndlg('Symbolic model file not found!','Pessoa Error')");\
									mexErrMsgTxt("\nERROR: Symbolic model file not found!.\n"); \
								} \
								else \
								{ \
									fclose((file));\
								}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/*
 *
 */
class shortestPathMex {
public:
	shortestPathMex();
	virtual ~shortestPathMex();
};






#endif /* SHORTESTPATHMEX_HH_ */
