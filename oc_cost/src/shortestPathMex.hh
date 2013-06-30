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
 * Contains main mex function to implement sp.
 *
 * 		No details yet.
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