/**
 * @file
 * @author Athanasios Tasoglou <A.Tasoglou@student.tudelft.nl>
 * @version 0.92
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
 * Contains...
 *
 * 		No details yet.
 */

#ifndef PESSOACOST2ADD_HH_
#define PESSOACOST2ADD_HH_





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
#include "dddmp.h"   // for storing DD into files.
// Pessoa Header
//#include "pessoa.h"
//
#include "lib/ShortestPath.hh" //TODO: IMPORTAT! ONLY NEEDED FOR filterCosts(). OTHERWISE REMOVE AND DELETE filterCosts().

//#define WAITBAR_ENABLE

#define KEEP_VALID_STATES 		1
#define KEEP_VALID_TRANSITIONS 	2


#define FILE_EXISTS(file) 	if ((file)==NULL) \
								{ \
									mexEvalString("warndlg('Symbolic model file not found!','Pessoa Error')");\
									mexErrMsgTxt("\nERROR: Symbolic model file not found!.\n"); \
								} \
								else \
								{ \
									fclose((file));\
								}

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

/*
 *
 */
class pessoa_cost2add {
private:

	int verbose;
	// Minimum/Maximum value of the state space. Needed to under-approximate
	// the state space if needed.
	double *wmin, *wmax;
	// This structure contains the parameters used in the construction of the symbolic abstraction.
	s_vector params_symb;
	// holds the MATLAB structure variable params_symb.
	mxArray *psv;
	// holds the MATLAB variable nbatch.
	long nbatch;
	// Number of states and inputs.
	unsigned int nstates, ninputs;
	// The name of the .add holding the States costs.
	char *SysStateCostADD_name;
	// The name of the .bdd holding the System representation.
	char *SysBDD_name; // delete if filtercosts() is removed.

	unsigned int total_add_states;

	int no_figure;

#ifdef WAITBAR_ENABLE
	double totloops;
#endif

	// CUDD Manager.
	Cudd *mgr;
	// ADD holding the System's State Costs
	ADD sys_state_cost;


	bool addBatchToADD(double *array, double *state_cost_array, s_vector* params_symb, long nbatch, int nnm);
	bool isInADD(double *array, double *cost);


public:
	pessoa_cost2add(mxArray **plhs, const mxArray **prhs, Cudd *mgr, int verbose);
	virtual ~pessoa_cost2add();

	void createSysStatesCost();

	ADD getSysStateCost();
	void setSysStateCost(ADD sysCosts);

	bool dumpSysStateCost();
	void dumpSysStateCostDot();

	void plotSysStateCost();

	//! Experimental. Filters out states depending on the mode.
	ADD filterCosts(int mode = KEEP_VALID_STATES);

};









#endif /* PESSOACOST2ADD_HH_ */
