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
 * Contains main function used to test the ShortestPath Class.
 *
 * 		No details yet.
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <sys/time.h>   /* time profiling */

#include <iostream> /* Writing to file */
#include <fstream>  /* Writing to file */

// Cudd
#include "cuddObj.hh"
//#include "util.h"
//#include "cudd.h"
//#include "util.h"

// Local
#include "ShortestPath.hh"

#define TRANSITION(xx,uu,xx_) 	BDD_transition(mgr,x,u,x_,no_state_vars,no_input_vars,xx,uu,xx_)


#define FILE_EXISTS(file) 	if ((file)==NULL) \
								{ \
									printf("\nERROR: Symbolic model file not found!.\n"); \
								} \
								else \
								{ \
									fclose((file));\
								}


void example_FW();
void example_DSP();
void example_NDSP();
void test_actual();

void get_S_xux(Cudd *mgr, BDD *T, int no_states, int no_inputs);
BDD BDD_transition(Cudd *mgr, BDD *x, BDD *u, BDD *x_, int no_states, int no_inputs, int xi, int ui, int xi_);
void get_S_cost_x(Cudd *mgr, ADD *C, int no_states, int no_inputs, int *costs);
BDD getTargetSet(Cudd *mgr, int no_states, std::vector<int> target_set);

unsigned int getNoBits(unsigned int number);
bool writeSysInfo_TBFile(int no_states, int no_inputs);
long long get_usec(void);

#endif /* MAIN_H_ */
