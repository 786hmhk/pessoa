/*
 * main.h
 *
 *  Created on: Apr 13, 2013
 *      Author: tanasaki
 */

#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include <pthread.h>

// Cudd
#include "cuddObj.hh"
//#include "util.h"
//#include "cudd.h"
//#include "util.h"

// Local
#include "ShortestPath.h"


void example_FW();
void example_DSP();

void get_S_xux(Cudd *mgr, BDD *T);
BDD BDD_transition(Cudd *mgr, BDD *x, BDD *u, BDD *x_, int no_states, int no_inputs, int xi, int ui, int xi_);
void get_S_cost_x(Cudd *mgr, ADD *T);

#endif /* MAIN_H_ */
