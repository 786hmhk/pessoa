/*
 * pessoacontrollerc.hh
 *
 *  Created on: Jul 14, 2013
 *      Author: tanasaki
 */

#ifndef PESSOA_CONTROLLER_HH_
#define PESSOA_CONTROLLER_HH_


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
#include "ShortestPath.hh"


/*
 *
 */
class pessoa_controller_c {
public:
	pessoa_controller_c();
	virtual ~pessoa_controller_c();
};

#endif /* PESSOACONTROLLERC_HH_ */
