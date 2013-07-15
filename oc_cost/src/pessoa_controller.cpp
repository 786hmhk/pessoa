/*
 * pessoacontrollerc.cpp
 *
 *  Created on: Jul 14, 2013
 *      Author: tanasaki
 */

#include "pessoa_controller.hh"

pessoa_controller_c::pessoa_controller_c() {
	// TODO Auto-generated constructor stub

}

pessoa_controller_c::~pessoa_controller_c() {
	// TODO Auto-generated destructor stub
}

//! Mex Function.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	int verbose;

	/* Initialize the Manager. */
	Cudd mgr(0, 0);


	mexPrintf("\n------------------ pessoa_controller: TEST BEGIN ----------------- \n\n");

	/* Check for proper number of arguments. */
	if(nrhs<1 || nrhs>4) {
		mexErrMsgTxt("Four inputs required: filename, min and max for the System's state space. Optional fourth input: verbose flag.");
		return;
	}
	if(nrhs==4)
		verbose = (int)mxGetScalar(prhs[3]);
	else
		verbose = 0;





	mexPrintf("\n\n------------------ pessoa_controller: TEST END   ----------------- \n");
}
