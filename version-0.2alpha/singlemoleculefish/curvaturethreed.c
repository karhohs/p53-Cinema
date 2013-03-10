/*=================================================================
 *
 * 	CURVATURETHREED.C
 *	        Solves the curvature at each point in a 3 dimensional matrix
 *
 * The calling syntax is:
 *
 *		[C] = yprime(I)
 *
 *
 *
 *
 *
 *
 *=================================================================*/
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments and conditions*/
    if (nrhs != 1) {
        mexErrMsgTxt("One input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if ( mxIsChar(prhs[0]) || mxIsClass(prhs[0], "sparse") ||
            mxIsComplex(prhs[0]) ){
        mexErrMsgTxt("the input must be real, full, and nonstring");
    }
    
}