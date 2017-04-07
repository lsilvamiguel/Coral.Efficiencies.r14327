#include "TAlgo.h"

/*!
  Inversion of a positive definite symmetric 5x5 matrix, 
  represented by it's lower triangle 
  \param a input matrix (array of 15)
  \return b - inverted matrix (array of 15)
*/


int TAlgo::Inv5(double* a, double* b)

{
    double x1, x2, x3, x4, x5, y1, y2, y3, y4, y5;

// ---------------------------------------------------------------------*
// Inversion of a positive definite symmetric matrix (5x5)              * 
// by a modification of the Gauss-Jordan method                         * 
//                                                                      * 
// Input parameters  : A(1/15) - the elements of the lower triangle of  * 
//                               the matrix to be inverted (double pr.) * 
//                               1                                      * 
//                               2   3                                  * 
//                               4   5  6                               * 
//                               7   8  9 10                            * 
//                               11 12 13 14 15                         * 
// Output parameters : B(1/15) - inverted matrix           (double pr.) * 
//                     return     = 0 O'K                               * 
//                     return     = 1 if A is not positive-definite     * 
//                                                                      * 
// Author : I.Gavrilenko                                                * 
// ---------------------------------------------------------------------*



    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */

    if (a[1] <= double(0)) {
	return 1;
    }
    x1 = double(1) / a[1];
    x2 = -a[2] * x1;
    x3 = -a[4] * x1;
    x4 = -a[7] * x1;
    x5 = -a[11] * x1;
    b[1] = a[3] + a[2] * x2;
    b[2] = a[5] + a[4] * x2;
    b[3] = a[6] + a[4] * x3;
    b[4] = a[8] + a[7] * x2;
    b[5] = a[9] + a[7] * x3;
    b[6] = a[10] + a[7] * x4;
    b[7] = a[12] + a[11] * x2;
    b[8] = a[13] + a[11] * x3;
    b[9] = a[14] + a[11] * x4;
    b[10] = a[15] + a[11] * x5;
    if (b[1] <= double(0)) {
	return 1;
    }
    y1 = double(1) / b[1];
    y2 = -b[2] * y1;
    y3 = -b[4] * y1;
    y4 = -b[7] * y1;
    y5 = x2 * y1;
    b[1] = b[3] + b[2] * y2;
    b[2] = b[5] + b[4] * y2;
    b[3] = b[6] + b[4] * y3;
    b[4] = b[8] + b[7] * y2;
    b[5] = b[9] + b[7] * y3;
    b[6] = b[10] + b[7] * y4;
    b[7] = x3 + x2 * y2;
    b[8] = x4 + x2 * y3;
    b[9] = x5 + x2 * y4;
    b[10] = x1 + x2 * y5;
    if (b[1] <= double(0)) {
	return 1;
    }
    x1 = double(1) / b[1];
    x2 = -b[2] * x1;
    x3 = -b[4] * x1;
    x4 = b[7] * x1;
    x5 = y2 * x1;
    b[1] = b[3] + b[2] * x2;
    b[2] = b[5] + b[4] * x2;
    b[3] = b[6] + b[4] * x3;
    b[4] = b[8] + b[7] * x2;
    b[5] = b[9] + b[7] * x3;
    b[6] = b[10] + b[7] * x4;
    b[7] = y3 + y2 * x2;
    b[8] = y4 + y2 * x3;
    b[9] = y5 + y2 * x4;
    b[10] = y1 + y2 * x5;
    if (b[1] <= double(0)) {
	return 1;
    }
    y1 = double(1) / b[1];
    y2 = -b[2] * y1;
    y3 = b[4] * y1;
    y4 = b[7] * y1;
    y5 = x2 * y1;
    b[1] = b[3] + b[2] * y2;
    b[2] = b[5] + b[4] * y2;
    b[3] = b[6] + b[4] * y3;
    b[4] = b[8] + b[7] * y2;
    b[5] = b[9] + b[7] * y3;
    b[6] = b[10] + b[7] * y4;
    b[7] = x3 + x2 * y2;
    b[8] = x4 + x2 * y3;
    b[9] = x5 + x2 * y4;
    b[10] = x1 + x2 * y5;
    if (b[1] <= double(0)) {
	return 1;
    }
    b[15] = double(1) / b[1];
    b[11] = b[2] * b[15];
    b[12] = b[4] * b[15];
    b[13] = b[7] * b[15];
    b[14] = y2 * b[15];
    b[1] = b[3] + b[2] * b[11];
    b[2] = b[5] + b[4] * b[11];
    b[3] = b[6] + b[4] * b[12];
    b[4] = b[8] + b[7] * b[11];
    b[5] = b[9] + b[7] * b[12];
    b[6] = b[10] + b[7] * b[13];
    b[7] = y3 + y2 * b[11];
    b[8] = y4 + y2 * b[12];
    b[9] = y5 + y2 * b[13];
    b[10] = y1 + y2 * b[14];
    return 0;
}






