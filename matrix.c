#include "matrix.h"

/** Create a matrix of real numbers
 *
 * The matrix is pre-filled with a value.
 *
 * Parameters:
 *  - d1: number of rows
 *  - d2: number of columns
 *  - defval: the value to put in every cells
 *
 * Return:
 *  The filled matrix.
 */
real **realMatrix(int d1, int d2, real defval) {
    int i, j;
    real **mat;
  
    smalloc(mat, d1 * sizeof(real*));
    for(i = 0; i < d1; i++){
        smalloc(mat[i], d2 * sizeof(real));
        for(j = 0; j < d2; j++) {
            mat[i][j] = defval;
        }
    }
    return mat;
}

/** Destroy a matrix of real numbers
 *
 * Parameters:
 *  - mat: the matrix to destroy
 *  - d1: the number of rows in the matrix
 */
void deleteRealMat(real **mat, int d1) {
    int i;
    for(i = 0; i < d1; i++) {
        sfree(mat[i]);
    }
    sfree(mat);
}

/** Create a matrix of integer
 *
 * The matrix is pre-filled with a value.
 *
 * Parameters:
 *  - d1: number of rows
 *  - d2: number of columns
 *  - defval: the value to put in every cells
 *
 * Return:
 *  The filled matrix.
 */
int **intMatrix(int d1, int d2, int defval) {
    int i, j;
    int **mat;
  
    smalloc(mat, d1 * sizeof(int*));
    for(i = 0; i < d1; i++){
        smalloc(mat[i], d2 * sizeof(real));
        for(j = 0; j < d2; j++) {
            mat[i][j] = defval;
        }
    }
    return mat;
}

/** Destroy a matrix of integer
 *
 * Parameters:
 *  - mat: the matrix to destroy
 *  - d1: the number of rows in the matrix
 */
void deleteIntMat(int **mat, int d1) {
    int i;
    for(i = 0; i < d1; i++) {
        sfree(mat[i]);
    }
    sfree(mat);
}


