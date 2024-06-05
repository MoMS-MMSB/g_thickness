//#ifndef _matrix_h
//#define _matrix_h

#include <gromacs/types/simple.h>
#include <gromacs/smalloc.h>

real **realMatrix(int d1, int d2, real defval);
void deleteRealMat(real **mat, int d1);
int **intMatrix(int d1, int d2, int defval);
void deleteIntMat(int **mat, int d1);

//#endif	[> _matrix_h <]
