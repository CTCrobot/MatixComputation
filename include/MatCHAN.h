#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define MAX_DIM 10 
#define EPS_MC 2.2204e-16

void mc_fnvMatrixPrint(double **dInMatrix, int nRow, int nColumn);
void mc_fnvVectorPrint(double *dInVector, int nNumber);
double mc_fnfMatrixDeterminant(double dInMatrix[MAX_DIM][MAX_DIM], int nRow, int nColumn);
void mc_fnvMatrixTransport(double **dInRotMatrix, double **dOutRotMatrix,int nSizeMax, int nRow, int nColuw, double dDet);
void mc_fnvMatrixInverse(double **dInRotMatrix, double **dOutRotMatrix, int n);
double mc_fndVectorNorm(double *dVector,int nNumber);
void mc_fnvMatrixMultp(double **dMatrix1,double **dMatrix2, double **dOutMatrix, int nRow1, int nColumn1, int nColumn2);
void mc_fnvEularToRotMatrix(double dEl[3], double dOutRotMatrix[3][3]);
void mc_fnvRotMatrixToOmega(double dMatrix[3][3], double dOutOmega[3]);
void mc_fnvMatrixVectorMultp(double **dMatrix, double *dVector, double *dOutVector, int nRow, int nColumn);
void mc_fnvDiagMatrix(double *dInVector, double **dOutMatrix, int nNumber);
void mc_fnvIdentityMatrix( double **dOutMatrix, int nNumber, double dMultNumber, double *dPlusNumber); 
double* mc_fndVectorCopy(double *dInVector, int nNumber);
double** mc_fndMatixCopy(double **dInMatrix, int nRow, int nColumn);