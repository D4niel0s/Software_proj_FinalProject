#ifndef SYMNMF_H
#define SYMNMF_H

/*INCLUDES*/
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/*DEFINITIONS*/

/*A struct representing a point*/
typedef struct P{
    double *coords;
    int dim;
}Point;

/*PROTOTYPES*/

double **computeSimMat(Point*,int);
double **computeDegMat(double**,int);
double **computeNormSimMat(double**,double**,int);
double **Hoptimization(double**,double**,int,int,int,double);
double **UpdateH(double**,double**,int,int,double);

double squaredFrobeniusNorm(double**,int,int);
double **mulMat(double**,double**,int,int,int);
double **transpose(double**,int,int);
double **matDiff(double**,double**,int,int);
double eucDist(Point,Point);

int main(int, char **);
int countCommas(char *);
char *skipUntilComma(char []);

#endif /*SYMNMF_H*/