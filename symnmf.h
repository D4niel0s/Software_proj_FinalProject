/*INCLUDES*/
//#include <stdlib.h>
//#include <math.h>
//#include <assert.h>

/*PROTOTYPES*/
double **computeSimMat(Point*,int);
double **computeDegMat(double**,int);
double **computeNormSimMat(double**,double**,int);
double **UpdateH(double**,double**,int,int,double);

double **mulMat(double**,double**,int,int,int);
double **transpose(double**,int,int);
double **matDiff(double**,double**,int,int);
double eucDist(Point,Point);

/*DEFINITIONS*/
typedef struct P{
    double *coords;
    int dim;
}Point;