/*INCLUDES*/
//#include <stdlib.h>
//#include <math.h>
//#include <assert.h>

/*PROTOTYPES*/
double **computeSimMat(Point*,int,int);
double **computeDegMat(double**,int);
double **computeNormSimMat(double**,double**,int);
double **mulMat(double**,double**,int);
double eucDist(Point,Point);

/*DEFINITIONS*/
typedef struct P{
    double *coords;
    int dim;
}Point;