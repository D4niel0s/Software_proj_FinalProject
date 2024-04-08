/*INCLUDES*/
//#include <stdlib.h>
//#include <math.h>

/*PROTOTYPES*/
double eucDist(Point,Point);
double **computeSimMat(Point*,int,int);
double **computeDegMat(double**,int);

/*DEFINITIONS*/
typedef struct P{
    double *coords;
    int dim;
}Point;