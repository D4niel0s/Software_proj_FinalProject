#ifndef SYMNMF_H
    #define SYMNMF_H
    #include "symnmf.h"
#endif

/*COMPUTE SIMILARITY MATRIX*/

/*Recieves data which is an array of n d-dimensional points*/
double **computeSimMat(Point *data, int n,int d){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);

        for(j=0; j<n; ++j){
            if(i != j){
                res[i][j] = exp((-0.5)*eucDist(data[i],data[j]));
            }else{
                res[i][j] = 0;
            }
        }
    }

    return res;
}

/*Recieves the nxn similarity matrix of the datapoints*/
double **computeDegMat(double **simMat, int n){
    double **res = (double **)malloc(sizeof(double *)*n);
    double sum;
    int i,j;


    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);
        sum = 0;

        for(j=0; j<n; ++j){ /*Compute d[i]*/
            sum += simMat[i][j];
        }

        for(j=0; j<n; ++j){
            if(i == j){
                res[i][j] = sum;
            }else{
                res[i][j] = 0;
            }
        }
    }

    return res;
}


double eucDist(Point a, Point b){
    int i;
    double res = 0;
    for(i=0; i<a.dim; ++i){
        res += pow(a.coords[i] - b.coords[i], 2);
    }

    return res;
}