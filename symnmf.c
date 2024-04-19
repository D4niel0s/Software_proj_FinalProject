#ifndef SYMNMF_H
    #define SYMNMF_H
    #include "symnmf.h"
#endif

/*COMPUTE SIMILARITY MATRIX*/

/*Recieves data which is an array of n d-dimensional points*/
double **computeSimMat(Point *data, int n,int d){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;
    assert(res);

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);
        assert(res[i]);

        for(j=0; j<n; ++j){
            if(i == j){
                res[i][j] = 0;
            }else{
                res[i][j] = exp((-0.5)*eucDist(data[i],data[j]));
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
    assert(res);

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);
        assert(res[i]);
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

/*A is the similarity matrix, and D is the diagonal degree matrix*/
double **computeNormSimMat(double **A,double **D, int n){
    double **normD = (double **)malloc(sizeof(double *)*n);
    int i,j;
    assert(normD);

    for(i=0; i<n; ++i){
        normD[i] = (double *)malloc(sizeof(double)*n);
        assert(normD[i]);

        for(j=0; j<n; ++j){
            if(i == j){
                normD[i][j] = 1.0/(sqrt(D[i][j]));
            }else{
                normD[i][j] = 0;
            }
        }
    }

    double **C = mulMat(normD,A,n);
    double **res = mulMat(C,normD,n);

    /*Free auxiliary allocated memory*/
    for(i=0; i<n; ++i){
        free(normD[i]);
        free(C[i]);
    }
    free(normD);
    free(C);

    return res;
}





/*Multiplies two nxn matrices, returned value is a new matrix: A*B*/
double **mulMat(double **A,double **B, int n){
    double **res = (double **)malloc(sizeof(double *)*n);
    int sum;
    int i,j,k;
    assert(res);

    for(i=0; i<n; ++i){ /*Iterate rows of A*/
        res[i] = (double *)malloc(sizeof(double)*n);
        assert(res[i]);

        for(j=0; j<n; ++j){ /*Iterate columns of B*/
            sum = 0;
            for(k=0; k<n; ++k){ /*multiply elements*/
                sum += A[i][k] * B[k][j];
            }

            res[i][j] = sum;
        }
    }

    return res;
}


/*Calculates euclidean distance between a and b, assumes matching dimensions*/
double eucDist(Point a, Point b){
    int i;
    double res = 0;
    for(i=0; i<a.dim; ++i){
        res += pow(a.coords[i] - b.coords[i], 2);
    }

    return res;
}