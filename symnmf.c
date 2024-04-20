#ifndef SYMNMF_H
    #define SYMNMF_H
    #include "symnmf.h"
#endif

/*Recieves data which is an array of n d-dimensional points*/
double **computeSimMat(Point *data, int n){
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

    double **C = mulMat(normD,A, n,n,n);
    double **res = mulMat(C,normD, n,n,n);

    /*Free auxiliary allocated memory*/
    for(i=0; i<n; ++i){
        free(normD[i]);
        free(C[i]);
    }
    free(normD);
    free(C);

    return res;
}


/*(H:nxk, W:nxn)*/
double **UpdateH(double **H, double **W, int n, int k,double beta){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;
    assert(res);

    double **WH = mulMat(W,H, n,n,k);
    double **H_t = transpose(H, n,k);
    double **HH_t = mulMat(H,H_t, n,k,n);
    double **HH_tH = mulmat(HH_t,H, n,n,k);

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*k);
        assert(res[i]);

        for(j=0; j<k; ++j){
            res[i][j] = H[i][j] * ((1-beta)+beta*(WH[i][j]/HH_tH[i][j]));
        }

    }
    /*Free all auxiliary memory allocation*/
    for(i=0; i<n; ++i){
        if(i<k){
            free(H_t[i]);
        }
        free(WH[i]);
        free(HH_t[i]);
        free(HH_tH[i]);
    }
    free(WH);
    free(H_t);
    free(HH_t);
    free(HH_tH);

    return res;
}   

/*Returns frobenius norm of a real nxm matrix A*/
double squaredFrobeniusNorm(double **A, int n, int m){
    double **A_t = transpose(A, n,m); /*The Hermitian conjucate of A is also the transpose because A is a real matrix*/
    double **A_tA = mulMat(A_t,A, m,n,m);
    int i;
    double res = 0;

    for(i=0; i<m; ++i){
        res += A_tA[i][i];
    }
    
    /*Free auxiliary memory allocation*/
    for(i=0; i<m; ++i){
        free(A_t[i]);
        free(A_tA[i]);
    }
    free(A_t);
    free(A_tA);

    return res;
}

/*Multiplies two  matrices(A:nxm, B:mxl), returned value is a new matrix: A*B*/
double **mulMat(double **A,double **B, int n,int m, int l){
    double **res = (double **)malloc(sizeof(double *)*n);
    int sum;
    int i,j,k;
    assert(res);

    for(i=0; i<n; ++i){ /*Iterate rows of A*/
        res[i] = (double *)malloc(sizeof(double)*l);
        assert(res[i]);

        for(j=0; j<l; ++j){ /*Iterate columns of B*/
            sum = 0;
            for(k=0; k<m; ++k){ /*multiply elements*/
                sum += A[i][k] * B[k][j];
            }

            res[i][j] = sum;
        }
    }

    return res;
}

/*Returns the transpose of a given nxm matrix*/
double **transpose(double **A, int n, int m){
    double **res = (double **)malloc(sizeof(double *)*m);
    int i,j;
    assert(res);

    for(i=0; i<m; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);
        assert(res[i]);

        for(j=0; j<n; ++j){
            res[i][j] = A[j][i];
        }
    }
    return res;
}

/*Returns the differrence between two nxm matrices. (A-B)*/
double **matDiff(double **A, double **B, int n, int m){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;
    assert(res);

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*n);
        assert(res[i]);

        for(j=0; j<m; ++j){
            res[i][j] = A[i][j] - B[i][j];
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