#ifndef SYMNMF_C /*Inclusion guard*/
#define SYMNMF_C

#include "symnmf.h"

#include <stdio.h>
#include <string.h>

#define MAX_LINE_LEN 500

int main(int argc, char **argv){
    char *fName, *goal, buf[MAX_LINE_LEN];
    Point *data;
    FILE *fp;
    int N=0 ,d, i,j;
    double **SIM,**DDG, **W, **OUT;
    int Dflag=0, Wflag=0;

    if(argc != 3){
        fprintf(stderr, "An error has occurred\n");
        exit(0);
    }

    goal = (char *)malloc(sizeof(char)*strlen(argv[1]));
    fName = (char *)malloc(sizeof(char)*strlen(argv[2]));
    assert(goal);
    assert(fName);

    strcpy(goal, argv[1]);
    strcpy(fName, argv[2]);

    if(strcmp(goal, "sym") != 0 && strcmp(goal, "ddg") != 0 &&strcmp(goal, "norm") != 0){
        fprintf(stderr, "An error has occurred\n");
        exit(0);
    }
    
    fp = fopen(fName, "r");

    /* Get N,d from file */
    while(fgets(buf, MAX_LINE_LEN, fp) != NULL){
        d = countCommas(buf) + 1;
        N++;
    }
    rewind(fp);
    
    /* Assign values to data */
    data = (Point *)malloc(sizeof(Point) * N);
    assert(data);

    for(i=0; i<N; ++i){
        fgets(buf, MAX_LINE_LEN, fp);
        
        data[i].coords = (double *)malloc(sizeof(double) *d);
        data[i].dim = d;
        assert(data[i].coords);
        for(j=0; j<d; ++j){
            sscanf(buf, "%lf,", &(data[i].coords[j]));
            strcpy(buf, skipUntilComma(buf));
        }
    }

    /* Perform calculation depending on input */
    SIM = computeSimMat(data, N);
    OUT = SIM;
    if(strcmp(goal, "ddg") == 0 || strcmp(goal, "norm") == 0){
        DDG = computeDegMat(SIM, N);
        Dflag = 1;
        OUT = DDG;
    }
    if(strcmp(goal, "norm") == 0){
        W = computeNormSimMat(SIM, DDG, N);
        Wflag = 1;
        OUT = W;
    }

    /* Print output */
    for(i=0;i<N;++i){
        for(j=0;j<N;++j){
            printf("%.4f",OUT[i][j]);
            if(j != N-1){
                printf(",");
            }
        }
        printf("\n");
    }

    /* Free all allocations */
    for(i=0;i<N;++i){
        if(Dflag){
            free(DDG[i]);
        }
        if(Wflag){
            free(W[i]);
        }
        free(data[i].coords);
    }
    if(Dflag){
        free(DDG);
    }
    if(Wflag){
        free(W);
    }
    free(goal);
    free(fName);
    free(data);

    return 1;
}


int countCommas(char *str){
    int i;
    char *s;
    for(s=str, i=0; s[i];){
        if(s[i] == ','){
            i++;
        }else{
            s++;
        }
    }

    return i;
}

char *skipUntilComma(char s[]){
    while(*s != ','){s++;}
    return ++s;
}



/*Recieves data which is an array of n points*/
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

        for(j=0; j<n; ++j){ /*Compute d[i] and assign zeroes*/
            sum += simMat[i][j];

            if(i != j){
                res[i][j] = 0;
            }
        }
        res[i][i] = sum;
    }

    return res;
}

/*A is the similarity matrix, and D is the diagonal degree matrix*/
double **computeNormSimMat(double **A,double **D, int n){
    double **normD = (double **)malloc(sizeof(double *)*n);
    double **C, **res;
    int i,j;
    assert(normD);

    for(i=0; i<n; ++i){
        normD[i] = (double *)malloc(sizeof(double)*n);
        assert(normD[i]);

        for(j=0; j<n; ++j){
            if(D[i][j] != 0){
                normD[i][j] = 1.0/(sqrt(D[i][j]));
            }else{
                normD[i][j] = 0;
            }
        }
    }

    C = mulMat(normD,A, n,n,n);
    res = mulMat(C,normD, n,n,n);
    /*Free auxiliary allocated memory*/
    for(i=0; i<n; ++i){
        free(normD[i]);
        free(C[i]);
    }
    free(normD);
    free(C);

    return res;
}

/*H is initialized randomly as said in the assignment, W is the laplacian (W:nxn, Hnxk)*/
double **Hoptimization(double **H, double **W,int n, int k, int max_iter, double eps){
    double **prevH,**curH,**diff;
    int iter,i,j;
    double diffNorm = INT_MAX;

    iter = 0;
    prevH = (double **)malloc(sizeof(double *)*n);
    assert(prevH);
    for(i=0;i<n;++i){
        prevH[i] = (double *)malloc(sizeof(double)*k);
        assert(prevH[i]);

        for(j=0;j<k;++j){
            prevH[i][j] = H[i][j];
        }
    }

    while(iter<max_iter && diffNorm>=eps){
        curH = UpdateH(prevH,W, n,k,0.5);
        diff = matDiff(curH,prevH, n,k);

        diffNorm = squaredFrobeniusNorm(diff, n,k);

        /*Free previos H and the differrence matrix (they become unused)*/
        for(i=0; i<n; ++i){
            free(prevH[i]);
            free(diff[i]);
        }
        free(prevH);
        free(diff);

        prevH = curH;
        iter++;
    }

    return curH;
}

/*Updates H matrix according to step 1.4.2 in the assignment (H:nxk, W:nxn)*/
double **UpdateH(double **H, double **W, int n, int k,double beta){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;

    double **WH = mulMat(W,H, n,n,k);
    double **H_t = transpose(H, n,k);
    double **HH_t = mulMat(H,H_t, n,k,n);
    double **HH_tH = mulMat(HH_t,H, n,n,k);
    assert(res);

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



/*Returns squared frobenius norm of a real nxm matrix A*/
/*Squared frobenius norm is given by trace(A^t*A) = sum of squares of all entries of A*/
double squaredFrobeniusNorm(double **A, int n, int m){
    int i,j;
    double res = 0;

    for(i=0; i<n; ++i){
        for(j=0;j<m;++j){
            res += pow(A[i][j],2);
        }
    }
    
    return res;
}

/**
 * Multiplies two  matrices(A:nxm, B:mxl), return value is A*B
 * Allocates a new matrix that should be free by the caller (free each row+the outer pointer)
 */
double **mulMat(double **A,double **B, int n,int m, int l){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j,k;
    assert(res);

    for(i=0;i<n;++i){
        res[i] = (double *)malloc(sizeof(double)*l);
        assert(res[i]);
        for(j=0;j<l;++j){
            res[i][j] = 0;
        }
    }

    for(i=0; i<n; ++i){ /*Iterate rows of A*/
        for(j=0; j<l; ++j){ /*Iterate columns of B*/
            for(k=0; k<m; ++k){ /*multiply elements*/
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return res;
}

/**
 * Returns the transpose of a given nxm matrix.
 * Allocates a new matrix that should be free by the caller (free each row+the outer pointer)
 */
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

/**
 * Returns the differrence between two nxm matrices. return value is A-B
 * Allocates a new matrix that should be free by the caller (free each row+the outer pointer)
 */
double **matDiff(double **A, double **B, int n, int m){
    double **res = (double **)malloc(sizeof(double *)*n);
    int i,j;
    assert(res);

    for(i=0; i<n; ++i){
        res[i] = (double *)malloc(sizeof(double)*m);
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

#endif /*SYMNMF_C*/