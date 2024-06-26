#ifndef SYMNMFMODULE_C /*Inclusion guard*/
#define SYMNMFMODULE_C

#define PY_SSIZE_T_CLEAN
/*#include <Python.h>*/
#include "C:\Users\Danik\AppData\Local\Programs\Python\Python311\include\Python.h"

#include "symnmf.c"

/* Function prototypes */
Point *parsePyData(PyObject *,int);
PyObject *buildPyMat(double **,int,int);
double **BuildCMat(PyObject *);

static PyObject *compSimMat(PyObject *,PyObject *);
static PyObject *compDegMat(PyObject *,PyObject *);
static PyObject *compNormMat(PyObject *,PyObject *);
static PyObject *OptimizeH(PyObject *,PyObject *);

/* Python functions */

static PyObject *compSimMat(PyObject *self, PyObject *args){
    PyObject *GivenData; 
    
    Point *data;
    int i,n;
    double **res;
    
    PyObject *OUT;

    if(!PyArg_ParseTuple(args, "O", &GivenData)){
        return NULL;
    }
    
    n = (int)PyList_Size(GivenData);
    
    data = parsePyData(GivenData, n);
    res = computeSimMat(data, n);
    OUT = buildPyMat(res, n,n);

    /*Free auxiliary memory allocations*/
    for(i=0; i<n; ++i){
        free(data[i].coords);
        free(res[i]);
    }
    free(data);
    free(res);

    return OUT;
}

static PyObject *compDegMat(PyObject *self, PyObject *args){
    double **sim, **deg;
    PyObject *OUT, *tmp;
    int n,i;

    tmp = compSimMat(self,args);
    n = (int)PyList_Size(tmp);

    sim = BuildCMat(tmp);
    deg = computeDegMat(sim, n);
    OUT = buildPyMat(deg, n,n);
    
    /*Free auxiliary memory allocations*/
    for(i=0; i<n; ++i){
        free(sim[i]);
        free(deg[i]);
    }
    free(sim);
    free(deg);

    return OUT;
}

static PyObject *compNormMat(PyObject *self, PyObject *args){
    double **sim, **deg, **W;
    PyObject *OUT, *tmp;
    int n,i;

    tmp = compSimMat(self,args);
    n = (int)PyList_Size(tmp);

    sim = BuildCMat(tmp);
    deg = computeDegMat(sim, n);
    W = computeNormSimMat(sim, deg, n);

    OUT = buildPyMat(W, n,n);
    
    /*Free auxiliary memory allocations*/
    for(i=0; i<n; ++i){
        free(sim[i]);
        free(deg[i]);
        free(W[i]);
    }
    free(sim);
    free(deg);
    free(W);

    return OUT;
}

static PyObject *OptimizeH(PyObject *self, PyObject *args){
    PyObject *GivenH, *GivenW;
    int i,n,k, DEF_MAX_ITER;
    double **H, **W,**OPT, EPS;

    /*Output*/
    PyObject *OUT;

    /*Default values are assigned as instructed in assignment*/
    DEF_MAX_ITER = 300;
    EPS = 0.0001;

    if(!PyArg_ParseTuple(args, "OO", &GivenH, &GivenW)){
        return NULL;
    }
    n = (int)PyList_Size(GivenH);
    k = (int)PyList_Size(PyList_GetItem(GivenH,0));

    H = BuildCMat(GivenH);
    W = BuildCMat(GivenW);
    
    OPT = Hoptimization(H,W,n,k,DEF_MAX_ITER, EPS);
    OUT = buildPyMat(OPT, n,k);

    /*free auxiliary memory allocations*/
    for(i=0;i<n;++i){
        free(H[i]);
        free(W[i]);
        free(OPT[i]);
    }
    free(H);
    free(W);
    free(OPT);

    return OUT;
}


/* Helper functions */

/*Converts python list[Point] into C Point array*/
Point *parsePyData(PyObject *GivenData, int n){
    Point *data;

    PyObject *arr;
    PyObject *OBJ;
    long i,j,d;

    data = (Point *)malloc(sizeof(Point)*n);
    assert(data);

    /*Get dimension of 1st point (assume all point have the same dimension)*/
    OBJ = PyList_GetItem(GivenData, 0);
    arr = PyObject_GetAttrString(OBJ, "dim");
    d = PyLong_AsLong(arr);

    for(i=0; i<n; ++i){
        OBJ = PyList_GetItem(GivenData, i);
        arr = PyObject_GetAttrString(OBJ, "coords");

        data[i].coords = (double *)malloc(sizeof(double)*d);
        assert(data[i].coords);
        data[i].dim = d;

        for(j=0; j<d; ++j){
            OBJ = PyList_GetItem(arr, j); /*OBJ is data[i].coords[j] from python*/
            data[i].coords[j] = PyFloat_AsDouble(OBJ);
        }
    }

    return data;   
}

/*Builds python matrix from given nxm matrix*/
PyObject *buildPyMat(double **mat, int n, int m){
    PyObject *arr,*OBJ, *OUT;
    int i,j;

    OUT = PyList_New(n);
    for(i=0; i<n; ++i){
        arr = PyList_New(m);
        for(j=0; j<m; ++j){
            OBJ = Py_BuildValue("d", mat[i][j]);
            PyList_SetItem(arr, j, OBJ);
        }
        PyList_SetItem(OUT, i, arr);
    }

    return OUT;
}

/*Build C double matrix from python matrix*/
double **BuildCMat(PyObject *mat){
    int n,m,i,j;
    PyObject *arr, *OBJ;
    double **res;

    n = (int)PyList_Size(mat);
    m = (int)PyList_Size(PyList_GetItem(mat,0));

    res =(double **)malloc(sizeof(double *)*n);
    assert(res);
    for(i=0;i<n;++i){
        res[i] = (double *)malloc(sizeof(double)*m);
        assert(res[i]);
    }

    for(i=0;i<n;++i){
        arr = PyList_GetItem(mat,i);
        for(j=0;j<m;++j){
            OBJ = PyList_GetItem(arr,j);
            res[i][j] = PyFloat_AsDouble(OBJ);
        }
    }

    return res;
}



/* Python module initialization */

static PyMethodDef symMeths[] = {
    {"sym",
        (PyCFunction)compSimMat,
        METH_VARARGS,
        PyDoc_STR("Recieves data (list[Point]) and returns similarity matrix\n")},
    {"ddg",
        (PyCFunction)compDegMat,
        METH_VARARGS,
        PyDoc_STR("Recieves data (list[Point]) and returns diagonal degree matrix\n")},
    {"norm",
        (PyCFunction)compNormMat,
        METH_VARARGS,
        PyDoc_STR("Recieves data (list[Point]) and returns laplacian matrix\n")},
    {"symnmf",
        (PyCFunction)OptimizeH,
        METH_VARARGS,
        PyDoc_STR("Recieves initial H and laplacian and returns optimized H\n")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symMeths
};

PyMODINIT_FUNC PyInit_symnmfmodule(void){
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if(!m){
        return NULL;
    }
    return m;
}

#endif /*SYMNMFMODULE_C*/