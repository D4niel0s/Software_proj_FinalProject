#define PY_SSIZE_T_CLEAN
/*#include <Python.h>*/
#include "C:\Users\Danik\AppData\Local\Programs\Python\Python311\include\Python.h"

#include "symnmf.c"

static PyObject *compSimMat(PyObject *self, PyObject *args){
    PyObject *GivenData; /*Input*/
    
    /*Local Vars*/
    Point *data;
    long i,j,n,d;
    double **res;
    PyObject *OBJ;
    PyObject *arr;

    /*Output*/
    PyObject *OUT;

    if(!PyArg_ParseTuple(args, "Oi", &GivenData, &n)){
        return NULL;
    }
    
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

        for(j=0; j<d; ++j){
            OBJ = PyList_GetItem(arr, j); /*OBJ is data[i].coords[j] from python*/
            data[i].coords[j] = PyFloat_AsDouble(OBJ);
        }
    }

    res = computeSimMat(data, n);
    
    OUT = PyList_New(n);
    for(i=0; i<n; ++i){
        arr = PyList_New(n);
        for(j=0; j<n; ++i){
            OBJ = Py_BuildValue("d", res[i][j]);
            PyList_SetItem(arr, j, OBJ);
        }
        PyList_SetItem(OUT, i, arr);
    }

    /*Free auxiliary memory allocations*/
    for(i=0; i<n; ++i){
        free(data[i].coords);
        free(res[i]);
    }
    free(data);
    free(res);

    return OUT;
}

static PyMethodDef symMeths[] = {
    {"sym",
        (PyCFunction)compSimMat,
        METH_VARARGS,
        PyDoc_STR("")},
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