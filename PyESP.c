#include <Python.h>
#include "simulate.h"

extern VAR_SIM s;
extern VAR_EXP v;

PyObject *PyESPError;
char err_str[100];

static PyObject *simulate(PyObject *self, PyObject *args)
{
    int i, np;

    /* Call the external C function to run the simulation. */
    if ((np = _simulate()) < 0) {
        PyErr_SetString(PyESPError, "Simulation failed on malloc.");
        return NULL;
    }

    PyObject* py_Pot = PyList_New(np);
    PyObject* py_Cur = PyList_New(np);
    for (i = 0; i < np; i++) {
        PyList_SetItem(py_Pot, i, Py_BuildValue("d", s.Pot[i]/1000.));
        PyList_SetItem(py_Cur, i, Py_BuildValue("d", s.Cur[i]));
    }
    PyObject* py_PotCur = PyList_New(2);
    PyList_SetItem(py_PotCur, 0, py_Pot);
    PyList_SetItem(py_PotCur, 1, py_Cur);
    //Py_DECREF(py_Pot);
    //Py_DECREF(py_Cur);

    return py_PotCur;
}

static PyObject *addSpecies(PyObject *self, PyObject *args)
{
    double conc, diff_const; int ref;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "dd", &conc, &diff_const))
        return NULL;
    /* Call the external C function to run the simulation. */
    if ((ref = _add_species(conc, diff_const)) == 0) {
        sprintf(err_str, "Cannot add more than %d species.", MAXSPEC);
        PyErr_SetString(PyESPError, err_str);
        return NULL;
    }
    return PyLong_FromLong(ref);
}

static PyObject *addRedox(PyObject *self, PyObject *args)
{
    double e, ke, a;
    int ox, red, n, ref;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "iiiddd", &ox, &red, &n, &e, &ke, &a))
        return NULL;
    /* Call the external C function to run the simulation. */
    if ((ref = _add_redox(ox, red, n, e, ke, a)) == 0) {
        sprintf(err_str, "Cannot add more than %d redox steps.", MAXREDOX);
        PyErr_SetString(PyESPError, err_str);
        return NULL;
    }
    return PyLong_FromLong(ref);
}

static PyObject *addChemical(PyObject *self, PyObject *args)
{
    double kf, kb;
    int a, b, c, d, ref;
    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "iiiidd", &a, &b, &c, &d, &kf, &kb))
        return NULL;
    /* Call the external C function to run the simulation. */
    if ((ref = _add_chemical(a, b, c, d, kf, kb)) == 0) {
        sprintf(err_str, "Cannot add more than %d chemical reactions.", MAXCHIM);
        PyErr_SetString(PyESPError, err_str);
        return NULL;
    }
    return PyLong_FromLong(ref);
}

static PyObject *set_params(PyObject* self, PyObject* args, PyObject* kwargs)
{
    static char* keywords[] = {"CP", "CT", "ET", "IP", "Mode", "SI", "ST", "V1", "V2",
        "VD", "FP", "NC", "SC", "AM", "ncycle", "T2", "FR", "PH", "IR", "DL", "RU", "AP", "TE", "AR", "WE", NULL};

    // default values by copying what's already in v:
    int CP=v.CP, IP=v.IP, SI=v.SI, V1=v.V1, V2=v.V2, FP=v.FP, NC=v.NC, SC=v.SC, AM=v.AM, ncycle=v.ncycle, T2=v.T2,
        PH=v.PH, IR=v.IR, AP=v.AP, WE=v.WE, Mode=v.Mode;
    double CT=v.CT, ET=v.ET, ST=v.ST, VD=v.VD, FR=v.FR, DL=v.DL, RU=v.RU, TE=v.TE, AR=v.AR;
    // read optional arguments (the pipe character "|" indicates that all arguments after it, are optional)
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iddiiidiidiiiiiidiiddiddi", keywords,
            &CP, &CT, &ET, &IP, &Mode, &SI, &ST, &V1, &V2,
                &VD, &FP, &NC, &SC, &AM, &ncycle, &T2, &FR, &PH, &IR, &DL, &RU, &AP, &TE, &AR, &WE))
        return NULL;
    // set all arguments
    if (_set_params(CP, CT, ET, IP, Mode, SI, ST, V1, V2,
        VD, FP, NC, SC, AM, ncycle, T2, FR, PH, IR, DL, RU, AP, TE, AR, WE) < 0) {
        PyErr_SetString(PyESPError, "Changing parameter(s) failed");
        return NULL;
    }
    return Py_None;
}

static PyObject *setup(PyObject *self, PyObject *args)
{
    /* Call the external C function to run the simulation. */
    if (_setup() < 0) {
        PyErr_SetString(PyESPError, "Setup failed, probably on malloc");
        return NULL;
    }
    return Py_None;
}

static PyObject *destroy(PyObject *self, PyObject *args)
{
    /* Call the external C function to run the simulation. */
    if (_destroy() < 0) {
        PyErr_SetString(PyESPError, "Destroy failed");
        return NULL;
    }
    return Py_None;
}

static PyMethodDef PyESPMethods[] = {
    {"setup",        setup,        METH_NOARGS,  "Setup default simulation"},
    {"destroy",      destroy,      METH_NOARGS,  "Free all memory"},
    {"set_params",   set_params,   METH_VARARGS | METH_KEYWORDS, "Set parameter(s)"},
    {"simulate",     simulate,     METH_VARARGS, "Run the simulation"},
    {"addSpecies",   addSpecies,   METH_VARARGS, "Add species to mechanism"},
    {"addRedox",     addRedox,     METH_VARARGS, "Add redox step to mechanism"},
    {"addChemical",  addChemical,  METH_VARARGS, "Add chemical reaction to mechanism"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef PyESPModule = {
    PyModuleDef_HEAD_INIT,
    "PyESP",     /* m_name */
    "This module provides a python interface to the ESP24 Do_Simul routine.",  /* m_doc */
    -1,                  /* m_size */
    PyESPMethods,        /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_PyESP(void)
{
    PyObject *m;

    Py_Initialize();

    m = PyModule_Create(&PyESPModule);

    if (m == NULL)
        return NULL;

    PyESPError = PyErr_NewException("PyESP.error", NULL, NULL);
    Py_INCREF(PyESPError);
    PyModule_AddObject(m, "error", PyESPError);

    return m;
}
