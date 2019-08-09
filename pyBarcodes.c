#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <zlib.h>
#include <stdarg.h>
#include "khash.h"
#define MINCLUSTERS 1000000
#define THRESHOLD 0.005
KHASH_MAP_INIT_STR(32, uint32_t)

#define pyBarcodesVersion "0.1.0"

static PyObject *pyGetStats(PyObject *self, PyObject *args);

static PyMethodDef barcodesMethods[] = {
    {"getStats", (PyCFunction) pyGetStats, METH_VARARGS,
"Get a dictionary of barcodes seen and their frequencies.\n\
\n\
Required arguments:\n\
    path: The path to the flow cell (it should contain a Data directory).\n\
    runType: One of HiSeq3000, HiSeq2500, NextSeq or MiSeq.\n\
    cycles:  The cycles containing the barcodes.\n\
\n\
Optional arguments:\n\
    lane:    The lane number (defaults to 1).\n\
\n\
Returns:\n\
    A dictionary with barcodes as keys and fractional prevalence as values.\n\
\n\
>>> from pyBarcodes import getStats\n\
>>> getStats(`/data/180215_J00182_0064_AHNVNGBBXX', 'HiSeq3000', range([77, 92]), 5)\n\
{'GGCAGAAAGAGGATA': 4.001095771789551, 'AGGCATGTCTACTCT': 4.900794982910156,\n\
'GACTCCTTCTACTCT': 6.861392974853516, 'CCTGAGCTCTACTCT': 8.011892318725586,\n\
'GGCAGAATCTACTCT': 6.203693866729736, 'NNNNNNNNNNNNNNN': 14.651585578918457,\n\
'AAGGCGATCTACTCT': 9.15669059753418, 'GACTCCTAGAGGATA': 4.463195323944092,\n\
'CCTGAGCAGAGGATA': 7.6103925704956055, 'GTACTAGTCTACTCT': 6.675393104553223,\n\
'AGGCATGAGAGGATA': 4.3499956130981445, 'GTACTAGAGAGGATA': 5.3411946296691895,\n\
'AAGGCGAAGAGGATA': 8.967890739440918}\n"},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
struct pyBarcodesmodule_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct pyBarcodesmodule_state*)PyModule_GetState(m))

static PyModuleDef pyBarcodesmodule = {
    PyModuleDef_HEAD_INIT,
    "pyBarcodes",
    "A python module for extracting barcode information from BCL files.",
    -1,
    barcodesMethods,
    NULL, NULL, NULL, NULL
};
#endif

#if PY_MAJOR_VERSION >= 3
PyObject *PyString_FromString(const char *v) {
    return PyUnicode_FromStringAndSize(v, strlen(v));
}
#endif

//NextSeq 500/550 and MiniSeq
gzFile *openNextSeqBCLs(char *basePath, int *cycles, int nCycles) {
    int i;
    char fname[16384];
    gzFile *o = NULL;
    o = calloc(nCycles, sizeof(FILE*));
    if(!o) return NULL;

    for(i=0; i<nCycles; i++) {
        sprintf(fname, "%s/Data/Intensities/BaseCalls/L001/%04i.bcl.bgzf", basePath, cycles[i]);
        o[i] = gzopen(fname, "r");
        if(!o[i] != Z_NULL) goto error;
    }
    return o;

error:
    for(i=0; i<=nCycles; i++) if(o[i] != Z_NULL && o[i] != NULL) gzclose(o[i]);
    free(o);
    return NULL;
}

//HiSeq 2000/2500/3000/4000/X single tile
gzFile *openHiSeq(char *basePath, int lane, int tile, int *cycles, int nCycles) {
    int i;
    char fname[16384];
    gzFile *o = NULL;
    o = calloc(nCycles, sizeof(FILE*));
    if(!o) return NULL;

    for(i=0; i<nCycles; i++) {
        sprintf(fname, "%s/Data/Intensities/BaseCalls/L00%i/C%i.1/s_%i_%i.bcl.gz", basePath, lane, cycles[i], lane, tile);
        o[i] = gzopen(fname, "r");
        if(!o[i] != Z_NULL) goto error;
    }
    return o;

error:
    for(i=0; i<=nCycles; i++) if(o[i] != Z_NULL && o[i] != NULL) gzclose(o[i]);
    free(o);
    return NULL;
}

//MiSeq runs are the same as HiSeq, except the bcl files aren't compressed
gzFile *openMiSeq(char *basePath, int tile, int *cycles, int nCycles) {
    int i;
    char fname[16384];
    gzFile *o = NULL;
    o = calloc(nCycles, sizeof(FILE*));
    if(!o) return NULL;

    for(i=0; i<nCycles; i++) {
        sprintf(fname, "%s/Data/Intensities/BaseCalls/L001/C%i.1/s_1_%i.bcl", basePath, cycles[i], tile);
        o[i] = gzopen(fname, "r");
        if(!o[i] != Z_NULL) goto error;
    }
    return o;

error:
    for(i=0; i<=nCycles; i++) if(o[i] != Z_NULL && o[i] != NULL) gzclose(o[i]);
    free(o);
    return NULL;
}

//NextSeq 500/550 and MiniSeq
FILE *openFilterNextSeq(char *basePath) {
    char fname[16384];
    sprintf(fname, "%s/Data/Intensities/BaseCalls/L001/s_1.filter", basePath);
    return fopen(fname, "r");
}

//HiSeq 2000/2500/3000/4000/X and MiSeq, single tile
FILE *openFilterHiSeq(char *basePath, int lane, int tile) {
    char fname[16384];
    sprintf(fname, "%s/Data/Intensities/BaseCalls/L00%i/s_%i_%i.filter", basePath, lane, lane, tile);
    return fopen(fname, "r");
}

void closeBCLs(gzFile *bcls, int nBCLs) {
    int i;
    for(i=0; i<nBCLs; i++) gzclose(bcls[i]);
    free(bcls);
}

//Returns 1 on error, 0 on success
int getSequence(uint32_t cluster, gzFile *bcls, uint32_t nBCLs, char *seq) {
    uint8_t byte;
    uint32_t i;

    for(i=0; i<nBCLs; i++) {
        gzseek(bcls[i], cluster + 3, SEEK_SET);
        if(gzread(bcls[i], (void*) &byte, 1) != 1) return 1;
        if(byte == 0) {
            seq[i] = 'N';
        } else {
            switch(byte & 3) {
                case 0:
                    seq[i] = 'A';
                    break;
                case 1:
                    seq[i] = 'C';
                    break;
                case 2:
                    seq[i] = 'G';
                    break;
                default:
                    seq[i] = 'T';
                    break;
            }
        }
    }

    return 0;
}

//Return the number of clusters passing filter (up to 1 million), -1 on error
int commonProcess(FILE *filterFile, gzFile *bcls, khash_t(32) *h, int nCycles) {
    khiter_t k;
    int i, ret, good = 0;
    uint32_t nClusters;
    uint8_t byte;
    char *seq = NULL;

    seq = malloc(nCycles + 1);
    if(!seq) return -1;
    seq[nCycles] = '\0';

    //Read in the header
    if(fread((void*) &nClusters, 4, 1, filterFile) != 1) goto error;
    if(fread((void*) &nClusters, 4, 1, filterFile) != 1) goto error;
    if(fread((void*) &nClusters, 4, 1, filterFile) != 1) goto error;

    for(i=0; i<nClusters; i++) {
        if(fread((void*) &byte, 1, 1, filterFile) != 1) goto error;
        if(byte & 1) {
            good++;
            if(getSequence(i, bcls, nCycles, seq)) goto error;

            //increment the counter
            k = kh_get(32, h, seq);
            if(k == kh_end(h)) {
                k = kh_put(32, h, seq, &ret);
                kh_value(h, k) = 0;
                seq = malloc(nCycles + 1);
                if(!seq) goto error;
                seq[nCycles] = '\0';
            }
            kh_value(h, k)++;
        }
        if(good > MINCLUSTERS) break;
    }

    free(seq);
    return good;

error:
    if(seq) free(seq);
    return -1;
}

//Handle NextSeq 500/550 and MiniSeq runs, will only look at lane 1
//Returns the number of values in *barcodes and *frequencies, which must both be free()d
int handleNextSeq(char *basePath, int nCycles, int *cycles, char ***barcodes, float **frequencies) {
    FILE *filterFile = NULL;
    gzFile *bcls = NULL;
    uint32_t i;
    int good = 0, nBarcodes = 0;
    khiter_t k;

    khash_t(32) *h = kh_init(32);
    filterFile = openFilterNextSeq(basePath);
    if(!filterFile) goto error;

    bcls = openNextSeqBCLs(basePath, cycles, nCycles);
    if(!bcls) goto error;

    good = commonProcess(filterFile, bcls, h, nCycles);
    if(good == -1) goto error;

    //Count the number of barcodes that will be output
    for(k = kh_begin(h); k != kh_end(h); k++) {
        if(kh_exist(h, k)) {
            if(kh_value(h, k) >= THRESHOLD * good) nBarcodes++;
        }
    }

    *barcodes = malloc(nBarcodes * sizeof(char**));
    *frequencies = malloc(nBarcodes * sizeof(float));
    if(!*barcodes) goto error;
    if(!*frequencies) goto error;
    for(k = kh_begin(h), i = 0; k != kh_end(h); k++) {
        if(kh_exist(h, k)) {
            if(kh_value(h, k) >= THRESHOLD * good) {
                (*frequencies)[i] = (100. * kh_value(h, k)) / good;
                (*barcodes)[i++] = (char *) kh_key(h, k);
            } else {
                free((char*) kh_key(h, k));
            }
        }
    }

    fclose(filterFile);
    closeBCLs(bcls, nCycles);
    kh_destroy(32, h);

    return nBarcodes;

error:
    if(bcls) closeBCLs(bcls, nCycles);
    if(filterFile) fclose(filterFile);
    kh_destroy(32, h);
    return -1;
}

// Returns the number of values in *barcodes and *frequencies, which must both be free()d
int handleHiSeq(char *basePath, int lane, int nCycles, int maxSwath, int maxTile, int *cycles, char ***barcodes, float **frequencies) {
    FILE *filterFile = NULL;
    gzFile *bcls = NULL;
    char *seq;
    uint32_t i, good = 0;
    int side, swath, tile, tileNum, rv, nBarcodes = 0;
    khiter_t k;

    seq = malloc(nCycles + 1);
    if(!seq) return -1;
    seq[nCycles] = '\0';

    //This hash will get reused until we've processed up to a million clusters
    khash_t(32) *h = kh_init(32);

    for(side=1; side<3; side++) {
        for(swath=1; swath<=maxSwath; swath++) {
            for(tile=1; tile<=maxTile; tile++) {
                tileNum = 1000 * side + 100 * swath + tile;
                filterFile = openFilterHiSeq(basePath, lane, tileNum);
                if(!filterFile) goto error;

                if(maxSwath > 1) bcls = openHiSeq(basePath, lane, tileNum, cycles, nCycles);
                else bcls = openMiSeq(basePath, tileNum, cycles, nCycles);
                if(!bcls) goto error;

                rv = commonProcess(filterFile, bcls, h, nCycles);
                if(rv == -1) goto error;
                good += rv;

                fclose(filterFile);
                closeBCLs(bcls, nCycles);
                bcls = NULL;
                filterFile = NULL;
                if(good > MINCLUSTERS) break;
            }
            if(good > MINCLUSTERS) break;
        }
        if(good > MINCLUSTERS) break;
    }

    //Count the number of barcodes that will be output
    for(k = kh_begin(h); k != kh_end(h); k++) {
        if(kh_exist(h, k)) {
            if(kh_value(h, k) >= THRESHOLD * good) nBarcodes++;
        }
    }

    *barcodes = malloc(nBarcodes * sizeof(char**));
    *frequencies = malloc(nBarcodes * sizeof(float));
    if(!*barcodes) goto error;
    if(!*frequencies) goto error;
    for(k = kh_begin(h), i = 0; k != kh_end(h); k++) {
        if(kh_exist(h, k)) {
            if(kh_value(h, k) >= THRESHOLD * good) {
                (*frequencies)[i] = (100. * kh_value(h, k)) / good;
                (*barcodes)[i++] = (char *) kh_key(h, k);
            } else {
                free((char*) kh_key(h, k));
            }
        }
    }

    kh_destroy(32, h);
    free(seq);

    return nBarcodes;

error:
    kh_destroy(32, h);
    if(seq) free(seq);
    if(bcls) closeBCLs(bcls, nCycles);
    if(filterFile) fclose(filterFile);
    return -1;
}

/********************************************************************
 *
 * Begin python wrapping stuff
 *
 ********************************************************************/
static PyObject *pyGetStats(PyObject *self, PyObject *args) {
    char *basePath = NULL;
    char *runType = NULL;
    PyObject *listObj = NULL, *item = NULL, *rv = NULL, *key = NULL, *value = NULL;
    int lane = 1;
    int *cycles = NULL, nCycles, i, nBarcodes = -1;
    char **barcodes = NULL;
    float *frequencies = NULL;

    if(!(PyArg_ParseTuple(args, "ssO|i", &basePath, &runType, &listObj, &lane))) {
        PyErr_SetString(PyExc_RuntimeError, "You must supply at least a path, a run type and a list of cycles.");
        return NULL;
    }

    if(strcmp(runType, "NextSeq") != 0 && \
       strcmp(runType, "HiSeq2500") != 0 && \
       strcmp(runType, "HiSeq3000") != 0 && \
       strcmp(runType, "NovaSeq") != 0 && \
       strcmp(runType, "MiSeq") != 0) {
        PyErr_SetString(PyExc_RuntimeError, "The run type must be one of NextSeq, HiSeq2500, HiSeq3000, NovaSeq or MiSeq");
        return NULL;
    }

    if(lane < 1 || lane > 8) {
        PyErr_SetString(PyExc_RuntimeError, "You have specified an illegal lane (only values between 1 and 8 are acceptable for currently existing machines");
        return NULL;
    }

    //set up the bounds
    if(PySequence_Check(listObj)) nCycles = PySequence_Size(listObj);
    else nCycles = PySequence_Size(listObj);
    cycles = malloc(nCycles * sizeof(int));
    if(!cycles) {
        PyErr_SetString(PyExc_RuntimeError, "Ran out of memory!");
        return NULL;
    }
    for(i=0; i<nCycles; i++) {
        if(PySequence_Check(listObj)) item = PySequence_GetItem(listObj, i);
        else item = PyList_GET_ITEM(listObj, i);
        if(!PyLong_Check(item)) goto error;
        cycles[i] = (int) PyLong_AsLong(item);
    }

    if(strcmp(runType, "NextSeq") == 0) nBarcodes = handleNextSeq(basePath, nCycles, cycles, &barcodes, &frequencies);
    else if(strcmp(runType, "HiSeq3000") == 0 || \
            strcmp(runType, "HiSeq4000") == 0 || \
            strcmp(runType, "NovaSeq") == 0 || \
            strcmp(runType, "HiSeqX") == 0) nBarcodes = handleHiSeq(basePath, lane, nCycles, 2, 28, cycles, &barcodes, &frequencies);
    else if(strcmp(runType, "HiSeq2500") == 0 || \
            strcmp(runType, "HiSeq2000") == 0) nBarcodes = handleHiSeq(basePath, lane, nCycles, 2, 16, cycles, &barcodes, &frequencies);
    else if(strcmp(runType, "MiSeq") == 0) nBarcodes = handleHiSeq(basePath, 1, nCycles, 1, 19, cycles, &barcodes, &frequencies);
    if(nBarcodes < 0) printf("The number of barcodes is %i\n", nBarcodes);
    if(nBarcodes < 0) goto error;

    // Create a dictionary with barcodes as keys and frequencies as views
    rv = PyDict_New();
    if(!rv) printf("No new dictionary\n");
    if(!rv) goto error;
    if(nBarcodes) {
        for(i=0; i<nBarcodes; i++) {
            key = PyString_FromString(barcodes[i]);
            if(!key) goto error;
            value = PyFloat_FromDouble((double) frequencies[i]);
            if(!value) goto error;
            if(PyDict_SetItem(rv, key, value)) goto error;
            free(barcodes[i]);
        }
        free(barcodes);
        free(frequencies);
    }
    free(cycles);

    return rv;

error:
    PyErr_SetString(PyExc_RuntimeError, "Received an error while parsing the BCL files!");
    if(cycles) free(cycles);
    if(barcodes) free(barcodes);
    if(frequencies) free(frequencies);
    if(rv) Py_DECREF(rv);
    if(key) Py_DECREF(key);
    if(value) Py_DECREF(value);

    return NULL;
}

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_pyBarcodes(void) {
#else
PyMODINIT_FUNC initpyBarcodes(void) {
#endif
    PyObject *res;
    errno = 0; //just in case

#if PY_MAJOR_VERSION >= 3
    res = PyModule_Create(&pyBarcodesmodule);
    if(!res) return NULL;
#else
    res = Py_InitModule3("pyBarcodes", barcodesMethods, "A module for ");
#endif

    PyModule_AddStringConstant(res, "__version__", pyBarcodesVersion);

#if PY_MAJOR_VERSION >= 3
    return res;
#endif
}
