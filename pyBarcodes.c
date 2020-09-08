#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <zlib.h>
#include <stdarg.h>
#include <unistd.h>
#include "khash.h"
#define MINCLUSTERS 1000000
#define THRESHOLD 0.005
KHASH_MAP_INIT_STR(32, uint32_t)

typedef struct CBCL CBCL;
struct CBCL {
    uint32_t nTiles;
    uint32_t *nClusters;
    uint32_t *uncompressedSize;
    uint32_t *compressedSize;
    uint64_t *offsets;
};

#define pyBarcodesVersion "0.2.0"

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

// NovaSeq 6000
FILE **openNovaSeq(char *basePath, int lane, int *cycles, int nCycles) {
    int i;
    char fname[16384];
    FILE **o = NULL;
    o = calloc(nCycles, sizeof(FILE*));  // Use only 1 surface
    if(!o) return NULL;

    for(i=0; i<nCycles; i++) {
        sprintf(fname, "%s/Data/Intensities/BaseCalls/L00%i/C%i.1/L00%i_1.cbcl", basePath, lane, cycles[i], lane);
        if(access(fname, F_OK) != 0) sprintf(fname, "%s/Data/Intensities/BaseCalls/L00%i/C%i.1/L00%i_2.cbcl", basePath, lane, cycles[i], lane);
        o[i] = fopen(fname, "r");
        if(!o[i] != Z_NULL) goto error;
    }

    return o;

error:
    for(i=0; i<=nCycles; i++) if(o[i] != NULL) fclose(o[i]);
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

void closeCBCLs(FILE **cbcls, int nBCLs) {
    int i;
    for(i=0; i<nBCLs; i++) fclose(cbcls[i]);
    free(cbcls);
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

//returns 1 on error, 0 on success
int initCBCL(CBCL *cbcl, uint32_t nTiles) {
    cbcl->nClusters = calloc(nTiles, sizeof(uint32_t));
    if(!cbcl->nClusters) return 1;
    cbcl->uncompressedSize = calloc(nTiles, sizeof(uint32_t));
    if(!cbcl->uncompressedSize) return 1;
    cbcl->compressedSize = calloc(nTiles, sizeof(uint32_t));
    if(!cbcl->compressedSize) return 1;
    cbcl->offsets = calloc(nTiles, sizeof(uint64_t));
    if(!cbcl->offsets) return 1;
    return 0;
}

//destroy a CBCL object
void destroyCBCL(CBCL *cbcl) {
    if(cbcl->nTiles) {
        if(cbcl->nClusters) free(cbcl->nClusters);
        if(cbcl->uncompressedSize) free(cbcl->uncompressedSize);
        if(cbcl->compressedSize) free(cbcl->compressedSize);
        if(cbcl->offsets) free(cbcl->offsets);
    }
    free(cbcl);
}

//Returns a NULL pointe on error
CBCL* loadCBCL(FILE *bcl) {
    uint8_t bitsPerBase, bitsPerQScore;
    uint16_t version;
    uint32_t i, j;
    uint32_t headerSize, QBins, nTiles;
    uint64_t blockOffset;
    CBCL *cbcl = NULL;
    cbcl = calloc(1, sizeof(CBCL));
    if(!cbcl) goto error;

    //CBCL header
    fseek(bcl, 0, SEEK_SET);
    if(fread((void*) &version, 2, 1, bcl) != 1) goto error;
    if(fread((void*) &headerSize, 4, 1, bcl) != 1) goto error;
    if(fread((void*) &bitsPerBase, 1, 1, bcl) != 1) goto error;
    if(fread((void*) &bitsPerQScore, 1, 1, bcl) != 1) goto error;
    if(fread((void*) &QBins, 4, 1, bcl) != 1) goto error;
    fseek(bcl, 8 * QBins, SEEK_CUR);  // Skip Q-value binning definition
    if(fread((void*) &nTiles, 4, 1, bcl) != 1) goto error;

    cbcl->nTiles = nTiles;
    if(initCBCL(cbcl, nTiles)) goto error;
    
    blockOffset = (uint64_t) ftell(bcl);
    blockOffset += 16 * nTiles + 1;  // Start just after the per-tile information
    for(i=0; i<nTiles; i++) {
        // Tile number
        if(fread((void*) &headerSize, 4, 1, bcl) != 1) goto error;
        // nClusters
        if(fread((void*) &headerSize, 4, 1, bcl) != 1) goto error;
        cbcl->nClusters[i] = headerSize;
        // uncompressedSize
        if(fread((void*) &headerSize, 4, 1, bcl) != 1) goto error;
        cbcl->uncompressedSize[i] = headerSize;
        // compressedSize
        if(fread((void*) &headerSize, 4, 1, bcl) != 1) goto error;
        cbcl->compressedSize[i] = headerSize;
        cbcl->offsets[i] = blockOffset;
        blockOffset += headerSize;
    }

    return cbcl;

error:
    if(cbcl) destroyCBCL(cbcl);
    return NULL;
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

char getCBCLBase(uint8_t *uncompressedTiles, uint32_t cluster) {
    uint8_t byte = uncompressedTiles[cluster/2];
    int offset = cluster % 2;

    if(offset) byte>>=4;
    switch(byte & 3) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 3:
            return 'T';
        case 2:
        default:
            return 'G';
    }
}

//TODO handle a return value of 0
uint32_t cbclTile(uint8_t **uncompressedTiles, int nCycles, uint32_t nClusters, khash_t(32) *h) {
    khiter_t k;
    char *seq = NULL;
    int cycle, ret;
    uint32_t cluster;

    seq = malloc(nCycles + 1);
    if(!seq) return 0;
    seq[nCycles] = '\0';

    for(cluster=0; cluster<nClusters; cluster++) {
        for(cycle=0; cycle<nCycles; cycle++) {
           seq[cycle] = getCBCLBase(uncompressedTiles[cycle], cluster);
        }

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

    free(seq);
    return nClusters;

error:
    if(seq) free(seq);
    return 0;
}


// Returns the number of clusters returning sequence, which is all of them
// This will load nCycles of data into memory for one tile!
int getCBCLSequence(CBCL** CBCLs, FILE **bcls, int tile, khash_t(32) *h, int nCycles) {
    int rv = 0, i;
    uint8_t **uncompressedTiles = NULL;
    uint8_t *compressedTile = NULL;
    uLongf destLen;
    uLong sourceLen;
    z_stream zs = {
        .zalloc = NULL,
        .zfree = NULL,
        .msg = NULL
    };

    uncompressedTiles = (uint8_t**) calloc(nCycles, sizeof(uint8_t*));
    if(!uncompressedTiles) goto error;

    for(i=0; i<nCycles; i++) {
        destLen = CBCLs[i]->uncompressedSize[tile];
        sourceLen = CBCLs[i]->compressedSize[tile];

        uncompressedTiles[i] = malloc(destLen);
        if(!uncompressedTiles[i]) goto error;

        compressedTile = malloc(sourceLen);
        if(!compressedTile) goto error;

        fseek(bcls[i], CBCLs[i]->offsets[tile], SEEK_SET);
        fread(compressedTile, sourceLen, 1,  bcls[i]);

        // Skip the 10 byte header
        zs.next_in = compressedTile + 10;
        zs.avail_in = sourceLen - 10;
        zs.next_out = uncompressedTiles[i];
        zs.avail_out = destLen;

        rv = inflateInit2(&zs, -15);
        rv = inflate(&zs, Z_FINISH);
        inflateEnd(&zs);
        free(compressedTile);
        compressedTile = NULL;
        if(rv != Z_STREAM_END) goto error;
        if(destLen != zs.total_out) goto error;
    }

    rv += (int) cbclTile(uncompressedTiles, nCycles, CBCLs[0]->nClusters[0], h);

    for(i=0; i<nCycles; i++) free(uncompressedTiles[i]);
    free(uncompressedTiles);

    return rv;

error:
    if(compressedTile) free(compressedTile);
    if(uncompressedTiles) {
        for(i=0; i<nCycles; i++) {
            if(uncompressedTiles[i]) free(uncompressedTiles[i]);
        }
        free(uncompressedTiles);
    }

    //TODO return -1 on error and handle that
    return 0;
}

//Return the number of clusters passing filter (up to 1 million), -1 on error
// No filter files, since the cbcl files have been filtered already
int CBCLProcess(FILE **bcls, khash_t(32) *h, int nCycles) {
    int i, good = 0;
    CBCL** CBCLs = NULL;

    CBCLs = calloc(nCycles, sizeof(CBCL));
    if(!CBCLs) goto error;

    //Open each file, reading each into a data structure with offsets and a vector of numbers of clusters
    for(i=0; i<nCycles; i++) {
        CBCLs[i] = loadCBCL(bcls[i]);
        if(!CBCLs[i]) goto error;
    }

    //Send each tile into getCBCLSequence
    for(i=0; i<CBCLs[0]->nTiles; i++) {
        good += getCBCLSequence(CBCLs, bcls, i, h, nCycles);
        if(good > MINCLUSTERS) break;
    }

    for(i=0; i<nCycles; i++) destroyCBCL(CBCLs[i]);
    free(CBCLs);

    return good;

error:
    if(CBCLs) {
        for(i=0; i<nCycles; i++) {
            if(CBCLs[i]) destroyCBCL(CBCLs[i]);
        }
        free(CBCLs);
    }
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

// Returns the number of values in *barcodes and *frequencies, which must both be free()d
// Unlike the other functions, this doesn't care about tiles since they're concatenated.
// Only a single side of each lane is used.
int handleNovaSeq(char *basePath, int lane, int nCycles, int *cycles, char ***barcodes, float **frequencies) {
    FILE **cbcls = NULL;
    uint32_t i;
    int rv, good, nBarcodes = 0;
    khiter_t k;

    //This hash will get reused until we've processed up to a million clusters
    khash_t(32) *h = kh_init(32);

    cbcls = openNovaSeq(basePath, lane, cycles, nCycles);
    if(!cbcls) goto error;

    good = CBCLProcess(cbcls, h, nCycles);
    if(good == -1) goto error;

    closeCBCLs(cbcls, nCycles);
    cbcls = NULL;

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

    return nBarcodes;

error:
    kh_destroy(32, h);
    if(cbcls) closeCBCLs(cbcls, nCycles);
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
    else if(strcmp(runType, "NovaSeq") == 0) nBarcodes = handleNovaSeq(basePath, lane, nCycles, cycles, &barcodes, &frequencies);
    else if(strcmp(runType, "HiSeq3000") == 0 || \
            strcmp(runType, "HiSeq4000") == 0 || \
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
