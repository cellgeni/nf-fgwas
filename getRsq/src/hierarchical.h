#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#include <cblas.h>

#include <zlib.h>
#include <math.h>

double* beta0;
double* gamma0;

typedef struct{
	long tid;  // thread id
        long npeaksperthread;  // number of peaks processed [tid*npeaks ... (tid+1)*npeaks-1]
	long npeaks;  // total peaks
	long nvars;
	long* cumloci;
	long* cumcats;
	long P;
	double* X;
	double* BF;

	double* beta;
	double* Pi;

	double* eta;
	double* pjk;
	double* z;
	double* Z1;
	double* w;
	double* y;
	
	double* Xt;
	double* Ie;
        double lkhd;
}HIERARCHICAL_MT;

