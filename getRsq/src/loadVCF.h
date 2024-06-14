//#include <config.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "htslib/tbx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"

#include <zlib.h>
#include <math.h>
#include "getLogBF.h"
#include "sort.h"

extern int exp_gt_gtdsgl;
extern int verbose_loadVCF;


#define FORMAT_GT 1
#define FORMAT_GL 2
#define FORMAT_AP 3
#define FORMAT_GP 4
#define FORMAT_AS 5
#define FORMAT_RD 6
#define FORMAT_BF 7
#define FORMAT_DS 8
#define FORMAT_PS 9
#define FORMAT_OTHER 0


#define VT_SNP 0
#define VT_INDEL 1
#define VT_SV 2
#define VT_OTHER 3

typedef struct{
        int VT;
        double RSQ;
        double AF;
        double CER;
}VCF_info;



int getBiVCFRsq(const char* fname, const char* reg, double** pds1, int* psamplesize, int* pnbivars, char** pchr, int** ppos, char*** prss, char*** pba0, char*** pba1, int** vtype, double** rsqs);
