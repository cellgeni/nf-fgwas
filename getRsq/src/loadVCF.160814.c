#include "loadVCF.h"

int startWith(const char *pre, const char *str)
{
    return strncmp(pre, str, strlen(pre));
}

void clearAs(double* x, int n, double val){
    int i;
    for(i=0; i<n; i++){
        x[i] = val;
    }
}

void clearAs0(double* x, int n){
    int i;
    for(i=0; i<n; i++){
        x[i] = 0.0;
    }
}


double nk_var(double* x, double* y, long N){
    double mx, my;
    int i;
    mx=my=0.0;
    double res=0.0;
    for(i=0; i<N; i++){
        mx += x[i]/(double)N;
        my += y[i]/(double)N;
    }
    for(i=0; i<N; i++){
        res += (x[i]-mx)*(y[i]-my);
    }
    return res/(double)N;
};

double nk_cor(double* x, double* y, int n){
    return nk_var(x, y ,n) / sqrt(nk_var(x, x, n)) / sqrt(nk_var(y, y, n));
}

int gzfdscanf(gzFile f, double** x){// N x P matrix
    char c;
    int n=0;
    gzseek(f, 0L, SEEK_SET);
    int maxnc=0, nc=0;
    while((c=gzgetc(f)) != EOF){
        if(c=='\n'||c==' '||c=='\t'||c==','){
            n++;
            if(nc>maxnc){maxnc = nc;}
            nc=0;
        }else{
            nc++;
        }
    }
    //fprintf(stderr, "%d values read\n", n);
    int i=0;
    gzseek(f, 0L, SEEK_SET);
    (*x) = calloc(n, sizeof(double));
    char* cell; cell = (char*)calloc(maxnc+1, sizeof(char));
    nc=0;
    while((c=gzgetc(f)) != EOF){
        if(c=='\n'||c==' '||c=='\t'||c==','){
            cell[nc]='\0';
            sscanf(cell, "%lf", (*x)+i);
            i++;
            nc=0;
        }else{
            cell[nc++] = c;
        }
    }
    return n;
}



double nk_lm2(double* X, int* id, double* y, int n, double* work){
    clearAs0(work, 9);
    double* xtx = work;
    double* xy  = work + 3;
    double* ms  = work + 6;
    int i=0, j;
    double dn = (double)n;
    
    // mean
    for(i=0; i<n; i++){
        ms[0] += X[id[0]*n+i] / dn;
        ms[1] += X[id[1]*n+i] / dn;
        ms[2] += y[i] / dn;
    }
    // var
    for(i=0; i<n; i++){
        xtx[0] += (X[id[0]*n+i]-ms[0])*(X[id[0]*n+i]-ms[0]);
        xtx[1] += (X[id[0]*n+i]-ms[0])*(X[id[1]*n+i]-ms[1]);
        xtx[2] += (X[id[1]*n+i]-ms[1])*(X[id[1]*n+i]-ms[1]);
        
        xy[0]  += (X[id[0]*n+i]-ms[0])*(y[i]-ms[2]);
        xy[1]  += (X[id[1]*n+i]-ms[1])*(y[i]-ms[2]);
        xy[2]  += (y[i]-ms[2])        *(y[i]-ms[2]);
    }
    // inverse
    double det = xtx[0]*xtx[2]-xtx[1]*xtx[1];
    if(det<=0.0){return 0.0;}
    double tmp =  xtx[2]/det;
    xtx[2]     =  xtx[0]/det;
    xtx[0]     =  tmp;
    xtx[1]     = -xtx[1]/det;
    
    double res = 0.0;
    for(i=0; i<2; i++){
        for(j=0; j<2; j++){
            res += xy[i]*xy[j]*xtx[i+j];
        }
    }
    return sqrt(res/xy[2]);
}

int parseFormat(char* str, int* formatID){
    int i;
    int nfield=1;
    for(i=0; i<strlen(str); i++){if(str[i]==':'){nfield++;}}
    char format[3];
    
    int ns=0;
    for(i=0; i<nfield; i++){
        if(i<nfield-1){
            sscanf(str+ns, "%[^:]:", format);
        }else{
            sscanf(str+ns, "%s", format);
        }
        ns += strlen(format)+1;
        //fprintf(stderr, "%d %s ", ns, format);
        if(strcmp(format,"GT")==0){
            formatID[i]=FORMAT_GT;
        }else if(strcmp(format,"GL")==0){
            formatID[i]=FORMAT_GL;
        }else if(strcmp(format,"AP")==0){
            formatID[i]=FORMAT_AP;
        }else if(strcmp(format,"GP")==0){
            formatID[i]=FORMAT_GP;
        }else if(strcmp(format,"PP")==0){
            formatID[i]=FORMAT_GP;
        }else if(strcmp(format,"AS")==0){
            formatID[i]=FORMAT_AS;
        }else if(strcmp(format,"RD")==0){
            formatID[i]=FORMAT_RD;
        }else if(strcmp(format,"BF")==0){
            formatID[i]=FORMAT_BF;
        }else if(strcmp(format,"DS")==0){
            formatID[i]=FORMAT_DS;
        }else{
            formatID[i]=FORMAT_OTHER;
        }
    }
    return nfield;
}

int doseFormatExist(int* formatID, int nfield, int formatID1){
    int i;
    for(i=0; i<nfield; i++){
        if(formatID[i]==formatID1){return 1;}
    }
    return 0;
}

int parseInfo(char* str, char* infostr, VCF_info* vinfo){
    int i;
    int nfield=1;
    for(i=0; i<strlen(str); i++){if(str[i]==';'){nfield++;}}
    int ns=0;
    vinfo->VT=-100;
    vinfo->RSQ = -1.0;
    for(i=0; i<nfield; i++){
        sscanf(str+ns, "%[^;];", infostr);
        if(strcmp(infostr,"VT=SNP")==0){
            vinfo->VT=VT_SNP;
        }else if(strcmp(infostr,"VT=INDEL")==0){
            vinfo->VT=VT_INDEL;
        }else if(strcmp(infostr,"VT=SV")==0){
            vinfo->VT=VT_SV;
        }else if(startWith("RSQ=",infostr)==0){
            sscanf(infostr, "RSQ=%lf", &(vinfo->RSQ));
        }else if(startWith("AF=",infostr)==0){
            sscanf(infostr, "AF=%lf", &(vinfo->AF));
        }else if(startWith("IMP2=",infostr)==0){
            sscanf(infostr, "IMP2=%lf,%lf,%lf", &(vinfo->AF), &(vinfo->RSQ), &(vinfo->CER));
        }
        ns += strlen(infostr)+1;
    }
    return 0;
}



int isSnp(int* allen, int nal){
    int i;
    for(i=0; i<nal; i++){
        if(allen[i]>1){return 0;}
    }
    return 1;
}

void printV(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf ", x[i]);}
    fprintf(stderr, "%lf\n", x[i]);
}

double nk_dsum(double* x, int n, int ldx){
    int i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

int nk_isum(int* x, int n, int ldx){
    int i;
    int res=0;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

void gt2ap(int* gt, int n, double* ap){
    int i;
    double* p; 
    double* q;
    clearAs0(ap, 2*n);
    p = ap;
    q = ap+n;
    p[gt[0]] = 1.0;
    q[gt[1]] = 1.0;
}

void rsq2ap(int* gt, double rsq, double* afs, int n, int samplesize, double* ap){// all samples with afs
    int i;
    double tafs = 0.0;
    for(i=0; i<n; i++){
        tafs += afs[i]*(1.0-afs[i]);
    }
    double ep = (1.0-rsq) * tafs / (2.0*((double)(n-1)));
    double* p; 
    double* q;
    clearAs(ap, 2*n*samplesize, ep);
    for(i=0; i<samplesize; i++){
        p = ap+i*n*2;
        q = ap+i*n*2+n;
        p[gt[i*2+0]] = 1.0-ep*(double)(n-1);
        q[gt[i*2+1]] = 1.0-ep*(double)(n-1);
    }
}

void ds2ap(int* gt, double* ds, int n, double* ap){
    if(gt[0]==gt[1]){
        int k;
        for(k=0; k<n; k++){
            ap[k]=ap[k+n]=ds[k]/2.0;
        }
    }else{
        clearAs0(ap, n*2);
        double* p; p = ap;
        double* q; q = ap+n;
        p[gt[0]] = ds[gt[0]];
        q[gt[1]] = ds[gt[1]];
        if(p[gt[0]]>1.0){p[gt[0]]=1.0; q[gt[0]]=ds[gt[0]]-1.0;}
        if(q[gt[1]]>1.0){q[gt[1]]=1.0; p[gt[1]]=ds[gt[1]]-1.0;}
        
        double tp = 1.0 - p[gt[0]] - p[gt[1]];
        double tq = 1.0 - q[gt[0]] - q[gt[1]];
        
        int k;
        for(k=0; k<n; k++){
            if(k!=gt[0] && k!=gt[1]){
                p[k] = tp/(tp+tq)*ds[k];
                q[k] = tq/(tp+tq)*ds[k];
            }
        }
    }
}

void gl2ap(int* gt, double* gl, int n, double* ap, double* d){
    clearAs(ap, n*2, 0.1/((double)(n-1)));
    double* p; p = ap;
    double* q; q = ap+n;
    p[gt[0]] = 0.9;
    q[gt[1]] = 0.9;
    int itr, j, k, l;
    double denom, lkhd, lkhd0=-1.0e10;
    for(itr=0; itr<100; itr++){
        clearAs0(d, n*n);
        l=0;
        lkhd=0.0;
        for(k=0; k<n; k++){
            for(j=0; j<=k; j++){
                denom = p[j]*q[k] + p[k]*q[j];
                if(denom>1.0e-20){
                    lkhd     += gl[l]*log(denom);
                    d[j*n+k] += gl[l]*p[j]*q[k] / denom;
                    d[k*n+j] += gl[l]*p[k]*q[j] / denom;
                }
                l++;
            }
        }
        for(j=0; j<n; j++){
            p[j] = nk_dsum(d+j*n, n, 1);
            q[j] = nk_dsum(d+j,   n, n);
        }
        if(fabs(lkhd-lkhd0)<1.0e-8){ break; }else{ lkhd0=lkhd; }
    }
}

int choose(int n, int k){
    int i;
    int res = 1;
    for(i=n; i>=n-k+1; i--){
        res *= i;
    }
    for(i=k; i>=1; i--){
        res /= i;
    }
    return res;
}
int achoose(int n){
    int i;
    int res=0;
    for(i=1; i<=n-1; i++){
        res += choose(n, i);
    }
    return res/2;
}

int getCombAlk(int n, int k, int* v, int id, int idmax, char* a0, char* a1, char** ba0, char** ba1){
    if(k==n){
        if(id==-1 || id>=idmax){return id+1;}
        int i, i0, i1, allen=0;
        memcpy(ba0[id], a0, strlen(a0));
        i0=strlen(a0);
        i1=0;
        k=1; // kth allele for a1
        for(i=0; i<strlen(a1)+1; i++){
            if(a1[i]==',' || a1[i]=='\0'){
                if(v[k]==1){
                    if(i1==0){
                        memcpy(ba1[id]+i1, a1+(i-allen), allen);
                        i1 += allen;
                    }else{
                        ba1[id][i1++]=',';
                        memcpy(ba1[id]+i1, a1+(i-allen), allen);
                        i1 += allen;
                    }
                }else{
                    ba0[id][i0++]=',';
                    memcpy(ba0[id]+i0, a1+(i-allen), allen);
                    i0 += allen;
                }
                k++;
                allen=0;
            }else{
                allen++;
            }
        }
        ba0[id][i0]='\0';
        ba1[id][i1]='\0';
        
        return id+1;
    }
    v[k]=0;
    id = getCombAlk(n, k+1, v, id, idmax, a0, a1, ba0, ba1);
    v[k]=1;
    id = getCombAlk(n, k+1, v, id, idmax, a0, a1, ba0, ba1);
}


int getCombk(int n, int k, int* v, int id, int idmax, double* ds, double* bds, int samplesize){
    if(k==n){
        if(id==-1 || id>=idmax){return id+1;}
        int i;
        //printf("%d ", id); for(i=0; i<n; i++){printf("%d ", v[i]);} 
        for(i=0; i<n; i++){
            if(v[i]>0){
                bds[id*samplesize]+=ds[i];
            }
        }
        //printf("%lf\n", bds[id]);
        
        return id+1;
    }
    v[k]=0;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize);
    v[k]=1;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize);
}


int countFields(char* alStr, char sep){
    int i;
    int n=0;
    for(i=0; i<strlen(alStr); i++){
        if(alStr[i]==sep){
            n++;
        }
    }
    return n+1;
}

int parseGT(char* cell, int* gt){
    int nchar1;
    sscanf(cell, "%d|%d%n", gt, gt+1, &nchar1);
    return nchar1 + 1;
}

int parseDS(char* cell, int n, double* ds){
    int nchar1, nchar=0, i;
    for(i=1; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", ds+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", ds+i, &nchar1);
    
    ds[0] = 1.0 - nk_dsum(ds+1, n-1, 1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int skipFormat(char* cell){
    int l=strlen(cell);
    int i;
    for(i=0; i<l; i++){
        if(cell[i]=='\t' || cell[i]==':' || cell[i]=='\n' || cell[i]=='\0'){return i;}
    }
}

int parseGP(char* cell, int n, double* gl){
    int nchar1, nchar=0, i;
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(cell+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1 + 1;
    return nchar;
}

int parseGL(char* cell, int n, double* gl){
    int nchar1, nchar=0, i;
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(cell+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1 + 1;
    for(i=0; i<n*(n-1)/2+n; i++){
        gl[i] = pow(10.0, gl[i]);
    }
    return nchar;
}

int parseAP(char* cell, int n, double* ap){
    int nchar1, nchar=0, i;
    for(i=1; i<n; i++){
        sscanf(cell+nchar, "%lf,%n", ap+i, &nchar1);
        nchar += nchar1;
    }
    for(i=n+1; i<2*n-1; i++){
        sscanf(cell+nchar, "%lf,%n", ap+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", ap+i, &nchar1);
    
    ap[0] = 1.0 - nk_dsum(ap+1,   n-1, 1);
    ap[n] = 1.0 - nk_dsum(ap+n+1, n-1, 1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int parseAS(char* cell, int n, double* as){
    int nchar1, nchar=0, i;
    for(i=0; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", as+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", as+i, &nchar1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int parseCell(char* cell, int nal, int* gt, double* gl, double* ap, double* asc, double* ds, int* formatID){
    int i, offs=0, k=0;
    for(i=0; i<strlen(cell)+1; i++){
        if(cell[i]==':' || cell[i]=='\0'){
            if(formatID[k]==FORMAT_GT){
                offs += parseGT(cell+offs, gt);
            }else if(formatID[k]==FORMAT_GL){
                offs += parseGL(cell+offs, nal, gl);
            }else if(formatID[k]==FORMAT_GP){
                offs += parseGP(cell+offs, nal, gl);
            }else if(formatID[k]==FORMAT_AP){
                offs += parseAP(cell+offs, nal, ap);
            }else if(formatID[k]==FORMAT_AS){
                offs += parseAS(cell+offs, nal, asc);
            }else if(formatID[k]==FORMAT_DS){
                offs += parseDS(cell+offs, nal, ds);
            }else{
                offs += skipFormat(cell+offs);
            }
            k++;
        }
    }
    return offs;
}




int parseBody3(char* body, int n, int* gt, double* ds, double* gl, double* ap, double* d){
    int i;
    int nchar=0, nchar1;
    
    // GT
    sscanf(body+nchar, "%d|%d:%n", gt, gt+1, &nchar1);
    nchar += nchar1;
    
    // DS
    for(i=1; i<n-1; i++){
        sscanf(body+nchar, "%lf,%n", ds+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%lf:%n", ds+i, &nchar1);
    nchar += nchar1;
    
    // GL
    for(i=0; i<n*(n-1)/2+n-1; i++){
        sscanf(body+nchar, "%lf,%n", gl+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%lf%n", gl+i, &nchar1);
    nchar += nchar1;
    
    return nchar+1;
}

int parseBody1(char* body, int n, int* gt, double* ds, double* gl, double* ap, double* d){
    int nchar1;
    sscanf(body, "%d|%d%n", gt, gt+1, &nchar1);
    return nchar1+1;
}

int parseBody(char* body, int n, int* gt){
    int i;
    int nchar=0, nchar1;
    for(i=0; i<n-1; i++){
        sscanf(body+nchar, "%d|%d\t%n", gt+2*i, gt+2*i+1, &nchar1);
        nchar += nchar1;
    }
    sscanf(body+nchar, "%d|%d", gt+2*i, gt+2*i+1);
    return i+1;
}

void gt2dsgl(int* gt, int nal, double* ds, double* gl){
    int j, k, l;
    clearAs0(ds, nal);
    clearAs0(gl, choose(nal,2)+nal);
    ds[gt[0]]++;
    ds[gt[1]]++;
    l=0;
    for(k=0; k<nal; k++){
        for(j=0; j<=k; j++){
            if((j==gt[0] && k==gt[1]) || (k==gt[0] && j==gt[1])){gl[l]++; break;}
            l++;
        }
    }
}

int parseCigarPoint(int start, char* c1, char* seq, int* at, int K, char*** als, int** asc, int* nals, int** allen){
    //char* op; // operation
    //op = (char*)calloc(100, sizeof(char));
    char op[10];
    int len;// length in cigar
    int start0=start;
    int nseq=0;
    char* pc1;
    int nc1=0;
	int ali, ai, k;
    //pc1 = c1;
    int end;
	int flag=0;
	int flagin=0;
	int nchar=0;
    int allelematch;
    char prevop;
    while((sscanf(c1+nc1, "%d%[MIDNSH=PX]%n", &len, op, &nchar))>0){
        nc1 += nchar;
		flagin=0;
        if(op[0]=='M' || op[0]=='=' || op[0]=='X'){
            end = start+len-1;
			long nseq0=nseq;
			for(k=0; k<K; k++){
				nseq=nseq0;
                if(start <= at[k] && at[k] <= end){
					flagin++;
                    nseq += (at[k] - start);
                    for(ai=0; ai<nals[k]; ai++){
#ifdef INDEL
                        if(end-at[k]+1 >= allen[k][ai]){
                            allelematch=1;
                            for(ali=0; ali<allen[k][ai]; ali++){
                                if(seq[nseq+ali]!=als[k][ai][ali]){
                                    allelematch=0;
                                    break;
                                }
                            }
                        }else{// allele exceeds match region
                            allelematch=0;
                        }
                        asc[k][ai] += allelematch;
#else
                        if(seq[nseq]==als[k][ai][0] && allen[k][ai]==1){
                            asc[k][ai]++;
                        }
#endif
                    }
                }
			}
            start += len;
            nseq = nseq0+len;
        }else if(op[0]=='D'){
#ifdef INDEL
            end = start+len-1;
            for(k=0; k<K; k++){
                if(start-1 == at[k] && end < at[k]+allen[k][0]){
                    fprintf(stderr, "%s\n", seq+nseq-1);
					flagin++;
                    for(ai=1; ai<nals[k]; ai++){// only alt allele(s)
                        if(allen[k][0]-allen[k][ai]==len){
                            //fprintf(stderr, "%d %d\n", allen[k][0], allen[k][ai]);
                            asc[k][ai]++;
                        }
                    }
                }
			}
#endif
            start += len;
        }else if(op[0]=='N' || op[0]=='P'){
            start += len;
        }else{// S, H and I
            nseq += len;
        }
		if(flagin>0){flag++;}
        prevop = op[0];
    }
    
    return flag;
}

void countAS(const char* fname, const char* reg){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nal;
    
    // count row num
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        body = str.s+nchar;
        nvars++;
    }
    tbx_itr_destroy(itr);
    
    // load VCF
    char*** als; als   = (char***)calloc(nvars, sizeof(char**));
    int** allen; allen =   (int**)calloc(nvars, sizeof(int*));
    int* poss;   poss  =    (int*)calloc(nvars, sizeof(int));
    int* nals;   nals  =    (int*)calloc(nvars, sizeof(int));
    int** asc;   asc   =   (int**)calloc(nvars, sizeof(int*));
    
    double* gl; gl = (double*)calloc(210, sizeof(double));
    double* ds; ds = (double*)calloc(20, sizeof(double));
    double* ap; ap = (double*)calloc(40, sizeof(double));
    double* d;  d  = (double*)calloc(400, sizeof(double));
    int*    gt; gt = (int*)calloc(2, sizeof(int));
    char sep;
    int l=0;
    itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        
        poss[l] = pos;
        
        nal = countFields(a1, ',')+1;
        nals[l] = nal;
        
        // memoty alloc
        allen[l] = (int*)calloc(nal, sizeof(int));
        als[l]   = (char**)calloc(nal, sizeof(char*));
        asc[l]   = (int*)calloc(nal, sizeof(int));
        
        // get allele length and copy
        allen[l][0] = strlen(a0);
        als[l][0]   = (char*)calloc(allen[l][0]+1, sizeof(char));
        strcpy(als[l][0], a0);
        k=1;
        for(i=0; i<strlen(a1)+1; i++){
            if(a1[i]==','){
                als[l][k]   = (char*)calloc(allen[l][k]+1, sizeof(char));
                memcpy(als[l][k], a1+(i-allen[l][k]), allen[l][k]);
                k++;
            }else if(a1[i]=='\0'){
                als[l][k]   = (char*)calloc(allen[l][k]+1, sizeof(char));
                memcpy(als[l][k], a1+(i-allen[l][k]), allen[l][k]);
            }else{
                allen[l][k]++;
            }
        }
        
        //fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%s\t%s\tAS\t%d\t", chr, poss[l], rs, a0, a1, qual, filter, nals[l]);
        
        //for(k=0; k<nals[l]; k++){
        //    fprintf(stdout, "%s %d ", als[l][k], allen[l][k]);
        //}
        //fprintf(stdout, "\n");
        
        l++; // variant ID
    }
    tbx_itr_destroy(itr);
    
    // parse fragment
    int left, right;
    char* seq; seq = (char*)calloc(1000, sizeof(char));
	char* c1;  c1 =  (char*)calloc(1000, sizeof(char));
	int pivot=0;
	int isAS, numOfFrags=0, numOfASFrags=0, nvInFrag;
	while(scanf("%s\t%d\t%d\t%s\t%s", chr, &left, &right, c1, seq)!=EOF){
		nvInFrag=0;
		isAS=0;
		for(l=pivot; l<nvars; l++){
			if(poss[l]>right){break;}
			if(poss[l]>=left){nvInFrag++;}else{pivot=l+1;}
		}
		if(nvInFrag>0){
			isAS = parseCigarPoint(left, c1, seq, poss+pivot, nvInFrag, als+pivot, asc+pivot, nals+pivot, allen+pivot);
		}
		if(isAS>0){
			numOfASFrags++;
		}
		numOfFrags++;
	}
    for(l=0; l<nvars; l++){
        //printf("%d\t%d\t%d\t", poss[l], nals[l], nk_isum(allen[l], nals[l], 1));
        // only snp counts
#ifdef INDEL
        
#else
        if(isSnp(allen[l], nals[l])==0){for(k=0; k<nals[l]; k++){ asc[l][k]=0.0; }}
#endif
        for(k=0; k<nals[l]; k++){
            if(k<nals[l]-1){sep=',';}else{sep='\n';}
            printf("%d%c", asc[l][k], sep);
		}
	}
	fprintf(stderr, "%d out of %d overlapped with SNPs\n", numOfASFrags, numOfFrags);
    
    
    
    
}



void createDataBase(const char* fname, const char* reg){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nrs;
    char* rs1; rs1 = (char*)calloc(1000, sizeof(char));
    int rspos;
    char* rschr; rschr = (char*)calloc(1000, sizeof(char));
    int geta;
    char* num; num = (char*)calloc(100, sizeof(char));
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        nrs = countFields(rs, ';');
        if(nrs==1){
            if(strlen(rs)>3){
                sscanf(rs, "%[^0-9]%d", rschr, &rspos);
                geta=0;
                if(rspos>500000000){geta = sprintf(num, "%d", rspos/100000000); rspos = rspos%100000000;}
                rs[4+geta] = '\0';
                fprintf(stdout, "%s\t%d\t%s\t%d\t%s\t%s\n", rs, rspos, chr, pos, a0, a1);
            }
        }else{
            int rslen=0;
            for(i=0; i<strlen(rs)+1; i++){
                if(rs[i]==';' || rs[i]=='\0'){
                    if(rslen>3){
                        memcpy(rs1, rs+(i-rslen), rslen);
                        rs1[rslen]='\0';
                        sscanf(rs1, "%[^0-9]%d", rschr, &rspos);
                        geta=0;
                        if(rspos>500000000){geta = sprintf(num, "%d", rspos/100000000); rspos = rspos%100000000;}
                        rs1[4+geta] = '\0';
                        fprintf(stdout, "%s\t%d\t%s\t%d\t%s\t%s\n", rs1, rspos, chr, pos, a0, a1);
                    }
                    rslen=0;
                }else{
                    rslen++;
                }
            }
        }
    }
}



int nrowBed(const char* fname, const char* reg){
    
    int i;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos1, pos2;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    int nrow=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d", chr, &pos1, &pos2);
        nrow++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}

int loadBed(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2){
    
    int i;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int nrow = nrowBed(fname, reg);
    
    int gstart=270000000, gend=0;
    (*pchr)=(char*)calloc(1000, sizeof(char));
    (*ppos1)=(int*)calloc(nrow, sizeof(int));
    (*ppos2)=(int*)calloc(nrow, sizeof(int));
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    i=0;
    int nchar;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d%n", *pchr, (*ppos1)+i, (*ppos2)+i, &nchar);
        i++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}





int nrowVCF(const char* fname, const char* reg, int binarize, int* pnvars, int* psamplesize){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nal;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        
        if(nvars==0){nsamples = countFields(body, '\t');}
        
        nvars++;
        nal = countFields(a1, ',')+1;
        
        nbivars += achoose(nal);
    }
    tbx_itr_destroy(itr);
    (*pnvars) = nvars;
    (*psamplesize) = nsamples;
    if(binarize>0){
        return nbivars;
    }else{
        return nvars;
    }
}


int getBiVCF(const char* fname, const char* reg, double** pds1, int* psamplesize, int* pnbivars, char** pchr, int** ppos, char*** prss, char*** pba0, char*** pba1, int** vtype){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    (*pchr)   =(char*)calloc(1000, sizeof(char));
    char* chr; chr = (*pchr);
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info=(char*)calloc(1000, sizeof(char));
    VCF_info vinfo;
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    double* gl; gl = (double*)calloc(210, sizeof(double));
    double* ds; ds = (double*)calloc(20, sizeof(double));
    double* asc; asc = (double*)calloc(20, sizeof(double));
    double* ap; ap = (double*)calloc(40, sizeof(double));
    int*    gt; gt = (int*)calloc(2, sizeof(int));
    double* d;  d  = (double*)calloc(400, sizeof(double));
    
    int* formatID; formatID=(int*)calloc(100, sizeof(int));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0, nbivars1;
    int nal;
    int nfields;
    
    int ncol_ds1 = nrowVCF(fname, reg, 1, &nvars, &nsamples);
    (*pnbivars) = ncol_ds1;
    (*psamplesize) = nsamples;
    (*pds1) = (double*)calloc(nsamples*ncol_ds1, sizeof(double));
    (*pba0) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*pba1) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*prss) = (char**)calloc(ncol_ds1, sizeof(char*));
    (*ppos) = (int*)calloc(ncol_ds1, sizeof(int));
    (*vtype) = (int*)calloc(ncol_ds1, sizeof(int));
    
    int* work; work=(int*)calloc(100, sizeof(int));
    char* cwork; cwork=(char*)calloc(1000, sizeof(char));
    int offs;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    int l=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        body = str.s+nchar;
        nal = countFields(a1, ',')+1;
        
        nfields = parseFormat(format, formatID);
        parseInfo(info, cwork, &vinfo);
        
        nbivars1 = achoose(nal);
        for(i=nbivars; i<nbivars+nbivars1; i++){
            (*pba0)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*pba1)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*prss)[i] = (char*)calloc(strlen(rs)+1, sizeof(char));
            strcpy((*prss)[i], rs);
            (*ppos)[i] = pos;
            (*vtype)[i] = vinfo.VT;
        }
        getCombAlk(nal, 0, work, -1, nbivars1, a0, a1, (*pba0)+nbivars, (*pba1)+nbivars);
        for(i=nbivars; i<nbivars+nbivars1; i++){
            if((*vtype)[i]<0){
                if(strlen((*pba0)[i])==1 && strlen((*pba1)[i])==1){
                    (*vtype)[i] = VT_SNP;
                }else{
                    (*vtype)[i] = VT_INDEL;
                }
            }
        }
        offs = 0;
        for(i=0; i<nsamples; i++){
            offs += parseCell(body+offs, nal, gt, gl, ap, asc, ds, formatID);
            if(doseFormatExist(formatID, nfields, FORMAT_DS)==0){
                if(doseFormatExist(formatID, nfields, FORMAT_GP)==1){
                    gl2ap(gt, gl, nal, ap, d);;
                }
                for(k=0; k<nal; k++){
                    ds[k] = ap[k]+ap[k+nal];
                }
            }
            getCombk(nal, 0, work, -1, nbivars1, ds, (*pds1)+nbivars*nsamples+i, nsamples);
        }
        nbivars += nbivars1;
        l++;
    }
    tbx_itr_destroy(itr);
    
    
    for(i=0; i<nbivars; i++){
        //fprintf(stderr, "%s %s %s\n", rss[i], ba0[i], ba1[i]);
    }
    
}

int getRsq(const char* fname, const char* reg1, const char* reg2, const char* fid, int ld2){
    double* ds1;
    double* ds2;
    int nbivars1, nbivars2, samplesize;
    char* chr;
    int* pos1;
    int* pos2;
    char** ba0_1;
    char** ba1_1;
    char** ba0_2;
    char** ba1_2;
    char** rss1;
    char** rss2;
    int* vt1;
    int* vt2;
    
    double* work; work=(double*)calloc(9, sizeof(double));
    int id[2];
    
    getBiVCF(fname, reg1, &ds1, &samplesize, &nbivars1, &chr, &pos1, &rss1, &ba0_1, &ba1_1, &vt1);
    getBiVCF(fname, reg2, &ds2, &samplesize, &nbivars2, &chr, &pos2, &rss2, &ba0_2, &ba1_2, &vt2);
    
    int i, j, k;
    
    
    for(i=0; i<nbivars2; i++){
        for(j=0; j<nbivars1; j++){
            printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t", fid, rss2[i], ba0_2[i], ba1_2[i], rss1[j], ba0_1[j], ba1_1[j]);
            printf("%lf\n", nk_cor(ds1+samplesize*j, ds2+samplesize*i, samplesize));
        }
    }
   
    if(ld2>0){
        for(i=0; i<nbivars2; i++){
            for(j=0; j<nbivars1-1; j++){
                for(k=j+1; k<nbivars1; k++){
                    id[0] = j; id[1] = k;
                    printf("%s\t%s\t%s\t%s\t%s;%s\t%s;%s\t%s;%s\t", fid, rss2[i], ba0_2[i], ba1_2[i], rss1[j], rss1[k], ba0_1[j], ba0_1[k], ba1_1[j], ba1_1[k]);
                    printf("%lf\n", nk_lm2(ds1, id, ds2+samplesize*i, samplesize, work));
                }
            }
        }
    }
    
}


void expandVCF(const char* fname, const char* reg){
    
    int i, k;
    
    verbose_loadVCF=0;
    
    char* regchr; regchr=(char*)calloc(1000, sizeof(char));
    int regstart, regend;
    sscanf(reg, "%[^:]:%d-%d", regchr, &regstart, &regend);
    
    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) fprintf(stderr, "Could not read tabixed file %s\n", fname);
    //enum htsExactFormat format = hts_get_format(fp)->format;
    
    char *fnidx = calloc(strlen(fname) + 5, 1);
    strcat(strcpy(fnidx, fname), ".tbi");
    
    //regidx_t *reg_idx = NULL;
    
    tbx_t *tbx = tbx_index_load(fnidx);
    if ( !tbx ) fprintf(stderr, "Could not load .tbi index of %s\n", fnidx);
    
    kstring_t str = {0,0,0};
    
    //int nseq;
    //const char **seq = NULL;
    //if ( reg_idx ) seq = tbx_seqnames(tbx, &nseq);
    
    int gstart=270000000, gend=0;
    char* chr;    chr   =(char*)calloc(1000, sizeof(char));
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(1000, sizeof(char));
    char* a1;     a1    =(char*)calloc(1000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nal;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        
        if(nvars==0){nsamples = countFields(body, '\t');}
        
        nvars++;
        nal = countFields(a1, ',')+1;
        
        for(k=1; k<=nal-1; k++){
            nbivars += choose(nal, k);
        }
    }
    nbivars /= 2;
    tbx_itr_destroy(itr);
    
    itr = tbx_itr_querys(tbx, reg);
    double* gl; gl = (double*)calloc(210, sizeof(double));
    double* ds; ds = (double*)calloc(20, sizeof(double));
    double* ap; ap = (double*)calloc(40, sizeof(double));
    double* d;  d  = (double*)calloc(400, sizeof(double));
    int*    gt; gt = (int*)calloc(2, sizeof(int));
    char sep;
    
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        nal = countFields(a1, ',')+1;
        
        //parseBody(body, nsamples, gt);
        int nchar1=0;
        if(exp_gt_gtdsgl>0){
            fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s:DS:GP\t", chr, pos, rs, a0, a1, qual, filter, info, format);
        }else{
            fprintf(stdout, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", chr, pos, rs, a0, a1, qual, filter, info, format);
        }
        for(i=0; i<nsamples; i++){
            if(exp_gt_gtdsgl>0){
                nchar1 += parseBody1(body+nchar1, nal, gt, ds, gl, ap, d);
            }else{
                nchar1 += parseBody3(body+nchar1, nal, gt, ds, gl, ap, d);
            }
            // separator
            if(i<nsamples-1){sep='\t';}else{sep='\n';};
            //fprintf(stdout, "%d|%d%c", gt[0], gt[1], sep);
            
            // GT
            fprintf(stdout, "%d|%d:", gt[0], gt[1]);
            
            if(exp_gt_gtdsgl>0){
                gt2dsgl(gt, nal, ds, gl);
            }
            
            //gl2ap(gt, gl, nal, ap, d);
            //for(k=0; k<nal; k++){
            //    ds[k]=ap[k]+ap[k+nal];
            //}
            
            // DS
            for(k=1; k<nal-1; k++){
                fprintf(stdout, "%.2lf,", ds[k]);
            }
            fprintf(stdout, "%.2lf:", ds[k]);
            
            // GL
            for(k=0; k<choose(nal,2)+nal-1; k++){
                fprintf(stdout, "%.2lf,", gl[k]);
            }
            fprintf(stdout, "%.2lf%c", gl[k], sep);
        }
    }
    tbx_itr_destroy(itr);
}





int lm(int argc, char** argv){
    
    int i, j;
    const char* fname;  fname  = argv[1]; // VCF
    const char* reg;    reg    = argv[2];
    const char* fnamey; fnamey = argv[3]; // fpkm
    int fid = atoi(argv[4]);              // start point of fpkm, starting from 1
    
    int randomperm=0;
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-r")==0){randomperm=1;break;}}
    
    int start, end;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-s")==0){start=atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "-e")==0){end  =atoi(argv[i+1]); break;}}
    int tss = (start+end)/2;
    
    char* peakbed=NULL;
    int npeaks=0;
    int* pos1bed;
    int* pos2bed;
    int* pcents;
    int fid_bed=0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--peak-bed")==0){
        peakbed=argv[i+1]; 
        char* chrbed;
        npeaks = loadBed(peakbed, reg, &chrbed, &pos1bed, &pos2bed);
        pcents = (int*)calloc(npeaks, sizeof(int));
        for(j=0; j<npeaks; j++){
            pcents[j] = (pos1bed[j]+pos2bed[j])/2;
            if(pos1bed[j]<=tss && tss<=pos2bed[j]){
                fid_bed = j;
            }
        }
        break;
    }}
    
    
    double* ds;
    int nbivars, samplesize;
    char* chr;
    int* pos;
    char** ba0;
    char** ba1;
    char** rss;
    int* vt;
    
    double* work; work=(double*)calloc(9, sizeof(double));
    
    getBiVCF(fname, reg, &ds, &samplesize, &nbivars, &chr, &pos, &rss, &ba0, &ba1, &vt);
    
    double maf0 = 0.05;
    
    double* af; af=(double*)calloc(nbivars, sizeof(double));
    double* rsq; rsq=(double*)calloc(nbivars, sizeof(double));
    double* bf; bf=(double*)calloc(nbivars, sizeof(double));
    int* loccat;  loccat  =(int*)calloc(nbivars, sizeof(int));
    int* loccat2; loccat2 =(int*)calloc(nbivars, sizeof(int));
    
    // outcome
    FILE* fy;
    fy = fopen(fnamey, "rb");
    double* y;
    y = (double*)calloc(samplesize, sizeof(double));
    fseek(fy, samplesize*(fid-1)*sizeof(double), SEEK_SET);
    int fread_info = fread(y, sizeof(double), samplesize, fy);
    
    srand((unsigned)(time(NULL)+getpid()));
    if(randomperm>0){
        randomise(y, samplesize);
    }
    
    
    int k, maxid;
    double maxbf = -1.0e10;
    int ntested=0;
    
    for(j=0; j<nbivars; j++){
        af[j]  = nk_dsum(ds+j*samplesize, samplesize, 1)/2.0/(double)samplesize;
        rsq[j] = nk_var(ds+j*samplesize, ds+j*samplesize, samplesize)/af[j]/(1.0-af[j])/2.0;
        /*if(abs(pos[j]-tss)<10000){
            loccat[j] = 3;
        }else if(abs(pos[j]-tss)<50000){
            loccat[j] = 4;
        }else if(abs(pos[j]-tss)<100000){
            loccat[j] = 5;
        }else{
            loccat[j] = 0;
        }*/
        for(k=0; k<npeaks; k++){
            if(pos1bed[k]<pos[j]&&pos[j]<pos2bed[k]){loccat[j] = 2; break;}
        }
        if(start<pos[j]&&pos[j]<end){
            loccat[j] = 1;
        }
        if(af[j]>maf0 && af[j]<1.0-maf0 && rsq[j]>0.3){
            bf[j] = getLogBF(ds+j*samplesize, y, samplesize, 10.0, work);
            if(maxbf < bf[j]){maxid = j; maxbf=bf[j];}
            ntested++;
            printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%d\t%d\t%lf\n", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]-1, loccat[j], loccat2[j], bf[j]);
        }
    }
    //printf("%s\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%lf\t%d\n", fid, chr, pos[maxid], rss[maxid], ba0[maxid], ba1[maxid], af[maxid], rsq[maxid], vt[maxid], bf[maxid], ntested);
}









int mainds(){
    double ds[4] = {0.9, 0.8, 0.2, 0.1};
    double ap[8];
    int gt[2] ={0,1};
    ds2ap(gt, ds, 4, ap);
    printV(ds, 4);
    printV(ap, 4);
    printV(ap+4, 4);
}

int main(int argc, char** argv){
    exp_gt_gtdsgl = 0;
#ifdef ASCOUNT
    fprintf(stderr, "\n\nAS count\n\n");
    countAS(argv[1], argv[2]);
    //createDataBase(argv[1], argv[2]);
#else
    fprintf(stderr, "\n\nRSQ\n\n");
    int ld2 = argc>5 ? atoi(argv[argc-1]) : 0;
    getRsq(argv[1], argv[2], argv[3], argv[4], ld2);
#endif
    //lm(argc, argv);
}

int main0(){
    double gl[15]={0.01,0.01,0,0.95,0.02,0,0,0,0,0,0,0,0,0,0};
    double* ds; ds = (double*)calloc(100, sizeof(double));
    double* ap; ap = (double*)calloc(100, sizeof(double));
    double* d;   d = (double*)calloc(100, sizeof(double));
    int gt[2];
    gt[0]=0; gt[1]=2;
    gl2ap(gt, gl, 5, ap, d);
    printV(ap, 5);
    printV(ap+5, 5);
}

int main1(){
    int j, k, j1, k1, nal=5;
    double* gl; gl = (double*)calloc(5000, sizeof(double));
    double* ds; ds = (double*)calloc(100, sizeof(double));
    int gt[2];
    for(k1=0; k1<nal; k1++){
        for(j1=0; j1<=k1; j1++){
            gt[0]=k1; gt[1]=j1;
            
            // GT
            fprintf(stdout, "%d|%d\t", gt[0], gt[1]);
            
            gt2dsgl(gt, nal, ds, gl);
            
            // DS
            for(k=1; k<nal-1; k++){
                fprintf(stdout, "%.0lf,", ds[k]);
            }
            fprintf(stdout, "%.0lf\t", ds[k]);
            
            // GL
            for(k=0; k<choose(nal,2)+nal-1; k++){
                fprintf(stdout, "%.0lf,", gl[k]);
            }
            fprintf(stdout, "%.0lf\n", gl[k]);
        }
    }
}









int mainCOmbdose(int argc, char** argv){
    int n = atoi(argv[1]);
    int v[n];
    double ds[4] = {0.1,0.8,0.7,0.4};
    int an = achoose(n);
    double bds[an];
    getCombk(n, 0, v, -1, achoose(n), ds, bds, 1);
}

















