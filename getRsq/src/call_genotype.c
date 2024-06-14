#include "loadVCF.h"
#include "util.h"
#include "getVCFHeader.h"
#include <gsl/gsl_sf_gamma.h>

int printHeader(char* fn){
    gzFile f = gzopen(fn, "rb6f");
    int nlines=0;
    char c;
    while((c=gzgetc(f)) != EOF){
        if(c=='\n'){nlines++;}
    	putchar(c);
    }
    gzclose(f);
    return nlines;
}

double lgamma(double x){
        return gsl_sf_lngamma(x);
}

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

void nk_scale(double* x, long N){
    double mx = 0.0;
    double vx = 0.0;
    int i;
    for(i=0; i<N; i++){
        mx += x[i];
        vx += x[i]*x[i];
    }
    mx /= (double)N;
    vx /= (double)N;
    vx -= mx*mx;
    double sx = sqrt(vx)*sqrt((double)N);
    for(i=0; i<N; i++){
        x[i] = (x[i]-mx)/sx;
    }
    return;
};

void nk_Cor(double* X, double* R, int n, int m){
    int i,j,k;
    for(i=0; i<m; i++){
        nk_scale(X+i*n, n);
    }
    for(i=0; i<m; i++){
        for(j=i; j<m; j++){
            if(i==j){
                if(isnan(X[j*n])==0)R[i+j*m] = 1.0;
            }else{
                R[i+j*m] = R[j+i*m] = 0.0;
                for(k=0; k<n; k++){R[i+j*m] += X[k+i*n]*X[k+j*n];}
                R[j+i*m] = R[i+j*m];
            }
        }
    }
    return;
}

double nk_cor(double* x, double* y, int n){
    return nk_var(x, y ,n) / sqrt(nk_var(x, x, n)) / sqrt(nk_var(y, y, n));
}

int bdfscanf1(char* filename, double** x, int n0){
    FILE* f; f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    int n = ftell(f)/8;
    fseek(f, 0, SEEK_SET);
    if(n0>0){n=n0;}
    (*x) = (double*)calloc(n, sizeof(double));
    n = fread(*x, sizeof(double), n, f);
    fclose(f);
    return n;
}

int bdfscanf(char* filename, double** x){
    return bdfscanf1(filename, x, -1);
}


int gzfdscanf(char* filename, double** x){// N x P matrix
    gzFile f = gzopen(filename, "rb6f");
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
    char format[100];
    
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
        }else if(strcmp(format,"PS")==0){
            formatID[i]=FORMAT_PS;
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
        }else if(startWith("DR2=",infostr)==0){
            sscanf(infostr, "DR2=%lf", &(vinfo->RSQ));
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
    for(i=0; i<n-1; i++){fprintf(stdout, "%lf ", x[i]);}
    fprintf(stdout, "%lf\n", x[i]);
}
void printV2(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf ", x[i]);}
    fprintf(stderr, "%lf\n", x[i]);
}
void printVint(int* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stdout, "%d ", x[i]);}
    fprintf(stdout, "%d\n", x[i]);
}
void printVlog(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf ", log(x[i]));}
    fprintf(stderr, "%lf\n", log(x[i]));
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


int getCombk(int n, int k, int* v, int id, int idmax, double* ds, double* bds, int samplesize, int* gt, int* bgt, double* gl, double* bgl){
    if(k==n){
        if(id==-1 || id>=idmax){return id+1;}
        int i;
        //printf("%d ", id); for(i=0; i<n; i++){printf("%d ", v[i]);} 
        for(i=0; i<n; i++){
            if(v[i]>0){
                bds[id*samplesize]+=ds[i];
            }
        }
        bgt[id*samplesize*2+0] = v[gt[0]];
        bgt[id*samplesize*2+1] = v[gt[1]];
        int j;
        for(i=0; i<n; i++){
            for(j=i; j<n; j++){
                bgl[id*samplesize*3+v[i]+v[j]] += gl[(j*(j+1)/2)+i];
            }
        }
        //printf("%lf\n", bds[id]);
        
        return id+1;
    }
    v[k]=0;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize, gt, bgt, gl, bgl);
    v[k]=1;
    id = getCombk(n, k+1, v, id, idmax, ds, bds, samplesize, gt, bgt, gl, bgl);
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

int isPhased(char* x, int l){
    int i;
    for(i=0; i<l; i++){
        //fprintf(stderr, "%c\n", x[i]);
        if(x[i]=='|'){
            return 1;
        }else if(x[i]=='/'){
            return 0;
        }
    }
}

int parseGT(char* cell, int* gt){
    if(cell[0]=='.' && cell[1]=='/' && cell[2]=='.'){gt[0]=gt[1]=-9999999; return 4;}
    if(cell[0]=='.'){gt[0]=gt[1]=-9999999; return 2;}
    int nchar1;
    if(isPhased(cell, strlen(cell))>0){
        sscanf(cell, "%d|%d%n", gt, gt+1, &nchar1);
    }else{
        sscanf(cell, "%d/%d%n", gt, gt+1, &nchar1);
    }
    return nchar1 + 1;
}

int parsePS(char* cell, int n){
    int nchar1, nchar=0, i;
    double ds;
    for(i=1; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", &ds, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", &ds, &nchar1);
    
    nchar += nchar1 + 1;
    return nchar;
}


int parseDS(char* cell, int n, double* ds){
    if(cell[0]=='.'){ds[0]=-9999999; return 2;}
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
        if(cell[i]=='\t' || cell[i]==':' || cell[i]=='\n' || cell[i]=='\0'){return i+1;}
    }
}

int parseGP(char* cell, int n, double* gl){
    if(cell[0]=='.'){gl[0]=gl[1]=gl[2]=-9999999; return 2;}
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
    if(cell[0]=='.'){gl[0]=gl[1]=gl[2]=-9999999; return 2;}
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
    if(cell[0]=='.'){ap[0]=ap[1]=-9999999; return 2;}
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
    if(cell[0]=='.'){as[0]=as[1]=-9999999; return 2;}
    int nchar1, nchar=0, i;
    for(i=0; i<n-1; i++){
        sscanf(cell+nchar, "%lf,%n", as+i, &nchar1);
        nchar += nchar1;
    }
    sscanf(cell+nchar, "%lf%n", as+i, &nchar1);
    
    nchar += nchar1 + 1;
    return nchar;
}

int parseCell(char* cell, int nal, int* gt, double* gl, double* ap, double* asc, double* ds, int* formatID, int nfields){
    int i, offs=0, k=0;
    for(k=0; k<nfields; k++){
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
    char* a0;     a0    =(char*)calloc(100000, sizeof(char));
    char* a1;     a1    =(char*)calloc(100000, sizeof(char));
    
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
    char* a0;     a0    =(char*)calloc(100000, sizeof(char));
    char* a1;     a1    =(char*)calloc(100000, sizeof(char));
    
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

int loadBed1(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val, int nrow){
    
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
    
    int isalloc=1;  // memory has been allocated
    if(nrow<0){
        isalloc=0;
        nrow = nrowBed(fname, reg);
        (*pchr)=(char*)calloc(1000, sizeof(char));
        (*ppos1)=(int*)calloc(nrow, sizeof(int));
        (*ppos2)=(int*)calloc(nrow, sizeof(int));
    }
    
    int gstart=270000000, gend=0;
    
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    i=0;
    int nchar, ncharval;
    int nval=0;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d%n", *pchr, (*ppos1)+i, (*ppos2)+i, &nchar);
        if(i==0){
            for(k=nchar; k<strlen(str.s); k++){
                if(str.s[k]=='\t'){nval++;}
            }
            if(isalloc==0)(*val) = (double*)calloc(nrow*nval, sizeof(double));
        }
        for(k=0; k<nval; k++){
            sscanf(str.s+nchar, "\t%lf%n", (*val)+i+k*nrow, &ncharval);
            nchar += ncharval;
        }
        i++;
    }
    tbx_itr_destroy(itr);
    return nrow;
}

int loadBed(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val){
    return loadBed1(fname, reg, pchr, ppos1, ppos2, val, -1);
}

double getCovFromBedG(const char* fname, char* chr, int pos){
    char* reg; reg = (char*)calloc(1000, sizeof(char));
    sprintf(reg, "%s:%d-%d", chr, pos, pos+1);
    //fprintf(stderr, "%s %s\n", fname, reg);
    int n = nrowBed(fname, reg);
    int i;
    char* chrs; chrs = (char*)calloc(n, sizeof(char));
    int* pos1; pos1 = (int*)calloc(n, sizeof(int));
    int* pos2; pos2 = (int*)calloc(n, sizeof(int));
    double* val; val = (double*)calloc(n, sizeof(double));
    loadBed1(fname, reg, &chrs, &pos1, &pos2, &val, n);
    double res = val[0];
    free(reg);
    free(chrs);
    free(pos1);
    free(pos2);
    free(val);
    return res;
}

double cov2eta(double x, double beta0, double* beta1, double* xk, int nk, int EXPIT){
    double res = beta0 + beta1[0]*x;
    int i;
    //fprintf(stderr, "%lf %lf ", x, beta0);
    for(i=0; i<nk; i++){
        //fprintf(stderr, "%lf ", beta1[i+1]);
        res += beta1[i+1] * rk1(x, xk[i]);
    }
    //fprintf(stderr, "\n");
    return EXPIT>0 ? expit(res) : res;
}

int nrowVCF(const char* fname, const char* reg, int* pbinarize, int* pnvars, int* psamplesize){
    
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
    char* a0;     a0    =(char*)calloc(100000, sizeof(char));
    char* a1;     a1    =(char*)calloc(100000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info  =(char*)calloc(1000, sizeof(char));
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nbivars=0;
    int nbivars1=0;
    int nal;
    int maxnbivars=0;
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", 
               chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        
        body = str.s+nchar;
        
        if(nvars==0){nsamples = countFields(body, '\t');}
        
        nvars++;
        nal = countFields(a1, ',')+1;
        
        nbivars1 = achoose(nal);
        nbivars += nbivars1;
        if(maxnbivars<nbivars1){ maxnbivars=nbivars1; }
        
    }
    tbx_itr_destroy(itr);
    (*pnvars) = nvars;
    (*psamplesize) = nsamples;
    if((*pbinarize)>0){
        pbinarize[0] = maxnbivars; 
        return nbivars;
    }else{
        return nvars;
    }
}


int call_genotype(const char* fname, const char* reg, double** pth, int* outform, double mincov, double callrate0){
    
    int i, k, printOrig=0, j, jform;
    
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
    char* chr; chr = (char*)calloc(1000, sizeof(char));;
    int pos;
    char* rs;     rs    =(char*)calloc(1000, sizeof(char));
    char* a0;     a0    =(char*)calloc(100000, sizeof(char));
    char* a1;     a1    =(char*)calloc(100000, sizeof(char));
    
    char* qual;   qual  =(char*)calloc(1000, sizeof(char));
    char* filter; filter=(char*)calloc(1000, sizeof(char));
    char* info;   info=(char*)calloc(1000, sizeof(char));
    VCF_info vinfo;
    char* format; format=(char*)calloc(1000, sizeof(char));
    
    double* gl; gl = (double*)calloc(210, sizeof(double));
    double* ds; ds = (double*)calloc(20, sizeof(double));
    double* af; af = (double*)calloc(20, sizeof(double));
    double* asc; asc = (double*)calloc(20, sizeof(double));
    double* ap; ap = (double*)calloc(40, sizeof(double));
    int*    gt; gt = (int*)calloc(2, sizeof(int));
    double* d;  d  = (double*)calloc(400, sizeof(double));
    
    int* formatID; formatID=(int*)calloc(100, sizeof(int));
    
    char* body;
    
    int nsamples;
    int nchar;  // n of characters of body
    int nvars=0;
    int nfields;
    
    int maxnbivars = 1;
    
    int* work; work=(int*)calloc(100, sizeof(int));
    char* cwork; cwork=(char*)calloc(1000, sizeof(char));
    int offs;
    char formatlab[5][3] = {"GT", "GL", "GP", "AS", "DS"};
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        body = str.s+nchar;
        nfields = parseFormat(format, formatID);
        parseInfo(info, cwork, &vinfo);
	offs = 0;
	nsamples = countFields(body, '\t');
	// fprintf(stderr, "nsamples=%d\n", nsamples);
	double callrate = 0.0;
        for(i=0; i<nsamples; i++){ offs += parseCell(body+offs, 2, gt, gl, ap, asc, ds, formatID, nfields); if(asc[0]+asc[1]>mincov){callrate++;}}
        offs = 0;
	if(callrate>callrate0*(double)nsamples && vinfo.AF>0.0){
            printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t", chr, pos, rs, a0, a1, qual, filter, info);
            jform=0;
            for(j=0; j<5; j++){
                if(outform[j]==1){if(jform==0){printf("%s",formatlab[j]);}else{printf(":%s",formatlab[j]);}; jform++;}
            }
            for(i=0; i<nsamples; i++){
                offs += parseCell(body+offs, 2, gt, gl, ap, asc, ds, formatID, nfields);
                if(asc[0]+asc[1]<0.5){
                    jform = 0; for(j=0; j<5; j++){if(outform[j]==1){if(jform==0){printf("\t.");}else{printf(":.");}; jform++;}}
                }else{
                    double th = (*pth)[i];
                    double a=0.99*th; 
                    double b=0.01*th; 
                    double n = asc[0]+asc[1];
                    double k = asc[0];
                    int *gt2; gt2 = (int*)calloc(3, sizeof(int));
                    double *gl3; gl3 = (double*)calloc(3, sizeof(double));
                    double *gp3; gp3 = (double*)calloc(3, sizeof(double));
                    gl3[0] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + lgamma(k+a) + lgamma(n-k+b) - lgamma(n+a+b) + lgamma(a+b) - lgamma(a) - lgamma(b);
                    a = b = 0.5*th;
                    gl3[1] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + lgamma(k+a) + lgamma(n-k+b) - lgamma(n+a+b) + lgamma(a+b) - lgamma(a) - lgamma(b);
                    b=0.99*th; a=0.01*th;
        	    gl3[2] = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + lgamma(k+a) + lgamma(n-k+b) - lgamma(n+a+b) + lgamma(a+b) - lgamma(a) - lgamma(b);
                    double geta = max(gl3[0], max(gl3[1], gl3[2]));
                    double scalingfactor = log10(exp(gl3[0]-geta)+exp(gl3[1]-geta)+exp(gl3[2]-geta));
                    gl3[0] = (gl3[0]-geta)/log(10.) - scalingfactor;
                    gl3[1] = (gl3[1]-geta)/log(10.) - scalingfactor;
                    gl3[2] = (gl3[2]-geta)/log(10.) - scalingfactor;
    		    //printf("\t%lf,%lf,%lf:%lf,%lf", gl[0], gl[1], gl[2],asc[0],asc[1]);
    		    //printf("\t%lf,%lf,%lf", gl3[0], gl3[1], gl3[2]);
                    gp3[0] = exp(gl3[0])/(exp(gl3[0])+exp(gl3[1])+exp(gl3[2]));
                    gp3[1] = exp(gl3[1])/(exp(gl3[0])+exp(gl3[1])+exp(gl3[2]));
                    gp3[2] = exp(gl3[2])/(exp(gl3[0])+exp(gl3[1])+exp(gl3[2]));
                    if(gp3[2]>gp3[1] && gp3[2]>gp3[0]){gt2[0]=gt2[1]=1;}else if(gp3[1]>gp3[2] && gp3[1]>gp3[0]){gt2[0]=0;gt2[1]=1;}else{gt2[0]=gt2[1]=0;}
                    jform = 0;
                    if(outform[0]==1){if(jform==0){printf("\t%d/%d",gt2[0],gt2[1]);}else{printf(":%d/%d",gt2[0],gt2[1]);}; jform++;}
                    if(outform[1]==1){if(jform==0){printf("\t%lf,%lf,%lf",gl3[0], gl3[1], gl3[2]);}else{printf(":%lf,%lf,%lf",gl3[0], gl3[1], gl3[2]);}; jform++;}
                    if(outform[2]==1){if(jform==0){printf("\t%lf,%lf,%lf",gp3[0], gp3[1], gp3[2]);}else{printf(":%lf,%lf,%lf",gp3[0], gp3[1], gp3[2]);}; jform++;}
                    if(outform[3]==1){if(jform==0){printf("\t%.0lf,%.0lf",asc[0],asc[1]);}else{printf(":%.0lf,%.0lf",asc[0],asc[1]);}; jform++;}
                    if(outform[4]==1){if(jform==0){printf("\t%lf",gp3[1]+2.*gp3[2]);}else{printf(":%lf",gp3[1]+2.*gp3[2]);}; jform++;}
                }
            }
	    printf("\n");
        }
    }
    tbx_itr_destroy(itr);
}


int main(int argc, char** argv){
    int header = 0;
    int i;
    for(i=0; i<argc; i++){
        if(strcmp(argv[i], "-h")==0){header=1;}
        if(strcmp(argv[i], "-H")==0){header=2;}
    }
    int* outform; outform = (int*)calloc(5, sizeof(int));
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-gt")==0){outform[0]=1;}}
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-gl")==0){outform[1]=1;}}
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-gp")==0){outform[2]=1;}}
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-as")==0){outform[3]=1;}}
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-ds")==0){outform[4]=1;}}
    double* ths;
    double mincov = 10.0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--min-cov")==0){mincov=(double)atof(argv[i+1]);break;}}
    double callrate0 = 0.9;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--call-rate")==0){callrate0=(double)atof(argv[i+1]);break;}}
    int nth = gzfdscanf(argv[3], &ths);
    exp_gt_gtdsgl = 0;
    printHeader(argv[4]);
    call_genotype(argv[1], argv[2], &ths, outform, mincov, callrate0);
}








