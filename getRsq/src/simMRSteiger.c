#include "loadVCF.h"
#include "util.h"




double fisherZ(double r){
    return log( (1.+r)/(1.-r) )/2.0;
}
double steigerZ(double rgx, double rgy, double rxy, int N){
    
    double Zgx = fisherZ(rgx);
    double Zgy = fisherZ(rgy);
    
    double rm2 = (rgx*rgx + rgy*rgy)/2.0;
    
    double f = (1.-rxy)/(1.-rm2)/2.0;
    
    double h = (1.-f*rm2)/(1.-rm2);
    
    return (Zgx-Zgy)*sqrt(((double)N-3.)/(1.-rxy)/2./h);
}
double getMRTZ(double* g, double* x, double* y, int N, double* pZ){
    int i, j, k;
    double dn = (double)N;
    double g1, g2, x1, x2, y1, y2, gx, gy, xy;
    double s2 = 10.;
    g1=g2=x1=x2=y1=y2=gy=xy=gx=0.0;
    for(i=0; i<N; i++){
        x1 += x[i];
        x2 += x[i]*x[i];
        y1 += y[i];
        y2 += y[i]*y[i];
        g1 += g[i];
        g2 += g[i]*g[i];
        
        gx += x[i]*g[i];
        xy += x[i]*y[i];
        gy += g[i]*y[i];
    }
    x1 /= dn;
    x2 /= dn;
    y1 /= dn;
    y2 /= dn;
    g1 /= dn;
    g2 /= dn;
    gx /= dn;
    gy /= dn;
    xy /= dn;
    x2 -= x1*x1;
    y2 -= y1*y1;
    g2 -= g1*g1;
    
    gx -= g1*x1;
    gy -= g1*y1;
    xy -= x1*y1;
    
    gx /= sqrt(g2*x2);
    gy /= sqrt(g2*y2);
    xy /= sqrt(x2*y2);
    if(gx<0.0 && gy<0.0){
        pZ[0] = steigerZ(-gx, -gy,  xy, N);
    }else if(gx<0.0 && gy>0.0){
        pZ[0] = steigerZ(-gx,  gy, -xy, N);
    }else if(gx>0.0 && gy<0.0){
        pZ[0] = steigerZ( gx, -gy, -xy, N);
    }else{
        pZ[0] = steigerZ(gx, gy, xy, N);
    }
    return sqrt((dn-2.)*gy*gy/(1.-xy*xy+(gy/gx-xy)*(gy/gx-xy)));
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
    for(i=0; i<n-1; i++){fprintf(stdout, "%lf,", x[i]);}
    fprintf(stdout, "%lf\n", x[i]);
}
void printV2(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf,", x[i]);}
    fprintf(stderr, "%lf\n", x[i]);
}
void printVint(int* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stdout, "%d,", x[i]);}
    fprintf(stdout, "%d\n", x[i]);
}
void printVlog(double* x, int n){
    int i;
    for(i=0; i<n-1; i++){fprintf(stderr, "%lf,", log(x[i]));}
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

int loadUKBB(const char* fname, const char* reg, char** pchr, int** ppos1, int** ppos2, double** val){
    
    int i, k;
    int nrow;
    
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
    
    
    nrow = nrowBed(fname, reg);
    (*pchr)=(char*)calloc(1000, sizeof(char));
    (*ppos1)=(int*)calloc(nrow, sizeof(int));
    (*ppos2)=(int*)calloc(nrow, sizeof(int));
    (*val) = (double*)calloc(nrow, sizeof(double));
    char* al1; al1 = (char*)calloc(100, sizeof(char));
    char* al2; al2 = (char*)calloc(100, sizeof(char));
    char* rsid; rsid = (char*)calloc(100, sizeof(char));
    int nSmp;
    double ac, ytx, beta, se, tstat, pval;
    
    int gstart=270000000, gend=0;
    
    // chr, start, end, allele 1, allele 2, rsid, nCompleteSamples, AC, ytx, beta, se, tstat, pval
    
    hts_itr_t *itr = tbx_itr_querys(tbx, reg);
    i=0;
    int nchar, ncharval;
    while (tbx_itr_next(fp, tbx, itr, &str) >= 0){
        sscanf(str.s, "%[^\t]\t%d\t%d\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", *pchr, (*ppos1)+i, (*ppos2)+i, al1, al2, rsid, &nSmp, &ac, &ytx, &beta, &se, &tstat, &pval);
        (*val)[i] = getLogWABFfromBetaSE(beta, se);
        i++;
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
    
    int i, k, printOrig=0;
    
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
    int nbivars=0, nbivars1;
    int nal;
    int nfields;
    
    int ncol_ds1 = nrowVCF(fname, reg, 1, &nvars, &nsamples);
    //fprintf(stderr, "getBiVCF %d\n", ncol_ds1);
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
        printOrig=0;
        //if(l%100==0)fprintf(stderr, "%d\n", l);
        sscanf(str.s, "%[^\t]\t%d\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%[^\t]\t%n", chr, &pos, rs, a0, a1, qual, filter, info, format, &nchar);
        body = str.s+nchar;
        nal = countFields(a1, ',')+1;
        
        nfields = parseFormat(format, formatID);
        parseInfo(info, cwork, &vinfo);
        
        nbivars1 = achoose(nal);
        //fprintf(stderr, "achoose=%ld\n", nbivars1);
        for(i=nbivars; i<nbivars+nbivars1; i++){
            (*pba0)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*pba1)[i] = (char*)calloc(strlen(a0)+strlen(a1), sizeof(char));
            (*prss)[i] = (char*)calloc(strlen(rs)+1, sizeof(char));
            strcpy((*prss)[i], rs);
            (*ppos)[i] = pos;
            (*vtype)[i] = vinfo.VT;
        }
        if(nal>2){
            getCombAlk(nal, 0, work, -1, nbivars1, a0, a1, (*pba0)+nbivars, (*pba1)+nbivars);
        }else{
            strcpy((*pba0)[nbivars], a0);
            strcpy((*pba1)[nbivars], a1);
        }
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
        //fprintf(stderr, "%s\n", body);
        for(i=0; i<nsamples; i++){
            //fprintf(stderr, "%s\n", body+offs);
            //int j; for(j=0; j<20; j++){fprintf(stderr, "%c", body[offs+j]);}fprintf(stderr, "\n");
            offs += parseCell(body+offs, nal, gt, gl, ap, asc, ds, formatID, nfields);
            //clearAs0(ds, nal);
            //ds[gt[0]]++;
            //ds[gt[1]]++;
            if(doseFormatExist(formatID, nfields, FORMAT_DS)==0){
                if(gt[0]<0){
                    ds[1]=-999999;
                }else{
                    if(doseFormatExist(formatID, nfields, FORMAT_GP)==1){
                        gl2ap(gt, gl, nal, ap, d);;
                    }else{
                        clearAs0(ap, nal*2);
                        ap[gt[0]]++;
                        ap[gt[1]+nal]++;
                    }
                    for(k=0; k<nal; k++){
                        ds[k] = ap[k]+ap[k+nal];
                    }
                }
            }
            if(nal>2){
                getCombk(nal, 0, work, -1, nbivars1, ds, (*pds1)+nbivars*nsamples+i, nsamples);
            }else{
                //getCombk(nal, 0, work, -1, nbivars1, ds, (*pds1)+nbivars*nsamples+i, nsamples);
                (*pds1)[nbivars*nsamples+i] = ds[1];
            }
        }
#ifdef PRINTVCF
        // allele freq and print
        for(i=0; i<nbivars1; i++){
            af[i] =  nk_dsum((*pds1)+(nbivars+i)*nsamples, nsamples, 1)/2.0/(double)nsamples;
            printOrig += af[i]>0.05 && af[i]<0.95 ? 1 : 0;
        }
        if(printOrig>0){printf("%s\n", str.s);}
#endif
        
        nbivars += nbivars1;
        l++;
    }
    //fprintf(stderr, "finish\n");
    tbx_itr_destroy(itr);
    
    
    for(i=0; i<nbivars; i++){
        //fprintf(stderr, "%s %s %s\n", rss[i], ba0[i], ba1[i]);
    }
    
}

int printVCF(const char* fname, const char* reg){ 
    double* ds;
    int samplesize;
    int nbivars;
    char* chr;
    int* pos;
    char** rss; 
    char** ba0; 
    char** ba1;
    int* vtype;
    int i, j;
    char sep;
    getBiVCF(fname, reg, &ds, &samplesize, &nbivars, &chr, &pos, &rss, &ba0, &ba1, &vtype);
    return 0;
}

int getDose(const char* fname, const char* reg){ 
    double* ds;
    int samplesize;
    int nbivars;
    char* chr;
    int* pos;
    char** rss; 
    char** ba0; 
    char** ba1;
    int* vtype;
    int i, j;
    char sep;
    getBiVCF(fname, reg, &ds, &samplesize, &nbivars, &chr, &pos, &rss, &ba0, &ba1, &vtype);
    for(j=0; j<nbivars; j++){
        double af = nk_dsum(ds+j*samplesize, samplesize, 1)/2.0/(double)samplesize;
        //if(af>0.05&&af<0.95){
            printf("%s:%s:%s\t", rss[j], ba0[j], ba1[j]);
            for(i=0; i<samplesize; i++){
                if(i==samplesize-1){sep='\n';}else{sep='\t';}
                printf("%lf%c", ds[j*samplesize+i], sep);
            }
        //}
    }
    return 0;
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
        fprintf(stderr, "2 variants LD\n");
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


// 0 .. npeaks in [-2W, 2W] for prior prob
// 0 .. M in [-W, W] for testing

int lm(int argc, char** argv){
    init_gsl_rand();
    double xk[4] = {0.1470588, 0.3823529, 0.6176471, 0.8529412};// for spline
    
    int i, j, k;
    
    int verbose=0;
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-v")==0){verbose=1; break;}}
    
     
    const char* fname = NULL;
    const char* fname2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--vcf")==0){fname = fname2 = argv[i+1]; break;}}
    if(fname==NULL){fprintf(stderr, "VCF file is missing.\n"); return 1;}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--vcf2")==0){fname2  = argv[i+1]; break;}}
    
    if(verbose>0){fprintf(stderr, "VCF : %s %s\n", fname, fname2);};
    
    int tss = 0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-center")==0){tss = atoi(argv[i+1]); break;}}
    if(tss==0){fprintf(stderr, "Feature center is missing.\n"); return 1;}
    
    char* chrom=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--chrom")==0){chrom = argv[i+1]; break;}}
    if(chrom==NULL){fprintf(stderr, "Chromosome is missing.\n"); return 1;}
    
    const char* fnamey  = NULL;
    const char* fnamey2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--fpkm")==0){fnamey = fnamey2 = argv[i+1]; break;}}
    if(fnamey==NULL){fprintf(stderr, "FPKM file is missing.\n"); return 1;}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--fpkm2")==0){fnamey2 = argv[i+1]; break;}}
    
    if(verbose>0){fprintf(stderr, "FPKM : %s %s\n", fnamey, fnamey2);};
    
    int fid = 0;      // start point of fpkm, starting from 1
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-id")==0){fid = atoi(argv[i+1]); break;}}
    if(fid==0){fprintf(stderr, "Feature ID is missing.\n"); return 1;}
    int fid2 = 0;      // pair id, starting from 1
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--feature-id2")==0){fid2 = atoi(argv[i+1]); break;}}
    
    if(verbose>0){fprintf(stderr, "FID : %d %d\n", fid, fid2);};
    
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--pairwise")==0){fprintf(stderr, "OLD PROGRAM!\n"); return 1;}}
    
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--bf2")==0){fid2=0;}}
    
    int atac=1;
    for(i=0; i<argc; i++){if(strcmp(argv[i], "--expression")==0){atac=0; break;}}
    double phi0=-1.0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--phi0")==0){phi0=(double)atof(argv[i+1]); break;}}
    double del0=-1.0;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--del0")==0){del0=(double)atof(argv[i+1]); break;}}
    
    
    int printBF=1;
    gzFile outBFf; 
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--outputBF")==0){printBF=2; outBFf = gzopen(argv[i+1], "ab6f");}}
    
    double* beta = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--coefficients")==0){
        if(verbose>0){fprintf(stderr, "Beta1 : %s\n", argv[i+1]);};
        //int nbeta=gzfdscanf(argv[i+1], &beta);
        int nbeta=bdfscanf1(argv[i+1], &beta, 10);
        if(verbose>0)printV(beta, nbeta); 
        printBF=0;
        break;}
    }
    double* beta2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--coefficients2")==0){
        if(verbose>0){fprintf(stderr, "Beta2 : %s\n", argv[i+1]);};  
        //int nbeta2=gzfdscanf(argv[i+1], &beta2); 
        int nbeta2=bdfscanf1(argv[i+1], &beta2, 10); 
        if(verbose>0)printV(beta2, nbeta2);
        break;}
    }
    if(beta2 == NULL){beta2 = beta; }
    
    double* gamma = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--gamma")==0){
        if(verbose>0){fprintf(stderr, "Gamma1 : %s\n", argv[i+1]);}; 
        int nbeta=bdfscanf1(argv[i+1], &gamma, 6);
        if(verbose>0)printV(gamma, nbeta); 
        break;}
    }
    //if(gamma==NULL){fprintf(stderr, "Gamma is null!\n"); return 0;}
    double* gamma2 = NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--gamma2")==0){
        if(verbose>0){fprintf(stderr, "Gamma2 : %s\n", argv[i+1]);}; 
        int nbeta=bdfscanf1(argv[i+1], &gamma2, 6);
        if(verbose>0)printV(gamma2, nbeta); 
        break;}
    }
    if(gamma2==NULL){gamma2 = gamma;}
    
    // fit pairwise model to get posterior prob by solving the causality
    int print_posterior=0;
    double* beta_psi;
    double* Psi;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--psi")==0){
        Psi=(double*)calloc(4, sizeof(double));
        print_posterior=1;
        if(verbose>0){fprintf(stderr, "Psi1 : %s\n", argv[i+1]);};
        bdfscanf1(argv[i+1], &beta_psi, 12);
        if(verbose>0)printV(beta_psi, 12);
        break;}
    }
    
    
    for(i=0; i<argc; i++){if(strcmp(argv[i],"--force-print")==0){printBF=3;}}
    
    int randomperm=0;
    for(i=0; i<argc; i++){if(strcmp(argv[i], "-r")==0){randomperm=1; break;}}
    
    char* reg; reg = (char*)calloc(1000, sizeof(char));
    
    int wsize=500000;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--window-size")==0){wsize = atoi(argv[i+1])/2; break;}}
    
    int hypo=0;
    for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--hypo")==0){hypo = atoi(argv[k+1]);}}
    int bothsides = strcmp(fname, fname2)==0 ? 0 : 1;
    //bothsides=1; 
    if(print_posterior>0){ // fit pairwise model to get posterior prob by solving the consality for both sides
        // temporary comment in
        if(hypo==3 || hypo==5 || hypo==30 )bothsides=1;
    }
    
    int regstart = (fid2>0 && bothsides==1) ? tss-2*wsize : tss-wsize;
    if(regstart<=0){regstart = 1;}
    if(tss<0){tss=1;}
    
    sprintf(reg, "%s:%d-%d", chrom, regstart, (fid2==0) ? tss+wsize : tss+2*wsize);
    
    //fprintf(stderr, "reg=%s\n", reg);
    
    char* peakbed=NULL;
    int npeaks=0;
    int* pos1bed;
    int* pos2bed;
    int* pcents;
    double* midp;
    int fid_bed=-1;
    int fid2_bed=0;
    int fid_bed_end=0;// end + 1
    int fid_bed_sta=0;
    double* ph; // relative peak height
    int M=0; // num of peaks for pairewise tests tss < peaks < tss+wsize are tested.
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--peak-bed")==0){
        peakbed=argv[i+1]; 
        char* chrbed;
        npeaks = loadBed(peakbed, reg, &chrbed, &pos1bed, &pos2bed, &ph);
        pcents = (int*)calloc(npeaks, sizeof(int));
        midp   = (double*)calloc(npeaks, sizeof(double));
        for(j=0; j<npeaks; j++){
            pcents[j] = (pos1bed[j]+pos2bed[j])/2;
            midp[j] = ((double)pos1bed[j]+(double)pos2bed[j])/2.0;
            if(pos2bed[j]<tss-wsize){fid_bed_sta=j+1;}
            if(pos1bed[j]<=tss && tss<=pos2bed[j]){// target peak
                fid_bed = j;
            }else if(pos2bed[j]<tss){
                fid_bed = j;
            }
            if(pos1bed[j]<tss+wsize){fid_bed_end=j+1;}
        }
        if(fid2>0){
            if(bothsides==0){
                M = fid_bed_end - fid_bed - 1;
                fid2_bed=fid_bed+1;
            }else{
                M = fid_bed_end - fid_bed_sta;
                fid2 -= (fid_bed-fid_bed_sta) + 1;
                //fid_bed = fid_bed_sta;
                fid2_bed = fid_bed_sta;
                //fprintf(stderr, "fid_bed_sta=%d fid_bed_end=%d fid_bed=%d fid2_bed=%d fid2=%d", fid_bed_sta, fid_bed_end, fid_bed, fid2_bed, fid2);
            }
        }
        break;
    }}
    if(npeaks==0){fprintf(stderr, "no peak found\n"); return 0;}
    if(verbose>0){fprintf(stderr, "%d annotation peaks found.\n", npeaks);}
    
    double* ds;
    int nbivars, samplesize;
    char* chr;
    int* pos;
    char** ba0;
    char** ba1;
    char** rss;
    int* vt;
    
    double* work; work=(double*)calloc(9, sizeof(double));
    
    if(verbose>0){fprintf(stderr, "Loading genotype dose from %s in %s.\n", fname, reg);}
    getBiVCF(fname, reg, &ds, &samplesize, &nbivars, &chr, &pos, &rss, &ba0, &ba1, &vt);
    if(verbose>0){fprintf(stderr, "Genotype dose has been successfully loaded from %s in %s (%d x %d).\n", fname, reg, nbivars, samplesize);}
    
    double maf0 = 0.05;
    
    double* af;     af     = (double*)calloc(nbivars, sizeof(double));
    double* rsq;    rsq    = (double*)calloc(nbivars, sizeof(double));
    double* bf;     bf     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* bfmr;   bfmr   = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* bfmr2;  bfmr2  = (double*)calloc(nbivars * (M+1), sizeof(double));
    
    // posterior prob from pairwise model
    double* Zj;  Zj  = (double*)calloc(nbivars*2, sizeof(double));
    //double* Zk;  Zk  = (double*)calloc(nbivars * (M+1), sizeof(double));
    
    double* be;     be     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* se;     se     = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* eta;    eta    = (double*)calloc(nbivars * (M+1), sizeof(double));
    double* eta0;   eta0   = (double*)calloc(nbivars, sizeof(double));
    double* covs;   covs   = (double*)calloc(nbivars*2, sizeof(double));
    double* Pi1;    Pi1    = (double*)calloc(M+1, sizeof(double));
    if(gamma!=NULL){
        for(i=0; i<argc-1; i++){
            if(strcmp(argv[i], "--feature-coverage-quantile")==0){ Pi1[0] = cov2eta((double)atof(argv[i+1]), gamma[0], gamma+1, xk, 4, 1); }
        }
        //Pi1[0] = cov2eta(ph[fid_bed+npeaks], gamma[0], gamma+1, xk, 4, 1);
        if(verbose>0) fprintf(stderr, "Pi1=%lf\n", Pi1[0]);
        for(i=0;i<M;i++){ Pi1[i+1] = cov2eta(ph[fid2_bed+i+npeaks], gamma2[0], gamma2+1, xk, 4, 1); }
    }
    int*    loccat; loccat = (int*)calloc(nbivars * (M+1), sizeof(int));
    int*    loccatid; loccatid = (int*)calloc(nbivars, sizeof(int)); for(i=0;i<nbivars;i++){loccatid[i]=-1;}
    double* w;      w      = (double*)calloc(nbivars, sizeof(double)); // 1.0: effective loci; 0.0: unused (MAF==0 etc.)
    
    // vcf2
    double* ds2;
    int nbivars2, samplesize2;
    char* chr2;
    int* pos2;
    char** ba02;
    char** ba12;
    char** rss2;
    int* vt2;
    
    double* af2;
    double* rsq2;
    
    if(strcmp(fname, fname2)!=0){
        if(verbose>0){fprintf(stderr, "Loading genotype dose from %s in %s.\n", fname2, reg);}
        getBiVCF(fname2, reg, &ds2, &samplesize2, &nbivars2, &chr2, &pos2, &rss2, &ba02, &ba12, &vt2);
        if(verbose>0){fprintf(stderr, "Genotype dose has been successfully loaded from %s in %s (%d x %d).\n", fname2, reg, nbivars2, samplesize2);}
        
        af2     = (double*)calloc(nbivars, sizeof(double));
        rsq2    = (double*)calloc(nbivars, sizeof(double));
        if(nbivars != nbivars2){fprintf(stderr, "Different VCF conposition!\n"); return 1;}
    }else{
        ds2 = ds;
        samplesize2 = samplesize;
        nbivars2 = nbivars;
        chr2 = chr;
        pos2 = pos;
        ba02 = ba0;
        ba12 = ba1;
        rss2 = rss;
        vt2  = vt;
        af2  = af;
        rsq2 = rsq;
    }
    if(verbose>0){fprintf(stderr, "\n");}
    
    // outcome
    FILE* fy;
    fy = fopen(fnamey, "rb");
    double* y;  y  = (double*)calloc(samplesize, sizeof(double));
    double* y2;
    fseek(fy, samplesize*(fid-1)*sizeof(double), SEEK_SET);
    int fread_info = fread(y, sizeof(double), samplesize, fy);
    if(verbose>0){fprintf(stderr, "%d feature fpkms are loaded from %s.\n", samplesize, fnamey);}
    if(fid2>0){
        FILE* fy2;
        fy2 = fopen(fnamey2, "rb");
        y2 = (double*)calloc(samplesize2 * M, sizeof(double));
        fseek(fy2, samplesize2*(fid2-1)*sizeof(double), SEEK_SET);
        fread_info = fread(y2, sizeof(double), samplesize2 * M, fy2);
        if(verbose>0){fprintf(stderr, "%d x %d feature fpkms are loaded from %s.\n", M, samplesize2, fnamey2);}
    }
    if(verbose>0){fprintf(stderr, "\n");}
    
    srand((unsigned)(time(NULL)+getpid()));
    if(randomperm>0){
        randomise(y, samplesize);
    }
    
    int l, maxid;
    double maxbf = -1.0e10;
    int ntested=0;
    char* fcoverage=NULL; for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--coverage")==0){fcoverage = argv[i+1]; if(verbose>0){fprintf(stderr, "Coverage file: %s\n", fcoverage);}; break;}}
    //printBF=3;
    if(verbose>0) fprintf(stderr, "Prior calculation...");
    for(j=0; j<nbivars; j++){
        af[j]  = nk_dsum(ds+j*samplesize, samplesize, 1)/2.0/(double)samplesize;
        rsq[j] = nk_var(ds+j*samplesize, ds+j*samplesize, samplesize)/af[j]/(1.0-af[j])/2.0;
        if(beta!=NULL){
            eta[j] = eta0[j] = beta[vt[j]]; 
            for(l=0; l<M; l++){ eta[j+(l+1)*nbivars] = beta2[vt[j]]; }
        }
        //fprintf(stderr, "%lf %lf\n", af[j], rsq[j]);
        if(af[j]>maf0 && af[j]<1.0-maf0 && rsq[j]>0.3){
            //fprintf(stderr, ".");
            //Prior calculation
            for(k=0; k<npeaks; k++){// annnotation
                if(pos1bed[k]<=pos[j] && pos[j]<=pos2bed[k]){// in the kth peak
                    //fprintf(stderr,"+");
                    //double covAtj;
                    if(fcoverage!=NULL){ covs[j] = getCovFromBedG(fcoverage, chrom, pos[j]); } //fprintf(stderr, "%lf\n", covs[j]); }
                    if(beta!=NULL) eta0[j] += cov2eta(covs[j]/ph[k], beta[4], beta+5, xk, 4, 0);
                    if(k==fid_bed){// is target peak?
                        loccat[j] = 1; loccatid[j]=k;
                        if(beta!=NULL) eta[j] += cov2eta(covs[j]/ph[k], beta[3], beta+5, xk, 4, 0);
                    }else{
                        loccat[j] = 2; loccatid[j]=k;
                        if(beta!=NULL) eta[j] += cov2eta(covs[j]/ph[k], beta[4], beta+5, xk, 4, 0);
                    }
                    for(l=0; l<M; l++){// pairwise
                        if(k==l+fid2_bed){
                            loccat[j+(l+1)*nbivars] = 1;
                            if(beta2!=NULL) eta[j+(l+1)*nbivars] += cov2eta(covs[j]/ph[k], beta2[3], beta2+5, xk, 4, 0);
                        }else{
                            loccat[j+(l+1)*nbivars] = 2;
                            if(beta2!=NULL) eta[j+(l+1)*nbivars] += cov2eta(covs[j]/ph[k], beta2[4], beta2+5, xk, 4, 0);
                        }
                    }
                    break;
                }
            }
            
            // BF calculation
            //fprintf(stderr, "%s\t%lf\n", rss[j],ds[j*samplesize]);
            w[j] = 1.0;
            //bf[j] = getLogBF(ds+j*samplesize, y, samplesize, sqrt(10.0), work);
            bf[j] = getLogWABF(ds+j*samplesize, y, samplesize);
            be[j] = getBetaSigma(ds+j*samplesize, y, samplesize, se+j);
            if(fid2>0){// for pairwise
                af2[j]  = nk_dsum(ds2+j*samplesize2,       samplesize2, 1)/2.0/(double)samplesize2;
                rsq2[j] = nk_var( ds2+j*samplesize2, ds2+j*samplesize2, samplesize2)/af2[j]/(1.0-af2[j])/2.0;
                if(af2[j]>maf0 && af2[j]<1.0-maf0 && rsq2[j]>0.3 && af[j]>maf0 && af[j]<1.0-maf0 && rsq[j]>0.3){
                    for(k=0; k<M; k++){
                        //bf[j+(k+1)*nbivars] = getLogBF(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, 10.0, work);
                        bf[j+(k+1)*nbivars]    = getLogWABF(  ds2+j*samplesize2,    y2+k*samplesize2, samplesize2);
                        if(ds==ds2 || samplesize==samplesize2){
                            bfmr[j+(k+1)*nbivars]  = getLogWABFMR(ds2+j*samplesize2, y, y2+k*samplesize2, samplesize2);
                            bfmr2[j+(k+1)*nbivars] = getLogWABFMR(ds2+j*samplesize2, y2+k*samplesize2, y, samplesize2);
                        }else{
                            // fprintf(stderr, "atac-eqtl");
                            // expression - atac : not implemented
                            bfmr[ j+(k+1)*nbivars] = getLogWABFMR2(ds+ j*samplesize,  y,                samplesize,  ds2+j*samplesize2, y2+k*samplesize2, samplesize2);
                            bfmr2[j+(k+1)*nbivars] = getLogWABFMR2(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, ds+j*samplesize,   y,                samplesize );
                        }
                        be[j+(k+1)*nbivars]    = getBetaSigma(ds2+j*samplesize2, y2+k*samplesize2, samplesize2, se+(j+(k+1)*nbivars));
                    }
                }else{
                    w[j] = 0.0;
                }
            }
            //if(strcmp("rs2409780",rss[j])==0){w[j]=0.0;}
            //if(loccatid[j]==106){fprintf(stderr, "%s %d %d %lf\n", rss[j], loccat[j], loccatid[j], w[j]);}
            if(maxbf < bf[j]){maxid = j; maxbf=bf[j];}
            ntested++;
            if(printBF>0){
                if(printBF==1){
                    printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]);
                    for(k=0; k<M; k++){
                        printf("%d\t",  loccat[j+k*nbivars]);
                        printf("%lf\t%lf\t%lf\t%lf\t", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars], eta[j+k*nbivars]);
                    }
                    printf("%d\t", loccat[j+k*nbivars]);
                    printf("%lf\t%lf\t%lf\t%lf\t%lf\n", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars], eta[j+k*nbivars], covs[j]);
                }else if(printBF==3){// printing prior log odds
                    printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]);
                    for(k=1; k<2; k++){
                        //printf("\t%lf", eta[j+k*nbivars]);
                        printf("\t%lf\t%lf\t%lf\t%lf", bf[j], bf[j+k*nbivars], bfmr[j+k*nbivars], bfmr2[j+k*nbivars]);
                    }
                    printf("\n");
                }else{
                    gzprintf(outBFf, "%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t", fid, chr, pos[j], rss[j], ba0[j], ba1[j], af[j], rsq[j], vt[j]);
                    for(k=0; k<M; k++){
                        gzprintf(outBFf, "%d\t",  loccat[j+k*nbivars]);
                        gzprintf(outBFf, "%lf\t%lf\t%lf\t", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars]);
                    }
                    gzprintf(outBFf, "%d\t", loccat[j+k*nbivars]);
                    gzprintf(outBFf, "%lf\t%lf\t%lf\n", be[j+k*nbivars], se[j+k*nbivars], bf[j+k*nbivars]);
                }
            }
        }
    }
    if(verbose>0) fprintf(stderr, "Done\n");
    if(printBF==2){gzclose(outBFf);}
    
    // pairwise hierahical model for gwas
    for(i=0; i<argc-1; i++){if(strcmp(argv[i], "--bf2")==0){
        FILE* fbf2list; fbf2list = fopen(argv[i+1], "r");
        char* fbf2; fbf2  = (char*)calloc(1000, sizeof(char));
        
        while(fscanf(fbf2list, "%s\n", fbf2)!=EOF){
            //fprintf(stderr, "%s\n", fbf2);
        
            int nrowbf2 = nrowBed(fbf2, reg);
            char* pfbf2; for(k=strlen(fbf2)-1; k>=0; k--){if(fbf2[k]=='/'){pfbf2 = fbf2+k+1; break;}}
                //if(nrowbf2==0){fprintf(stderr, "no gwas bfs\n"); return 0;}
            if(nrowbf2<10){
                fprintf(stderr, "No line in %s in %s\n", fbf2, reg);
            }else{
                char* chrbf2;    chrbf2  = (char*)calloc(1000, sizeof(char));
                int* posbf2;     posbf2  = (int*)calloc(nrowbf2, sizeof(int));
                double* bf2orig; bf2orig = (double*)calloc(nrowbf2, sizeof(double));
                double* bf2;     bf2     = (double*)calloc(nbivars, sizeof(double));
                int* pos2bf2; pos2bf2 = (int*)calloc(nrowbf2, sizeof(int));
                //loadBed(fbf2, reg, &chrbf2, &posbf2, &pos2bf2, &bf2orig);
                loadUKBB(fbf2, reg, &chrbf2, &posbf2, &pos2bf2, &bf2orig);
                
                //printV(bf2orig, nrowbf2);
                clearAs0(w, nbivars);
                expand(pos2bf2, bf2orig, nrowbf2, pos, bf2, nbivars, w);
                //fprintf(stderr, "N ol vars=%lf\n", nk_dsum(w,nbivars,1)); 
                
                //double* pp3;  pp3  = (double*)calloc(3,  sizeof(double));
                double* pp12; pp12 = (double*)calloc(12, sizeof(double));
                //pwhm13(bf, bf2, bf, vt, loccat, w, nbivars, beta, pp3, pp12, &phi0, &del0);
                // Param 5
                //pwhmNewAtacGwas(bf, bf2, Pi1[0], w, nbivars, pp12);
                //printf("%d", fid);
                //for(j=0; j<6; j++){printf("\t%lf", log(pp12[j]));}printf("\n");
                
                
                gzFile outf=NULL;
                for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
                if(outf==NULL){
                    printf("%d\t%s\t%d", nrowbf2, pfbf2, fid);
                    for(k=0; k<6; k++){printf("\t%lf", log(pp12[k]));} printf("\n"); 
                }else{
                    gzprintf(outf, "%d\t%s\t%d", nrowbf2, pfbf2, fid);
                    for(k=0; k<6; k++){gzprintf(outf, "\t%lf", log(pp12[k]));} gzprintf(outf, "\n"); 
                    gzclose(outf);
                }
            }
        }
        return 0;
    }}
    
    
    int* cvs; cvs=(int*)calloc(6,sizeof(int));
    double* y_r;  y_r  = (double*)calloc(samplesize, sizeof(double));
    double* y2_r; y2_r = (double*)calloc(samplesize, sizeof(double));
    double* ya_r; ya_r = (double*)calloc(samplesize, sizeof(double));
    double* ds0;  ds0  = (double*)calloc(samplesize, sizeof(double));
    double* rsqOI      = (double*)calloc(4, sizeof(double));
    
    
    // pairewise hierahical model for atac & eqtl
    if(printBF==0 && fid2>0){
        if(verbose>0)fprintf(stderr, "Pairwise model\n");
        double* pp2;  pp2  = (double*)calloc(2,  sizeof(double));
        double* pp13; pp13 = (double*)calloc(12, sizeof(double));
        //double phi0 = 0.98451478; //1M caqtl   0.96687324;// 500Kb caqtl
        //double phi0 = 0.98883017;// eqtl-caqtl
        //double phi0 = 0.98478229; //eqtl
        for(j=0; j<M; j++){
            int nloci = 0; // #vars <- 500K - P1 - p2 - 500K ->
            int geta = 0;
            int cis_sta = mini(tss, pcents[fid2_bed+j]) - wsize; if(cis_sta<0){cis_sta=1;}
            int cis_end = maxi(tss, pcents[fid2_bed+j]) + wsize;
            for(k=0; k<nbivars; k++){
                if(pos[k]  <  cis_sta){geta++;}
                if(cis_sta <= pos[k] && pos[k] <= cis_end){nloci++;}
                if(cis_end <  pos[k]){break;}
            }
            /*if(fid2+j==84759){
                //printV(bf+geta, nloci);
                //printV(bf+(j+1)*nbivars+geta, nloci);
                printVint(loccat+geta, nloci);
                printVint(loccat+(j+1)*nbivars+geta, nloci);
            }*/
            
//fprintf(stderr, "%d\t%d\t%d\t%d\n", fid, fid2+j, pcents[j+fid2_bed]-tss, pos[k]-pos[geta]);
            //printf("%d\n", nloci);
            clearAs0(pp13, 13);
            
            double* steigerZ; steigerZ = calloc(8, sizeof(double));
            double* mrt;      mrt      = calloc(8, sizeof(double));
            
            if(strcmp(fname, fname2)==0){ // atac-atac
                if(fid!=fid2+j && fid2+j>fid){
                    
/*// pwhmfm
int hoge=260000000;
//fprintf(stderr, "%d\n", hoge);
if(fid2+j==fid+1){for(k=0; k<nloci; k++){if(w[geta+k]>0)printf("%d %d %lf %lf %lf %lf %lf %lf %s\n", pos[geta+k], loccat[geta+k], bf[geta+k], bf[(j+1)*nbivars+geta+k], bfmr[(j+1)*nbivars+geta+k], bfmr2[(j+1)*nbivars+geta+k], eta[geta+k], eta[(j+1)*nbivars+geta+k], rss[geta+k]);}}
*/

                    if(print_posterior>0 && atac>0){// posterior calculation after fitting
                        
                        
                        
                        rcausal1(eta0+geta, eta+geta, eta+(j+1)*nbivars+geta, w+geta, nloci, cvs);
                        
                        fprintf(stderr, "%s\n", rss[cvs[1]+geta]);
                        
                        // Linkage
                        rnorm(ds+(cvs[1]+geta)*samplesize, samplesize, 0.0, be[cvs[1]+geta],               se[cvs[1]+geta],               y_r);
                        rnorm(ds+(cvs[2]+geta)*samplesize, samplesize, 0.0, be[cvs[2]+geta+(j+1)*nbivars], se[cvs[2]+geta+(j+1)*nbivars], y2_r);
                        rsqOI[0] = nk_cor(ds+ (cvs[1]+geta)*samplesize, ds+ (cvs[2]+geta)*samplesize, samplesize);
                        
                        mrt[0] = getMRTZ(ds+(cvs[1]+geta)*samplesize, y_r, y2_r, samplesize, steigerZ+0);
                        //mrt[1] = getMRTZ(ds+(cvs[1]+geta)*samplesize, y2_r, y_r, samplesize, steigerZ+1);
                        
                        // Pleio
                        rnorm(ds+(cvs[0]+geta)*samplesize, samplesize, 0.0, be[cvs[0]              +geta], se[cvs[0]              +geta], y_r);
                        rnorm(ds+(cvs[0]+geta)*samplesize, samplesize, 0.0, be[cvs[0]+(j+1)*nbivars+geta], se[cvs[0]+(j+1)*nbivars+geta], y2_r);
                        
                        mrt[1] = getMRTZ(ds+(cvs[0]+geta)*samplesize, y_r, y2_r, samplesize, steigerZ+1);
                        //mrt[3] = getMRTZ(ds+(cvs[0]+geta)*samplesize, y2_r, y_r, samplesize, steigerZ+3);
                        
                        // Causal
                        double mrbe=be[cvs[1]+(j+1)*nbivars+geta]/be[cvs[1]+geta];
                        double mrse=getMRSE(y, y2+j*samplesize, samplesize, mrbe);
                        rnorm(ds+(cvs[1]+geta)*samplesize, samplesize, 0.0, be[cvs[1]+geta], se[cvs[1]+geta], y_r);
                        rnorm(y_r,                         samplesize, 0.0, mrbe,            mrse,            y2_r);
                        
                        mrt[2] = getMRTZ(ds+(cvs[1]+geta)*samplesize, y_r, y2_r, samplesize, steigerZ+2);
                        //mrt[5] = getMRTZ(ds+(cvs[1]+geta)*samplesize, y2_r, y_r, samplesize, steigerZ+5);
                        
                        // Causal2
                        mrbe=be[cvs[2]+geta]/be[cvs[2]+(j+1)*nbivars+geta];
                        mrse=getMRSE(y2+j*samplesize, y, samplesize, mrbe);
                        rnorm(ds+(cvs[2]+geta)*samplesize, samplesize, 0.0, be[cvs[2]+(j+1)*nbivars+geta], se[cvs[2]+(j+1)*nbivars+geta], y2_r);
                        rnorm(y2_r,                        samplesize, 0.0, mrbe,                          mrse,                          y_r);
                        
                        mrt[3] = getMRTZ(ds+(cvs[2]+geta)*samplesize, y_r, y2_r, samplesize, steigerZ+3);
                        //mrt[7] = getMRTZ(ds+(cvs[2]+geta)*samplesize, y2_r, y_r, samplesize, steigerZ+7);
                        
                        
                        gzFile outf=NULL;
                        for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
                        gzprintf(outf, "%d\t%d\t%lf", fid, fid2+j, rsqOI[0]);
                        
                        int jj=0;
                        for(jj=0; jj<4; jj++){
                            gzprintf(outf, "\t%lf\t%lf", mrt[jj], steigerZ[jj]);
                        }
                        gzprintf(outf, "\n");
                        gzclose(outf);
                        
                    }else{

                        if(atac>0){
                            //pwhm( bf+geta, bf+(j+1)*nbivars+geta, bfmr+(j+1)*nbivars+geta, bfmr2+(j+1)*nbivars+geta, vt+geta, loccat+geta, loccat+(j+1)*nbivars+geta, w+geta, nloci, beta, pp2, pp13, &phi0);
                            
                            // temporary comment out
                            pwhmnew(bf+geta, bf+(j+1)*nbivars+geta, bfmr+(j+1)*nbivars+geta, bfmr2+(j+1)*nbivars+geta, eta0+geta, eta+geta, eta+(j+1)*nbivars+geta, Pi1[0], Pi1[j+1], w+geta, nloci, pp2, pp13, loccatid+geta);
                            
                        }else{ // eqtl-eqtl
                            // not implemented
                            pwhm1(bf+geta, bf+(j+1)*nbivars+geta, vt+geta, loccat+geta, loccat+(j+1)*nbivars+geta, w+geta, nloci, beta, pp2, pp13, &phi0);
                        }
                        
                        
                        
                        // temporary comment out
                        gzFile outf=NULL;
                        for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
                        if(outf==NULL){
                            printf("%d\t%d", fid, fid2+j);
                            if(atac>0){ for(k=0; k<10; k++){printf("\t%lf", log(pp13[k]));} printf("\n"); }else{ for(k=0; k<5; k++){printf("\t%lf", pp13[k]);}printf("\n"); }
                        }else{
                            gzprintf(outf, "%d\t%d", fid, fid2+j);
                            if(atac>0){ for(k=0; k<10; k++){gzprintf(outf, "\t%lf", log(pp13[k]));} gzprintf(outf, "\n"); }else{ for(k=0; k<5; k++){printf("\t%lf", pp13[k]);}printf("\n"); }
                            gzclose(outf);
                        }
                    
                    }
                    
                }
            }else{ // atac-eqtl
                // 2 samples method
                //pwhmnewataceqtl(bf+geta, bf+(j+1)*nbivars+geta, eta0+geta, eta+(j+1)*nbivars+geta, Pi1[0], Pi1[j+1], w+geta, nloci, pp13, loccatid+geta, rss+geta, fid2+j);
                pwhmnewataceqtlAllParam(bf+geta, bf+(j+1)*nbivars+geta, eta0+geta, eta+(j+1)*nbivars+geta, Pi1[0], Pi1[j+1], w+geta, nloci, pp13, loccatid+geta, rss+geta, fid2+j);

                // joint obs.
                //pwhm12(bf+geta, bf+(j+1)*nbivars+geta, bfmr+(j+1)*nbivars+geta, bfmr2+(j+1)*nbivars+geta, vt+geta, loccat+geta, loccat+(j+1)*nbivars+geta, w+geta, nloci, beta, beta2, pp2, pp13, &phi0, rss+geta);
                //if(phi0<0.0){
                //    printf("%d\t%d\t", fid, fid2+j);
                //    for(k=0; k<2; k++){printf("%lf ", pp2[k]);}printf("\n");
                //}else{
                gzFile outf=NULL; 
                for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
                if(outf==NULL){
                    printf("%d\t%d", fid, fid2+j);
                    for(k=0; k<6; k++){printf("\t%lf", log(pp13[k]));}printf("\n");
                }else{
                    gzprintf(outf, "%d\t%d", fid, fid2+j);
                    for(k=0; k<6; k++){gzprintf(outf, "\t%lf", log(pp13[k]));}
                    gzprintf(outf, "\n");
                    gzclose(outf);
                }
                //}
            }
        }
        print_posterior=0;
        if(print_posterior>0){ // pwhmfm()
            int cis_sta = tss - wsize; if(cis_sta<0){cis_sta=1;}
            int cis_end = tss + wsize;
            double totzj=0.0, totzjnom=0.0;
            double* varloc; varloc=(double*)calloc(3,sizeof(double));
            for(i=0; i<nbivars; i++){
                if(cis_sta <= pos[i] && pos[i] <= cis_end){varloc[loccat[i]] += Zj[i]/totzj;}
            }
            
            // no pairwise solution
            pwhmfm0(bf, eta, Pi1[0], w, nbivars, Zj+nbivars);
            for(i=0; i<nbivars; i++){
                if(cis_sta <= pos[i] && pos[i] <= cis_end){
			totzj += Zj[i]; 
			totzjnom += Zj[i+nbivars];
		}
            }

            gzFile outf=NULL;
            for(k=0; k<argc-1; k++){if(strcmp(argv[k],"--output")==0){outf = gzopen(argv[k+1], "ab6f");}}
            if(outf==NULL){
                if(1==1){
                    for(k=0; k<nbivars; k++){
                        if(af[k]>maf0 && af[k]<1.0-maf0 && rsq[k]>0.3){
                            printf("%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n", fid, chr, pos[k], rss[k], ba0[k], ba1[k], af[k], rsq[k], vt[k], loccat[k], be[k], se[k], log(Zj[k]/totzj), log(Zj[k+nbivars]/totzjnom));
                        }
                    }
                }else{
                    printf("%d\t%lf\t%lf\t%lf\n", fid, varloc[0], varloc[1], varloc[2]);
                }
            }else{
                if(1==1){
                    for(k=0; k<nbivars; k++){
                        if(af[k]>maf0 && af[k]<1.0-maf0 && rsq[k]>0.3){
                            gzprintf(outf, "%d\t%s\t%d\t%s\t%s\t%s\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\n", fid, chr, pos[k], rss[k], ba0[k], ba1[k], af[k], rsq[k], vt[k], loccat[k], be[k], se[k], log(Zj[k]/totzj), log(Zj[k+nbivars]/totzjnom));
                        }
                    }
                    gzclose(outf);
                }else{
                    gzprintf(outf, "%d\t%lf\t%lf\t%lf\n", fid, varloc[0], varloc[1], varloc[2]);
                    //for(i=0; i<nbivars; i++){ if(cis_sta <= pos[i] && pos[i] <= cis_end){gzprintf(outf, "%lf\n", Zj[i]/totzj);} }
                    gzclose(outf);
                }
            }
            
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
#ifdef EXPANDVCF
    exp_gt_gtdsgl = 1;
    expandVCF(argv[1], argv[2]);
#elif ASCOUNT
    fprintf(stderr, "\n\nAS count\n\n");
    countAS(argv[1], argv[2]);
    //createDataBase(argv[1], argv[2]);
#elif GETRSQ
    int ld2 = argc>5 ? atoi(argv[5]) : 0;
    verbose_loadVCF = 0;
    if(verbose_loadVCF>0)fprintf(stderr, "\ngetRSQ\n\n");
    if(verbose_loadVCF>0)fprintf(stderr, " %s %s %s fid=%s 2SNP=%d\n\n", argv[1], argv[2], argv[3], argv[4], ld2);
    getRsq(argv[1], argv[2], argv[3], argv[4], ld2);
#elif BAYESLM
    if(verbose_loadVCF>0)fprintf(stderr, "\n\nbayesLm\n\n");
    lm(argc, argv);
#elif GETDOSE
    getDose(argv[1], argv[2]);
#elif PRINTVCF
    printVCF(argv[1], argv[2]);
#endif
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

















