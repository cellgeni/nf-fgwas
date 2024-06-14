#include "hierarchical.h"
#include "util.h"

double nk_mean(double* x, int n, int ldx){
    int i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx] / dn;}
    return res;
}

double nk_dsum(double* x, int n, int ldx){
    int i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

double nk_lsum2(double* x, double* p, int n, int ldx){
    int i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx]*log(p[i*ldx]);}
    return res;
}


char dim(gzFile f, int* pN, int* pP, int skip){// N x P matrix

    int N=0;
    int P=0;
    int i;
    char sep;
    int ret;
    char c;
    while((c=gzgetc(f)) != EOF){
        if(N==0+skip && (c=='\n'||c==' '||c=='\t'||c==',')){
            if(c!='\n'){sep=c;} 
            P++;
        }
        if(c=='\n'){N++;}
    }
    gzseek(f, 0L, SEEK_SET);
    (*pN)=N;
    (*pP)=P;
    return sep;
}

void qr(double* A, double* R, int M, int N, int lda){
    int lwork = N*N;
    int info;
    double* work; work = (double*)calloc(lwork, sizeof(double));
    double* tau;  tau  = (double*)calloc(N,     sizeof(double));
    dgeqrf_(&M, &N, A, &lda, tau, work, &lwork, &info);
    //fprintf(stderr, "dgeqrf=%d\n", info);
    
    int K = N;
    int i,j;
    for(i=0; i<N; i++){
        for(j=i; j<N; j++){
            R[i+j*N] = A[i+j*lda];
        }
    }
    //fprintf(stderr, "%d %d %d %d %d %d\n", M, N, K, lda, lwork, info);
    dorgqr_(&M, &N, &K, A, &lda, tau, work, &lwork, &info);
    //fprintf(stderr, "dorgqr=%d\n", info);
    free(work);
}

/*void clearAs0(double* x, int n){
    int i;
    for(i=0; i<n; i++){
        x[i] = 0.0;
    }
}*/

void BackSolve(double* R, int p, double* y, double* x){
    // To solve y = R x
    int i, j;
    clearAs0(x, p);
    
    x[p-1] = y[p-1]/R[p*p-1];
    
    double tm = 0.0;
    for(i=p-2; i>=0; i--){
        tm = 0.0;
        for(j=p-1; j>i; j--){
            tm += R[p*j+i]*x[j];
        }
        x[i] = (y[i]-tm)/R[p*i+i];
    }
}

int readTableInt(gzFile f, int* Y, int nrow, int ncol){
    int i, j, k=0;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t'){
                    cell[k] = '\0';
                    sscanf(cell, "%d", Y+i+j*nrow);
                    k = 0;
                    break;
                }else{
                    cell[k++] = c;
                }
            }
            if(c=='\0'){
                cell[k] = '\0';
                sscanf(cell, "%d", Y+i+j*nrow);
                k = 0;
                break;
            }
        }
    }
    return 0;
}

int readTable(gzFile f, int* Y, double* bf, int nrow, int ncol){
    int i, j, k=0;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol-1; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t'){
                    cell[k] = '\0';
                    sscanf(cell, "%d", Y+i+j*nrow);
                    k = 0;
                    break;
                }else{
                    cell[k++] = c;
                }
            }
        }
        while((c=gzgetc(f)) != EOF){
            if(c=='\n' || c=='\t'){
                cell[k] = '\0';
                sscanf(cell, "%lf", bf+i);
                bf[i] = exp(bf[i]);
                k = 0;
                break;
            }else{
                cell[k++] = c;
            }
        }
        if(c=='\0'){
            cell[k] = '\0';
            sscanf(cell, "%lf", bf+i);
            bf[i] = exp(bf[i]);
            k = 0;
            break;
        }
    }
    return 0;
}



void* Estep(void *args){
    HIERARCHICAL_MT* pmt = (HIERARCHICAL_MT *)args;
    
    int tid = pmt->tid;
    int npeaksperthread = pmt->npeaksperthread;
    
    int nvars; nvars = pmt->nvars;
    int* cumcats; cumcats = pmt->cumcats;
    int P = cumcats[nvars];
    
    int npeaks = pmt->npeaks;
    int* cumloci; cumloci = pmt->cumloci;
    int L = cumloci[npeaks];
    
    //fprintf(stderr, "L=%d P=%d nvars=%d tid=%d nppt=%d\n", L, P, nvars, tid, npeaksperthread);
    
    double* X;  X  = pmt->X;
    double* BF; BF = pmt->BF;
    
    double* beta; beta = pmt->beta;
    double* Pi;   Pi   = pmt->Pi;
    
    double* eta; eta = pmt->eta;
    double* pjk; pjk = pmt->pjk;
    double* z;   z   = pmt->z;
    double* Z1;  Z1  = pmt->Z1;
    double* w;   w   = pmt->w;
    double* y;   y   = pmt->y;
    double* Xt;  Xt  = pmt->Xt;
    
    double tot;
    
    double* xb; xb = (double*)calloc(P, sizeof(double));
    
    int j, k, l, pstart = tid*npeaksperthread, pend;
    pend = ((tid+1)*npeaksperthread < npeaks) ? (tid+1)*npeaksperthread : npeaks;
    for(j=pstart; j<pend; j++){
        //fprintf(stderr, "%d-%d\n", cumloci[j], cumloci[j+1]);
        //fprintf(stderr, "%dth peak %d %d\n", j, pstart, pend);
        // eta <- X, beta
        cblas_dgemv(CblasColMajor, CblasNoTrans, cumloci[j+1]-cumloci[j], P, 1.0, X+cumloci[j], L, beta, 1, 0.0, eta+cumloci[j], 1);
        // pjk <- softmax(eta)
        softmax(eta+cumloci[j], pjk+cumloci[j], cumloci[j+1]-cumloci[j]);
        //for(k=cumloci[j]; k<cumloci[j]+10; k++){fprintf(stderr, "%lf ", pjk[k]);}fprintf(stderr, "\n");
        // z   <- pjk, BF, Pi; &
        // Z1  <- z
        tot = Pi[0];
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            z[k] = Pi[1]*pjk[k]*BF[k];
            tot += Pi[1]*pjk[k]*BF[k];
        }
        Z1[j] = 0.0;
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            z[k] /= tot;
            Z1[j] += z[k];
        }
        // w   <- pjk, Z1
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            w[k] = Z1[j] * pjk[k];
        }
        //if(j==5){for(k=cumloci[j]; k<cumloci[j+1]; k++){fprintf(stderr, "%lf %lf ", w[k], Z1[j]);}fprintf(stderr, "\n");}
        // xb  <- X, pjk
        cblas_dgemv(CblasColMajor, CblasTrans, cumloci[j+1]-cumloci[j], P, 1.0, X+cumloci[j], L, pjk+cumloci[j], 1, 0.0, xb, 1);
        // Xt  <- X, xb, w; &
        // y   <- w, eta, z
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            y[k] = sqrt(w[k])*eta[k] + z[k]/sqrt(w[k]);
            for(l=0; l<P; l++){
                Xt[k+l*(L+P)] = (X[k+l*L]-xb[l]) * sqrt(w[k]);
            }
        }
    }
    
    free(xb);
    
    pthread_exit(NULL);
    return args;
}

void printV(double* x, int n){
    int i;
    for(i=0; i<n; i++){
        fprintf(stderr, "%lf ", x[i]);
    }
    fprintf(stderr, "\n");
}
void Mstep(double* Xt, double* R, double* y, double* beta, int L, int P, double sigma){
    int offs = 1;
    // prior on beta
    int i, j;
    for(i=0; i<P-offs; i++){
        y[L+i]=0.0;
        for(j=0; j<P-offs; j++){
            if(i==j){
                Xt[L+i+(L+P)*(offs+j)] = 1.0/sigma; 
            }else{
                Xt[L+i+(L+P)*(offs+j)] = 0.0;
            }
        }
    }
    // QR decomposition
    qr(Xt+(L+P)*offs, R, L+P-offs, P-offs, L+P);
    // beta1 <- t(Xt) %*% y
    cblas_dgemv(CblasColMajor, CblasTrans, L+P-offs, P-offs, 1.0, Xt+(L+P)*offs, L+P, y, 1, 0.0, beta+P+offs, 1);
    // beta  <- backsolve(R, beta1)
    BackSolve(R, P-offs, beta+P+offs, beta+offs);
    double th = 40.;
    for(i=1; i<P; i++){if(beta[i]>th){beta[i]=th;}else if(beta[i]<(-th)){beta[i]=-th;}}
}

double getQval(double* Z1, double* Pi, double* z, double* pjk, double* BF, int L, int npeaks){
    double Z1t = nk_dsum(Z1, npeaks, 1);
    double nd  = (double) npeaks;
    return Z1t*log(Pi[0]) + (nd-Z1t)*log(Pi[1]) + nk_lsum2(z, pjk, L, 1) + nk_lsum2(z, BF, L, 1);
}

// nvars = 1 (intercept) + # categocial variables
void em(double* X, double* BF, int L, int P, int npeaks, int nvars, int* cumloci, int* cumcats){
    int i, j, k;
    
    fprintf(stderr, "em: L=%d P=%d npeaks=%d nvars=%d\n", L, P, npeaks, nvars);
    
    fprintf(stderr, "Model fitting started...\n\n");
    
    // parameters
    double sigma = 100.;
    double Pi[2] = {0.9, 0.1};
    double* beta; beta = (double*)calloc(P*2, sizeof(double)); // coef length P + working space for another P
    // coproducts
    double* eta;  eta  = (double*)calloc(L, sizeof(double)); // X beta
    double* Z1;   Z1   = (double*)calloc(npeaks, sizeof(double)); // 1-Z0j; j=1,...,npeaks
    double* z;    z    = (double*)calloc(L, sizeof(double)); // porsterior probs
    double* y;    y    = (double*)calloc(L, sizeof(double)); // pseudo data : W^1/2 %*% eta + W^1/2 %*% z
    double* w;    w    = (double*)calloc(L, sizeof(double)); // weights for IRLS : expand(Z1) * pjk
    double* pjk;  pjk  = (double*)calloc(L, sizeof(double));  // softmax prior
    double* Xt;   Xt   = (double*)calloc((L+P)*P, sizeof(double)); // normalized X & replaced by Q
    double* R;    R    = (double*)calloc(P*P, sizeof(double)); // upper tri
    fprintf(stderr, "Memory allocated...\n\n");
    
    int nthreads = 5;
    int tid;
    HIERARCHICAL_MT* pmt; pmt=(HIERARCHICAL_MT*)calloc(nthreads, sizeof(HIERARCHICAL_MT));
    pthread_t* pid;       pid=(pthread_t*)calloc(nthreads, sizeof(pthread_t));
    int npeaksperthread  = npeaks / nthreads + 1;
    
    int itr;
    beta[0] = 0.0;
    beta[1] = -0.612239;//-0.639883;
    beta[2] = 5.617711;//5.319518; // 4.3
    beta[3] = 8.390924;//7.517685; // 3.0
    beta[4] = 1.145478;//0.461211;

    //for(i=1; i<P; i++){beta[i] = (double)i/100.0;}
    
    double qval, qval1=-1.0e10;
    for(itr=0; itr<100; itr++){
        // Esteps
        for(tid=0; tid<nthreads; tid++){
            pmt[tid].tid    = tid;
            pmt[tid].npeaksperthread = npeaksperthread;
            pmt[tid].npeaks = npeaks;
            pmt[tid].nvars  = nvars;
            pmt[tid].cumloci = cumloci;
            pmt[tid].cumcats = cumcats;
            pmt[tid].X = X;
            pmt[tid].BF = BF;
            pmt[tid].beta = beta;
            pmt[tid].eta  = eta;
            pmt[tid].Pi   = Pi;
            pmt[tid].pjk  = pjk;
            pmt[tid].z    = z;
            pmt[tid].Z1   = Z1;
            pmt[tid].w    = w;
            pmt[tid].y    = y;
            pmt[tid].Xt   = Xt;
            
            //Estep(pmt);
            
            int pthflag;
            if( (pthflag = pthread_create(pid+tid, NULL, (void*)Estep, (void*)(pmt+tid))) != 0){
                fprintf(stderr, "Thread not created...aborted.\n");
                return;
            }
        }
        for(tid=0; tid<nthreads; tid++){
            int pthflag;
            if( (pthflag = pthread_join(pid[tid], NULL)) !=0 ){fprintf(stderr, "Thread not joined...aborted.\n"); return ;};
        }
        
        // Msteps
        Mstep(Xt, R, y, beta, L, P, sigma);
        Pi[1] = nk_mean(Z1, npeaks, 1);
        Pi[0] = 1.0 - Pi[1];
        
        // Print res
        fprintf(stderr, "[%d] ", itr);
        for(i=0; i<P; i++){fprintf(stderr, "%lf ", beta[i]);}
        fprintf(stderr, "\n");
        fprintf(stderr, "Pi0=%lf %lf\n", Pi[0], Pi[1]);
        cblas_dgemv(CblasColMajor, CblasTrans, L, P, 1.0, X, L, z, 1, 0.0, beta+P, 1);
        for(i=0; i<P; i++){fprintf(stderr, "%lf ", beta[i+P]);}fprintf(stderr, "\n");
        
        qval = getQval(Z1, Pi, z, pjk, BF, L, npeaks);
        if(fabs(qval-qval1)<1.0e-7){break;}else{fprintf(stderr, "qval=%lf %lf\n", qval1, qval); qval1=qval;}
        
    }
    gzFile postf; postf = gzopen("posterior.gz", "wb6f");
    gzFile zf;    zf    = gzopen("Z1.gz", "wb6f");
    gzFile paraf; paraf = gzopen("param.gz", "wb6f");
    gzFile propf; propf = gzopen("proportion.gz", "wb6f");
    for(j=0; j<npeaks; j++){
        gzprintf(zf, "%lf\n", Z1[j]);
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            gzprintf(postf, "%d\t%lf\t%lf\t%lf\n", j+1, pjk[k], BF[k], z[k]);
        }
    }
    gzprintf(paraf, "%lf\t", Pi[1]);
    for(j=1; j<P-1; j++){ gzprintf(paraf, "%lf\t", beta[j]); }; gzprintf(paraf, "%lf\n", beta[P-1]);
    for(j=1; j<nvars; j++){
        char sep = '\t';
        if(j==nvars-1){sep = '\n';}
        double tot=0.0;
        for(k=cumcats[j]; k<cumcats[j+1]; k++){
            tot += beta[P+k];
            gzprintf(propf, "%lf\t", beta[P+k]);
        }
        gzprintf(propf, "%lf%c", beta[P]-tot, sep);
    }
    gzclose(paraf);
    gzclose(postf);
    gzclose(propf);
    gzclose(zf);
}

int main0(){
    double X[15] = {-0.09127699,-0.01028609,1.138874,0.4202284,-0.7107422,-0.5183098,-0.3752468,1.404647,0.306122,-0.5591557,0.2567613,-1.369183,0.6036138,1.117551,-0.4345506};
    double y[5] = {-0.3576127,-4.868328,5.75901,4.385124,-3.132706};
    double R[9] = {0,0,0,0,0,0,0,0,0};
    double beta[6] = {0,0,0,0,0,0};
    Mstep(X, R, y, beta, 5, 3, 100.0);
    int i,j;
    for(i=0; i<5; i++){
        for(j=0; j<3; j++){
            printf("%lf ", X[i+5*j]);
        }
        printf("\n");
    }
    printf("\n");
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            printf("%lf ", R[i+3*j]);
        }
        printf("\n");
    }
    printf("\n");
    printf("%lf %lf %lf\n", beta[0], beta[1], beta[2]);
    printf("%lf %lf %lf\n", beta[3], beta[4], beta[5]);
    return 0;
}

int main(int argc, char** argv){
    
    gzFile f = gzopen(argv[1], "rb6f");
    int nrow, ncol;
    dim(f, &nrow, &ncol, 0);
    fprintf(stderr, "\nFile name: %s\n", argv[1]);
    fprintf(stderr, "File size: %d x %d\n\n", nrow, ncol);
    
    
    int* X0=NULL; X0 = (int*)calloc(nrow*ncol,sizeof(int));
    double* bf=NULL; bf = (double*)calloc(nrow, sizeof(double));
    if(X0!=NULL && bf!=NULL){
        fprintf(stderr, "Memory allocated for data.\n\n");
    }else{fprintf(stderr, "Memory alloc error!\n"); return 1;}
    readTable(f, X0, bf, nrow, ncol);
    ncol--; // bayes factor column is diminished
    
    // counting levels for each var
    int i, j;
    fprintf(stderr, "Data file read.  First 30 rows:\n");
    int* ncats; ncats = (int*)calloc(ncol, sizeof(int));
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            if(i<30)fprintf(stderr, "%d\t", X0[i+j*nrow]);
            ncats[j] = ncats[j] < X0[i+j*nrow]+1 ? X0[i+j*nrow]+1 : ncats[j]; // number of levels
        }
        if(i<30)fprintf(stderr, "%lf\n", bf[i]);
    }
    ncats[0]--; // peak ID MUST start from 1
    
    for(i=0; i<nrow; i++){if(X0[i+2*nrow]>2){X0[i+2*nrow]=0;}}
    //ncats[2] = 3;
    
    // counting peak ids
    fprintf(stderr, "\n");
    fprintf(stderr, "Data: %d peaks\n", ncats[0]);
    int* nloci;   nloci   = (int*)calloc(ncats[0],   sizeof(int));
    int* cumloci; cumloci = (int*)calloc(ncats[0]+1, sizeof(int));
    for(i=0; i<nrow; i++){
        nloci[X0[i]-1]++;
    }
    for(i=0; i<ncats[0]; i++){
        cumloci[i+1] = cumloci[i]+nloci[i];
    }
    //for(i=0; i<ncats[0]+1; i++){ fprintf(stderr, "%d\t", cumloci[i]); }
    
    
    
    // prep for expanding model matrix
    int* cumcats; cumcats = (int*)calloc(ncol+1,     sizeof(int));
    cumcats[1]=1;
    for(j=1; j<ncol; j++){
        fprintf(stderr, "%dth variable: %d levels\n", j, ncats[j]);
        cumcats[j+1] = cumcats[j] + (ncats[j]-1);   // number of cumlative levels
    }
    fprintf(stderr, "\n");
    
    // model matrix encoding    
    double* X; X=(double*)calloc(cumcats[ncol] * nrow, sizeof(double));
    fprintf(stderr, "%d x %d matrix memory allocated.\n\n", nrow, cumcats[ncol]);
    for(i=0; i<nrow; i++){
        X[i] = 1.0;
        for(j=1; j<ncol; j++){   
            if(X0[i+j*nrow]>0){ X[i+(X0[i+j*nrow]-1+cumcats[j])*nrow] = 1.0; } // encoding
        }
    }
    //printing model matrix
    /*for(i=0; i<nrow; i++){// print X
        fprintf(stderr, "%d %d %d:\t", X0[i], nloci[X0[i]-1], cumloci[X0[i]-1]);
        for(j=0; j<cumcats[ncol]; j++){
            fprintf(stderr, "%.0lf\t", X[i+j*nrow]);
        }
        fprintf(stderr, "\n");
    }*/
    
    //fprintf(stderr, "main: L=%d P=%d npeaks=%d nvars=%d\n", nrow, cumcats[ncol], ncats[0], ncol);
    em(X, bf, nrow, cumcats[ncol], ncats[0], ncol, cumloci, cumcats);
    
}









