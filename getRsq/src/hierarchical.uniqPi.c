#include "hierarchical.h"
#include "util.h"

void printV(double* x, long n);

double rk1(double x, double z){
    return ((z-0.5)*(z-0.5)-1./12.)*((x-0.5)*(x-0.5)-1./12.)/4. - (pow(fabs(x-z)-0.5, 4.)-(fabs(x-z)-0.5)*(fabs(x-z)-0.5)/2. + 7./240.)/24.;
}
void rk(double* x0, double* x, double* xk, long n, long nk){
    long i, j;
    for(i=0; i<n; i++){
        if(x0[i]>=0.0){
            x[i] = x0[i];
            for(j=0; j<nk; j++){
                x[(j+1)*n+i] = rk1(x0[i], xk[j]);
            }
        }
    }
    
}

double nk_mean(double* x, long n, long ldx){
    long i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx] / dn;}
    return res;
}

double nk_dsum(double* x, long n, long ldx){
    long i;
    double res=0.0;
    double dn = (double)n;
    for(i=0; i<n; i++){res += x[i*ldx];}
    return res;
}

double nk_lsum2(double* x, double* p, long n, long ldx){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i*ldx]*log(p[i*ldx]);}
    return res;
}


char dim(gzFile f, long* pN, long* pP, long skip){// N x P matrix

    long N=0;
    long P=0;
    long i;
    char sep;
    long ret;
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

void qr(double* A, double* R, long M, long N, long lda){
    long lwork = N*N;
    long info;
    double* work; work = (double*)calloc(lwork, sizeof(double));
    double* tau;  tau  = (double*)calloc(N,     sizeof(double));
    dgeqrf_(&M, &N, A, &lda, tau, work, &lwork, &info);
    //fprintf(stderr, "dgeqrf=%d\n", info);
    
    long K = N;
    long i,j;
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

/*void clearAs0(double* x, long n){
    long i;
    for(i=0; i<n; i++){
        x[i] = 0.0;
    }
}*/

void BackSolve(double* R, long p, double* y, double* x){
    // To solve y = R x
    long i, j;
    clearAs0(x, p);
   
    /*fprintf(stderr, "\n");
    printV(y, p);
    fprintf(stderr, "\n");
    for(i=0; i<p; i++){printV(R+p*i, p);}
    fprintf(stderr, "\n");*/

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

long readTableInt(gzFile f, long* Y, long nrow, long ncol){
    long i, j, k=0;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t'){
                    cell[k] = '\0';
                    sscanf(cell, "%ld", Y+i+j*nrow);
                    k = 0;
                    break;
                }else{
                    cell[k++] = c;
                }
            }
            if(c=='\0'){
                cell[k] = '\0';
                sscanf(cell, "%ld", Y+i+j*nrow);
                k = 0;
                break;
            }
        }
    }
    return 0;
}

long readTable(gzFile f, long* Y, double* bf, long nrow, long ncol, long LOG){
    long i, j, k=0;
    char c;
    char* cell; cell = (char*)calloc(1000, sizeof(char));
    gzseek(f, 0L, SEEK_SET);
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol-1; j++){
            while((c=gzgetc(f)) != EOF){
                if(c=='\n' || c=='\t'){
                    cell[k] = '\0';
                    sscanf(cell, "%ld", Y+i+j*nrow);
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
//if(i<10)fprintf(stderr, "%s %lf\n", cell, bf[i]);
                if(LOG>0) bf[i] = exp(bf[i]);
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

void dgemv(double* A, long N, long M, long lda, double* x, double* y){
    long i, j;
    for(i=0; i<N; i++){
        y[i] = 0.0;
        for(j=0; j<M; j++){
            y[i] += A[i+j*lda]*x[j];
        }
    }
}


void dgemvT(double* A, long N, long M, long lda, double* x, double* y){
    long i, j;
    for(j=0; j<M; j++){
        y[j] = 0.0;
        for(i=0; i<N; i++){
            y[j] += A[i+j*lda]*x[i];
        }
    }
}


void* Estep(void *args){
    HIERARCHICAL_MT* pmt = (HIERARCHICAL_MT *)args;
    
    long tid = pmt->tid;
    long npeaksperthread = pmt->npeaksperthread;
    
    //long nvars; nvars = pmt->nvars;
    //long* cumcats; cumcats = pmt->cumcats;
    long P = pmt->P;//cumcats[nvars];
    
    long npeaks = pmt->npeaks;
    long* cumloci; cumloci = pmt->cumloci;
    long L = cumloci[npeaks];
    
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
    double* Ie;  Ie  = pmt->Ie;  clearAs0(Ie, P*P);
    
    double tot;
    
    double* xb; xb = (double*)calloc(P, sizeof(double));
    double* XjTZj; XjTZj = (double*)calloc(P, sizeof(double));
    
    long j, k, l, l2, pstart = tid*npeaksperthread, pend;
    pend = ((tid+1)*npeaksperthread < npeaks) ? (tid+1)*npeaksperthread : npeaks;
    for(j=pstart; j<pend; j++){
        //fprintf(stderr, "%d-%d\n", cumloci[j], cumloci[j+1]);
        //fprintf(stderr, "%dth peak %d %d\n", j, pstart, pend);
        // eta <- X, beta
#ifdef CBLAS
        fprintf(stderr, "blas\n");
        cblas_dgemv(CblasColMajor, CblasNoTrans, cumloci[j+1]-cumloci[j], P, 1.0, X+cumloci[j], L, beta, 1, 0.0, eta+cumloci[j], 1);
#else
        dgemv(X+cumloci[j], cumloci[j+1]-cumloci[j], P, L, beta, eta+cumloci[j]);
#endif
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
#ifdef CBLAS
        cblas_dgemv(CblasColMajor, CblasTrans, cumloci[j+1]-cumloci[j], P, 1.0, X+cumloci[j], L, pjk+cumloci[j], 1, 0.0, xb, 1);
#else
        dgemvT(X+cumloci[j], cumloci[j+1]-cumloci[j], P, L, pjk+cumloci[j], xb);
#endif
        // Xt  <- X, xb, w; &
        // y   <- w, eta, z
        // XjTZj <- X, z
        clearAs0(XjTZj, P);
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            // pseudo data
            y[k] = sqrt(w[k])*eta[k] + z[k]/sqrt(w[k]);
            for(l=0; l<P; l++){
                Xt[k+l*(L+P)] = (X[k+l*L]-xb[l]) * sqrt(w[k]);
                XjTZj[l]     += (X[k+l*L]-xb[l]) * z[k];
            }
        }
        // empirical hessian ~= Ic - Io
        for(l=0; l<P; l++){
            for(l2=0; l2<P; l2++){
                Ie[l+l2*P] += XjTZj[l]*XjTZj[l2];
            }
        }
    }
    
    free(XjTZj);
    free(xb);
    
    pthread_exit(NULL);
    return args;
}

void printV(double* x, long n){
    long i;
    for(i=0; i<n; i++){
        fprintf(stderr, "%lf ", x[i]);
    }
    fprintf(stderr, "\n");
}
int Itr;
void Mstep(double* Xt, double* R, double* y, double* beta, long L, long P, double sigma){
    double xk[4] = {0.1470588, 0.3823529, 0.6176471, 0.8529412};
    double lambda = 5.;
    long offs = 1;
    // prior on beta
    fprintf(stderr, "penalty ");
    long i, j;
    for(i=0; i<P-offs; i++){
        y[L+i]=0.0;
        for(j=0; j<P-offs; j++){
            if(i==j){
                Xt[L+i+(L+P)*(offs+j)] = 1.0/sigma; 
            }else{
                Xt[L+i+(L+P)*(offs+j)] = 0.0;
            }
            if(P>5){if(i>4 && j>4 && j>=i){
                Xt[L+i+(L+P)*(offs+j)] = lambda*rk1(xk[i-5],xk[j-5]);
            }}
        }
    }
    if(P>5){
        fprintf(stderr, "cholesky ");
        // Cholesky decomp
        char uplo = 'U';
        int n=4, lda=L+P, info;
        dpotrf_(&uplo, &n, Xt+(L+P)*6+L+5, &lda, &info);
    }
    // QR decomposition
    qr(Xt+(L+P)*offs, R, L+P-offs, P-offs, L+P);
    // beta1 <- t(Xt) %*% y
#ifdef CBLAS
    cblas_dgemv(CblasColMajor, CblasTrans, L+P-offs, P-offs, 1.0, Xt+(L+P)*offs, L+P, y, 1, 0.0, beta+P+offs, 1);
#else
    dgemvT(Xt+(L+P)*offs, L+P-offs, P-offs, L+P, y, beta+P+offs);
#endif
    // beta  <- backsolve(R, beta1)
    BackSolve(R, P-offs, beta+P+offs, beta+offs);
    double th = 15.;
    //for(i=1; i<P; i++){if(beta[i]>th){beta[i]=beta0[i];}else if(beta[i]<(-th)){beta[i]=beta0[i];}}
    if(Itr<2){for(i=1; i<P; i++){beta[i]=beta0[i];}}
}

double getQval(double* Z1, double* Pi, double* z, double* pjk, double* BF, long L, long npeaks){
    double Z1t = nk_dsum(Z1, npeaks, 1);
    double nd  = (double) npeaks;
    return Z1t*log(Pi[0]) + (nd-Z1t)*log(Pi[1]) + nk_lsum2(z, pjk, L, 1) + nk_lsum2(z, BF, L, 1);
}

// nvars = 1 (longercept) + # categocial variables
void em(double* X, double* BF, long L, long P, long npeaks, long nvars, long* cumloci, long* cumcats){
    long i, j, k, l, l2;
    
    fprintf(stderr, "em: L=%ld P=%ld npeaks=%ld nvars=%ld\n", L, P, npeaks, nvars);
    
    fprintf(stderr, "Model fitting started...\n\n");
    
    
    long nthreads = 5;
    long tid;
    HIERARCHICAL_MT* pmt; pmt=(HIERARCHICAL_MT*)calloc(nthreads, sizeof(HIERARCHICAL_MT));
    pthread_t* pid;       pid=(pthread_t*)calloc(nthreads, sizeof(pthread_t));
    long npeaksperthread  = npeaks / nthreads + 1;
    
    // parameters
    double sigma = 100.;
    double Pi[2] = {0.919806, 0.080194};
    double* beta; beta = (double*)calloc(P*2, sizeof(double)); // coef length P + working space for another P
    double* beta_old; beta_old = (double*)calloc(P, sizeof(double)); // coef length P + working space for another P
    beta0 = (double*)calloc(P*2, sizeof(double)); // coef length P + working space for another P
    // coproducts
    double* eta;  eta  = (double*)calloc(L, sizeof(double)); // X beta
    double* Z1;   Z1   = (double*)calloc(npeaks, sizeof(double)); // 1-Z0j; j=1,...,npeaks
    double* z;    z    = (double*)calloc(L, sizeof(double)); // porsterior probs
    double* y;    y    = (double*)calloc(L, sizeof(double)); // pseudo data : W^1/2 %*% eta + W^1/2 %*% z
    double* w;    w    = (double*)calloc(L, sizeof(double)); // weights for IRLS : expand(Z1) * pjk
    double* pjk;  pjk  = (double*)calloc(L, sizeof(double));  // softmax prior
    double* Xt;   Xt   = (double*)calloc((L+P)*P, sizeof(double)); // normalized X & replaced by Q
    double* R;    R    = (double*)calloc(P*P, sizeof(double)); // upper tri
    double* Ie;   Ie   = (double*)calloc(P*P*nthreads, sizeof(double));
    fprintf(stderr, "Memory allocated...\n\n");
    
    long itr;
    beta[0] = 0.0;//  0.000000 -0.685714 5.903557 9.346523 2.244446
    beta[1] = beta0[1] =-0.91;     //-0.639883;
    beta[2] = beta0[2] = 5.133316; //5.319518; // 4.3
    //beta[3] = beta0[3] = 1.828196; //7.517685; // 3.0
    beta[3] = beta0[3] = 9.828196; //7.517685; // 3.0
    if(P>4)beta[4] = beta0[4] = 2.1;      //0.461211;
    
    /*beta[1] = -0.827825;
    beta[2] = 5.177595;
    beta[3] = 9.308916;
    beta[3] = 2.0;
    if(P>4)beta[4] = 1.946503;*/
    
    
    if(P>5){
        beta0[1] = beta[1] = -0.482542;
        beta0[2] = beta[2] =  5.463498;
        beta0[3] = beta[3] =  7.017382;
        beta0[4] = beta[4] = -1.881284;
        beta0[5] = beta[5] =  4.303588 ;
        beta0[6] = beta[6] = -1107.495457;
        beta0[7] = beta[7] = -1400.326647;
        beta0[8] = beta[8] = -560.874219;
        beta0[9] = beta[9] = -1895.056297;
    }

    //for(i=1; i<P; i++){beta[i] = (double)i/100.0;}
    int again=0;
    double qval, qval1=-1.0e10;
    for(itr=0; itr<1000; itr++){
        Itr=itr;
        // Esteps
        for(tid=0; tid<nthreads; tid++){
            pmt[tid].tid    = tid;
            pmt[tid].npeaksperthread = npeaksperthread;
            pmt[tid].npeaks = npeaks;
            pmt[tid].nvars  = nvars;
            pmt[tid].cumloci = cumloci;
            pmt[tid].cumcats = cumcats;
            pmt[tid].P = P;
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
            pmt[tid].Ie   = Ie+P*P*tid;
            
            //Estep(pmt);
            
            long pthflag;
            if( (pthflag = pthread_create(pid+tid, NULL, (void*)Estep, (void*)(pmt+tid))) != 0){
                fprintf(stderr, "Thread not created...aborted.\n");
                return;
            }
        }
        for(tid=0; tid<nthreads; tid++){
            long pthflag;
            if( (pthflag = pthread_join(pid[tid], NULL)) !=0 ){fprintf(stderr, "Thread not joined...aborted.\n"); return ;};
        }
        
        qval = getQval(Z1, Pi, z, pjk, BF, L, npeaks);
        fprintf(stderr, "qval=%lf -> %lf\n", qval1, qval); 
        if(qval1-qval>1000.0 && again==0){
            fprintf(stderr, "devol!\n");
            fprintf(stderr, "Before: "); printV(beta, P);
            for(i=0; i<P; i++){ beta[i] = beta_old[i] + (beta[i]-beta_old[i])/100.0; }
            fprintf(stderr, "After : "); printV(beta, P);
            again=1;
        }else{
            
            // Msteps
            cblas_dcopy(P, beta, 1, beta_old, 1);
            fprintf(stderr, "Mstep\n");
            Mstep(Xt, R, y, beta, L, P, sigma);
            for(tid=1; tid<nthreads; tid++){for(l=0; l<P; l++){for(l2=0; l2<P; l2++){Ie[l+l2*P]+=Ie[l+l2*P+tid*P*P];}}}
            Pi[1] = nk_mean(Z1, npeaks, 1);
            Pi[0] = 1.0 - Pi[1];
            
            
            
            
            // Print res
            fprintf(stderr, "[%ld] ", itr);
            for(i=0; i<P; i++){fprintf(stderr, "%lf ", beta[i]);}
            fprintf(stderr, "\n");
            fprintf(stderr, "Pi0=%lf %lf\n", Pi[0], Pi[1]);
#ifdef CBLAS
            cblas_dgemv(CblasColMajor, CblasTrans, L, P, 1.0, X, L, z, 1, 0.0, beta+P, 1);
#else
            dgemvT(X, L, P, L, z, beta+P);
#endif
            double hoge=0.0;
            for(i=0; i<L; i++){hoge += X[i+(P-1)*L]*z[i]; if(hoge>(double)npeaks){fprintf(stderr, "Mul: %ld %lf %lf %lf\n", i, hoge, X[i+(P-1)*L], z[i]); break;}}
            for(i=0; i<P; i++){fprintf(stderr, "%lf ", beta[i+P]);}fprintf(stderr, "\n");
            
            // terminate
            if(fabs(qval-qval1)<1.0e-6 && again==0 && itr>10){break;}
            
            qval1=qval;
            again = 0;
            
        }
    }
    gzFile postf; postf = gzopen("posterior.gz", "wb6f");
    gzFile zf;    zf    = gzopen("Z1.gz", "wb6f");
    gzFile paraf; paraf = gzopen("param.gz", "wb6f");
    gzFile propf; propf = gzopen("proportion.gz", "wb6f");
    gzFile hessf; hessf = gzopen("ehess.gz", "wb6f");
    // posterior probs & Z1
    for(j=0; j<npeaks; j++){
        gzprintf(zf, "%lf\n", Z1[j]);
        for(k=cumloci[j]; k<cumloci[j+1]; k++){
            gzprintf(postf, "%ld\t%lf\t%lf\t%lf\n", j+1, pjk[k], BF[k], z[k]);
        }
    }
    // parameters
    gzprintf(paraf, "%lf\t", Pi[1]);
    for(j=1; j<P-1; j++){ gzprintf(paraf, "%lf\t", beta[j]); }; gzprintf(paraf, "%lf\n", beta[P-1]);
    // proportions
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
    // hessian
    for(l=1; l<P; l++){
        for(l2=1; l2<P; l2++){
            gzprintf(hessf, "%lf\t", Ie[l+l2*P]);
        }
        gzprintf(hessf, "\n");
    }
    gzclose(paraf);
    gzclose(postf);
    gzclose(propf);
    gzclose(zf);
    gzclose(hessf);
}

long main0(){
    double X[15] = {-0.09127699,-0.01028609,1.138874,0.4202284,-0.7107422,-0.5183098,-0.3752468,1.404647,0.306122,-0.5591557,0.2567613,-1.369183,0.6036138,1.117551,-0.4345506};
    double y[5] = {-0.3576127,-4.868328,5.75901,4.385124,-3.132706};
    double R[9] = {0,0,0,0,0,0,0,0,0};
    double beta[6] = {0,0,0,0,0,0};
    Mstep(X, R, y, beta, 5, 3, 100.0);
    long i,j;
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

long main(long argc, char** argv){
    long i, j;
    
    gzFile f = gzopen(argv[1], "rb6f");
    long nrow, ncol;
    dim(f, &nrow, &ncol, 0);
    fprintf(stderr, "\nFile name: %s\n", argv[1]);
    fprintf(stderr, "File size: %ld x %ld\n\n", nrow, ncol);
    
    
    long* X0=NULL; X0 = (long*)calloc(nrow*ncol,sizeof(long));
    double* bf=NULL; bf = (double*)calloc(nrow, sizeof(double));
    if(X0!=NULL && bf!=NULL){
        fprintf(stderr, "Memory allocated for data.\n\n");
    }else{fprintf(stderr, "Memory alloc error!\n"); return 1;}
    readTable(f, X0, bf, nrow, ncol, 1);
    ncol--; // bayes factor column is diminished
    
    // counting levels for each var
    fprintf(stderr, "Data file read.  First 30 rows:\n");
    long* ncats; ncats = (long*)calloc(ncol, sizeof(long));
    for(i=0; i<nrow; i++){
        for(j=0; j<ncol; j++){
            if(i<30)fprintf(stderr, "%ld\t", X0[i+j*nrow]);
            ncats[j] = ncats[j] < X0[i+j*nrow]+1 ? X0[i+j*nrow]+1 : ncats[j]; // number of levels
        }
        if(i<30)fprintf(stderr, "%lf\n", bf[i]);
    }
    ncats[0]--; // peak ID MUST start from 1
    
    //for(i=0; i<nrow; i++){if(X0[i+2*nrow]>2){X0[i+2*nrow]=0;}}
    //ncats[2] = 3;
    
    // counting peak ids
    fprintf(stderr, "\n");
    fprintf(stderr, "Data: %ld peaks\n", ncats[0]);
    long* nloci;   nloci   = (long*)calloc(ncats[0],   sizeof(long));
    long* cumloci; cumloci = (long*)calloc(ncats[0]+1, sizeof(long));
    for(i=0; i<nrow; i++){
        nloci[X0[i]-1]++;
    }
    for(i=0; i<ncats[0]; i++){
        cumloci[i+1] = cumloci[i]+nloci[i];
    }
    //for(i=0; i<ncats[0]+1; i++){ fprintf(stderr, "%d\t", cumloci[i]); }
    
    
    
    // prep for expanding model matrix
    long* cumcats; cumcats = (long*)calloc(ncol+1,     sizeof(long));
    cumcats[1]=1;
    for(j=1; j<ncol; j++){
        fprintf(stderr, "%ldth variable: %ld levels\n", j, ncats[j]);
        cumcats[j+1] = cumcats[j] + (ncats[j]-1);   // number of cumlative levels
    }
    fprintf(stderr, "\n");
    
    
    // spline smoothing vector
    double* xs;
    long ncolspline = 0;
    double xk[4] = {0.1470588, 0.3823529, 0.6176471, 0.8529412};
    for(i=0; i<argc-1; i++){
        if(strcmp(argv[i], "--spline")==0){
            ncolspline = 5;
            xs = (double*)calloc(nrow, sizeof(double));
            gzFile fs = gzopen(argv[i+1], "rb6f");
            readTable(fs, NULL, xs, nrow, 1, 0);
            printV(xs, 10);
        }
    }
    
    // model matrix encoding    
    double* X=NULL; X=(double*)calloc((cumcats[ncol]+ncolspline) * nrow, sizeof(double));
    if(X!=NULL) fprintf(stderr, "%ld x %ld matrix memory allocated.\n\n", nrow, cumcats[ncol]);
    for(i=0; i<cumcats[ncol]*nrow; i++){X[i] = 0.0;}
    for(i=0; i<nrow; i++){
        X[i] = 1.0;
        for(j=1; j<ncol; j++){   
            if(X0[i+j*nrow]>0){ X[i+(X0[i+j*nrow]-1+cumcats[j])*nrow] = 1.0; } // encoding
        }
    }
    if(ncolspline>0){
        rk(xs, X+cumcats[ncol]*nrow, xk, nrow, (long)(ncolspline-1));
        for(i=0; i<100; i++){for(j=0; j<cumcats[ncol]+ncolspline; j++){fprintf(stderr, "%lf ", X[i+j*nrow]);} fprintf(stderr, "\n");}
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
    //em(X, BF,  L, P,                        npeaks,   nvars,cumloci, cumcats)
    em(X, bf, nrow, cumcats[ncol]+ncolspline, ncats[0], ncol, cumloci, cumcats);
    
}









