#include "util.h"

const gsl_rng_type * rngT;
gsl_rng * rng;

int mini(int a, int b){return a<b ? a : b;}
int maxi(int a, int b){return a>b ? a : b;}

double logit(double x){
    return log(x/(1.-x));
}
double expit(double x){
    return 1./(1.+exp(-x));
}
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



void expand(int* pos, double* val, int n, int* pos2, double* val2, int n2, double* w){// n < n2
    int i=0, j=0;
    while(i<n){
        if(pos2[j]<pos[i]){
            j++;
        }else if(pos[i]<pos2[j]){
            i++;
        }else if(pos[i]==pos2[j]){
            val2[j] = val[i];
            w[j] = 1.0;
            i++;
            j++;
        }
    }
}


void softmax(double* eta, double* y, int n){
    int i;
    double offs=0.0, tot=0.0;
    for(i=0; i<n; i++){
        if(offs < eta[i]-20.0){offs = eta[i]-20.0;}
    }
    for(i=0; i<n; i++){
        y[i] = exp(eta[i] - offs);
        tot += y[i];
    }
    for(i=0; i<n; i++){
        y[i] /= tot;
    }
}

void wsoftmax(double* eta, double* y, double* w, int n){
    int i;
    double offs=0.0, tot=0.0;
    for(i=0; i<n; i++){if(w[i]>0.5){
        if(offs < eta[i]-20.0){offs = eta[i]-20.0;}
    }}
    for(i=0; i<n; i++){if(w[i]>0.5){
        y[i] = exp(eta[i] - offs);
        tot += y[i];
    }}
    for(i=0; i<n; i++){if(w[i]>0.5){
        y[i] /= tot;
    }else{y[i] = 0.0;}}
}

void clearAs0(double* x, int n){
    int i;
    for(i=0; i<n; i++){x[i]=0.0;}
}

void sumTo1(double* x, int n){
    double tot=0.0;
    int i;
    for(i=0; i<n; i++){tot += x[i];}
    for(i=0; i<n; i++){x[i] /= tot;}
}

void getIntrPost(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* Psi, double* w, int n, double* pp2, double* pp13, int* lci){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double Pi1_a = 0.07862747919984960920381;
    double Pi0_a = 1.0-Pi1_a;
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta0, p12,     w, n);
    wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Psi[0] * Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Psi[1] * Pi0_a       * exp(-2.*geta); 
    pp13[6] = Psi[2] * Pi0_1       * exp(-2.*geta); 
    pp13[8] = Psi[3] * Pi0_2       * exp(-2.*geta);
    for(j=0; j<n; j++){
        //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Psi[0] * Pi1_1 * Pi0_2 * p12[j+n]   * exp(lbf1[j]-geta); //
        pp13[2] += Psi[0] * Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); //
        //pp13[3] // linkage
        
        //pp13[4]  = Pi0_a;
        pp13[5] += Psi[1] * Pi1_a * p12[j]     * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
        //pp13[6]  = Pi0_1;
        pp13[7] += Psi[2] * Pi1_1 * p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta*2.); // 1 -> 2
        //pp13[8]  = Pi0_2;
        pp13[9] += Psi[3] * Pi1_2 * p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta*2.); // 2 -> 1
        
        // same var
        diaglikel += Psi[0] * Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j+n]*p12[j+n*2];
        /*if(lci[j]>=0){for(k=1; k<j+1; k++){// vars in the same peak
         if(lci[j]==lci[j-k]){
         diaglikel += Pi1_1*Pi1_2 *p12[j+n-k]*p12[j+n*2] * exp(lbf1[j-k]+lbf2[j] - geta*2.);
         diagprior += p12[j+n-k]*p12[j+n*2];
         diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2-k] * exp(lbf1[j]+lbf2[j-k] - geta*2.);
         diagprior += p12[j+n]*p12[j+n*2-k];
         }else{
         break;
         }
         }}*/
    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2/Psi[0] - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 10);
    
    free(p12);
}

void getLeadVar(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci, double* maxp, int* maxj){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;

    wsoftmax(eta0, p12,     w, n);
    wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    clearAs0(maxp, 7);
    for(j=0; j<n; j++){
        if(maxp[1] < p12[j+n]   * exp(lbf1[j]-geta)){maxp[1]=p12[j+n]   * exp(lbf1[j]-geta); maxj[0]=j;}
        if(maxp[2] < p12[j+n*2] * exp(lbf2[j]-geta)){maxp[2]=p12[j+n*2] * exp(lbf2[j]-geta); maxj[1]=j;}
        if(maxp[4] < p12[j]     * exp(lbf1[j]+lbf2[j]  -geta*2.)){maxp[4]=p12[j]     * exp(lbf1[j]+lbf2[j]  -geta*2.); maxj[4]=j;}; // Pleio
        if(maxp[5] < p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta*2.)){maxp[5]=p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta*2.); maxj[5]=j;}; // j -> k
        if(maxp[6] < p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta*2.)){maxp[6]=p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta*2.); maxj[6]=j;}; // k -> j
        for(k=0; k<n; k++){
            if(k!=j){if(maxp[3] < p12[j+n]*p12[k+n*2]*exp(lbf1[j]-geta+lbf2[k]-geta)){
                maxp[3] = p12[j+n]*p12[k+n*2]*exp(lbf1[j]-geta+lbf2[k]-geta); 
                maxj[2] = j;
                maxj[3] = k;
            }}
        }
    }
    
    free(p12);
}


void pwhmnewWithPia(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double Pi1_0, double* w, int n, double* pp2, double* pp13, int* lci){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double Pi0_0 = 1.0-Pi1_0;
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    
    double geta1 = 0;
    double geta2 = 0;
    for(j=0; j<n; j++){
        if(geta1<lbf1[j]){geta1=lbf1[j];}
        if(geta1<lbf2[j]){geta1=lbf2[j];}
        if(geta2<lbfmr1[j]){geta2=lbfmr1[j];}
        if(geta2<lbfmr2[j]){geta2=lbfmr2[j];}
    }
    //geta /= 2.0;
    geta1 -= 10.0;
    geta2 -= 10.0;
    wsoftmax(eta0, p12,     w, n);
    wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-geta1-geta2);
    pp13[4] = Pi0_0       * exp(-geta1-geta2); 
    pp13[6] = Pi0_1       * exp(-geta1-geta2); 
    pp13[8] = Pi0_2       * exp(-geta1-geta2);
    for(j=0; j<n; j++){
      //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1_1 * Pi0_2 * p12[j+n]   * exp(lbf1[j]-(geta1+geta2)/2.); //
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-(geta1+geta2)/2.); //
      //pp13[3] // linkage
        
      //pp13[4]  = Pi0_0;
        pp13[5] += Pi1_0 * p12[j]     * exp(lbf1[j]+lbf2[j]  -geta1-geta2); // Pleio
      //pp13[6]  = Pi0_1;
        pp13[7] += Pi1_1 * p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta1-geta2); // 1 -> 2
      //pp13[8]  = Pi0_2;
        pp13[9] += Pi1_2 * p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta1-geta2); // 2 -> 1
        
        // same var
        diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta1-geta2);
        diagprior += p12[j+n]*p12[j+n*2];
if(isnan(pp13[7])>0 ){fprintf(stderr,"%d nan %lf %lf %lf %lf\n", j, lbf1[j], lbf2[j], lbfmr1[j], lbfmr2[j]); break;}

    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-(geta1+geta2)/2.);
    pp13[2] *= exp(-(geta1+geta2)/2.);
   
    sumTo1(pp13, 10);
    
    free(p12);
}


void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double Pi0_0 = 1.0;//-0.076026;
    double Pi1_0 = 1.0;//0.076026;
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;

    wsoftmax(eta0, p12,     w, n);
    wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_0       * exp(-2.*geta); 
    pp13[6] = Pi0_1       * exp(-2.*geta); 
    pp13[8] = Pi0_2       * exp(-2.*geta);
    for(j=0; j<n; j++){
      //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1_1 * Pi0_2 * p12[j+n]   * exp(lbf1[j]-geta); //
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); //
      //pp13[3] // linkage
        
      //pp13[4]  = Pi0_0;
        pp13[5] += Pi1_0 * p12[j]     * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[6]  = Pi0_1;
        pp13[7] += Pi1_1 * p12[j+n]   * exp(lbf1[j]+lbfmr1[j]-geta*2.); // 1 -> 2
      //pp13[8]  = Pi0_2;
        pp13[9] += Pi1_2 * p12[j+n*2] * exp(lbf2[j]+lbfmr2[j]-geta*2.); // 2 -> 1
        
        // same var
        diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j+n]*p12[j+n*2];
        /*if(lci[j]>=0){for(k=1; k<j+1; k++){// vars in the same peak
            if(lci[j]==lci[j-k]){
                diaglikel += Pi1_1*Pi1_2 *p12[j+n-k]*p12[j+n*2] * exp(lbf1[j-k]+lbf2[j] - geta*2.);
                diagprior += p12[j+n-k]*p12[j+n*2];
                diaglikel += Pi1_1*Pi1_2 *p12[j+n]*p12[j+n*2-k] * exp(lbf1[j]+lbf2[j-k] - geta*2.);
                diagprior += p12[j+n]*p12[j+n*2-k];
            }else{
                break;
            }
        }}*/
    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 10);
    
    free(p12);
}

void colocWithGWAS(double* lbfj, double* lbfk, double* eta0, double Pi1_j, double Pi1_k, double* w, int n, double* pp13){
    int j, k;
    double* p12; p12 = (double*)calloc(n, sizeof(double));

    double Pi0_j = 1.0 - Pi1_j;
    double Pi0_k = 1.0 - Pi1_k;
    double Pi0_a=1.0, Pi1_a=1.0;
    Pi0_a = Pi0_j;
    Pi1_a = Pi1_j;

    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbfj[j]){geta=lbfj[j];}
        if(geta<lbfk[j]){geta=lbfk[j];}
    }
    geta /= 2.0;
    geta -= 10.0;

    wsoftmax(eta0, p12, w, n);

    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_j*Pi0_k * exp(-2.*geta);
    pp13[4] = Pi0_a       * exp(-2.*geta);
    for(j=0; j<n; j++){
        pp13[1] += Pi1_j * Pi0_k * p12[j]          * exp(lbfj[j]-geta); // eqtl
        pp13[2] += Pi0_j * Pi1_k * p12[j]          * exp(lbfk[j]-geta); // atac
        pp13[5] += Pi1_a         * p12[j]          * exp(lbfj[j]+lbfk[j]  -geta*2.); // Pleio
        diaglikel += Pi1_j*Pi1_k * p12[j] * p12[j] * exp(lbfj[j]+lbfk[j] - geta*2.);
        diagprior += p12[j]*p12[j];

    }

    pp13[3] = (pp13[1]*pp13[2]/Pi0_j/Pi0_k - diaglikel)/(1.-diagprior); // linkage

    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);

    sumTo1(pp13, 6);

    free(p12);
}



// Param 5
void pwhmnewataceqtlAllParam(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    Pi1_1 = Pi1_2 = 1.0;
    double Pi0_1 = Pi1_1;
    double Pi0_2 = Pi1_2;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta0, p12,     w, n);
    //wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_2       * exp(-2.*geta);
    for(j=0; j<n; j++){
        pp13[1] += Pi1_1 * Pi0_2 * p12[j]     * exp(lbf1[j]-geta); // eqtl
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); // atac
        pp13[5] += Pi1_2 * p12[j+n*2]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio 
        diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];

    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}

void pwhmnewataceqtl(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    double Pi0_1 = 1.0-Pi1_1;
    double Pi0_2 = 1.0-Pi1_2;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta0, p12,     w, n);
    //wsoftmax(eta1, p12+n,   w, n);
    wsoftmax(eta2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0_1*Pi0_2 * exp(-2.*geta);
    pp13[4] = Pi0_2       * exp(-2.*geta);
    //pp13[4] = (Pi0_1+Pi0_2)/2.0 * exp(-2.*geta);
    //pp13[4] = Pi0_1 * exp(-2.*geta);
    for(j=0; j<n; j++){
      //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1_1 * Pi0_2 * p12[j]     * exp(lbf1[j]-geta); // eqtl
        pp13[2] += Pi0_1 * Pi1_2 * p12[j+n*2] * exp(lbf2[j]-geta); // atac
      //pp13[3] // linkage
        
      //pp13[4]  = Pi0_0;
        pp13[5] += Pi1_2 * p12[j+n*2]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[5] += (Pi1_1+Pi1_2) * (p12[j]+p12[j+n*2])/4.0  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
      //pp13[5] += Pi1_1 * p12[j]  * exp(lbf1[j]+lbf2[j]  -geta*2.); // Pleio
        
        diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];
        
        // vars in the same peak
        /*if(lci[j]>=0){for(k=1; k<j+1; k++){
            //if(lci[j]==106){fprintf(stderr, "%d %d %s %d %s %d\n", j, j-k, rss[j], lci[j], rss[j-k], lci[j-k]);}
            if(lci[j]==lci[j-k]||lci[j-k]<0){
                if(fid2j==157756 &&strcmp("rs558245864",rss[j])==0 && strcmp("rs2409780",rss[j-k])==0){
                    fprintf(stderr, " %d %s %s %lf %lf %lf %lf\n", fid2j, rss[j-k],rss[j], p12[j-k], p12[j+n*2], p12[j], p12[j+n*2-k]);
                }
                //diaglikel += Pi1_1*Pi1_2 *p12[j-k]*p12[j+n*2] * exp(lbf1[j-k]+lbf2[j] - geta*2.);
                //diagprior += p12[j-k]*p12[j+n*2];
                //diaglikel += Pi1_1*Pi1_2 *p12[j]*p12[j+n*2-k] * exp(lbf1[j]+lbf2[j-k] - geta*2.);
                //diagprior += p12[j]*p12[j+n*2-k];
            }else{
                //if(strcmp("rs558245864",rss[j])==0||strcmp("rs558245864",rss[j-k])==0){fprintf(stderr, "  %s %s\n", rss[j-k],rss[j]);}
                break;
            }
        }}*/
    }
    
    //if(fid2j==157756)fprintf(stderr, "%lf %lf\n", diaglikel/(pp13[1]*pp13[2]/Pi0_1/Pi0_2), diagprior);
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0_1/Pi0_2 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}





void pwhmNewAtacGwas(double* lbf1, double* lbf2, double Pi1atac, double Pi1gwas,  double* w, int n, double* pp13){
    int j, k;
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* eta; eta = (double*)calloc(n,   sizeof(double));
    
    double Pi0atac = 1.0-Pi1atac;
    double Pi0gwas = 1.0-Pi1gwas;
    
    // Param 5
    //Pi0atac = Pi1atac = 1.0;
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(eta, p12,     w, n);
    wsoftmax(eta, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0atac * exp(-2.*geta);
    pp13[4] = Pi0atac * exp(-2.*geta);
    for(j=0; j<n; j++){
        //pp13[0]  = Pi0_1 * Pi0_2
        pp13[1] += Pi1atac * Pi0gwas * p12[j]     * exp(lbf1[j]-geta); // atac
        pp13[2] += Pi0atac * Pi1gwas * p12[j+n*2] * exp(lbf2[j]-geta); // gwas
        //pp13[3] // linkage
        
        //pp13[4]  = Pi0_0;
        pp13[5] += Pi1atac * Pi1gwas * p12[j] * exp(lbf1[j]+lbf2[j] - geta*2.); // Pleio
        
        diaglikel += Pi1atac * Pi1gwas * p12[j]*p12[j+n*2] * exp(lbf1[j]+lbf2[j] - geta*2.);
        diagprior += p12[j]*p12[j+n*2];
    }
    
    pp13[3] = (pp13[1]*pp13[2]/Pi0atac/Pi0gwas - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 6);
    
    free(p12);
}







void pwhmfm0(double* lbfj, double* etaj, double Pi1_j, double* w, int n, double* Zj0all){
    int l;
    
    double* p12; p12 = (double*)calloc(n, sizeof(double));
    double* zj; zj = (double*)calloc(n, sizeof(double));
    
    double Pi0_j = 1.0-Pi1_j;
    
    
    double geta = 0;
    for(l=0; l<n; l++){
        if(geta<lbfj[l]){geta=lbfj[l];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(etaj, p12,     w, n);
    
    
    double totj = 0.0;
    // merginal for j
    for(l=0; l<n; l++){
        zj[l]  = Pi1_j * p12[l]   * exp(lbfj[l]-geta);
        totj  += zj[l];
    }
    
    for(l=0; l<n; l++){
        Zj0all[l] = zj[l];
    }
    
    free(p12);
    free(zj);
}

void init_gsl_rand(){
    gsl_rng_env_setup();
    rngT = gsl_rng_default;
    rng = gsl_rng_alloc (rngT);
    gsl_rng_set (rng, (unsigned int)time( NULL )+getpid());
}

double runif(){
    return gsl_ran_flat (rng, 0.0, 1.0);
}

int nk_ran_multinomial1(int n, double* p){
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    double x = runif();
    int i;
    for(i=0; i<n; i++){
        if(x < p[i]){
            return i;
        }else{
            x-=p[i];
        }
    }
    return 0;
}


void nk_ran_multinomial2(int n, double* p, int* res){
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    double x = runif();
    int i;
    double tot=1.0;
    for(i=0; i<n; i++){
        if(x < p[i]){
            res[0] = i;
            tot -= p[i];
            p[i]=0.0;
            break;
        }else{
            x-=p[i];
        }
    }
    // normalize
    for(i=0; i<n; i++){
        p[i] /= tot;
    }
    // another random number
    x = runif();
    for(i=0; i<n; i++){
        if(x < p[i]){
            res[1] = i;
            break;
        }else{
            x-=p[i];
        }
    }
    return;
}


int rnorm(double* x, int n, double a, double b, double s, double* y){
    int i;
    for(i=0; i<n; i++){
        y[i] = a+b*x[i]+gsl_ran_gaussian(rng, s);
    }
}

int rnorm2(double* x1, double* x2, int n, double a, double b1, double b2, double s, double* y){
    int i;
    for(i=0; i<n; i++){
        y[i] = a+b1*x1[i]+b2*x2[i]+gsl_ran_gaussian(rng, s);
    }
}


int rcausal1(double* etaa, double* etaj, double* etak, double* w, int n, int* cvs){
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    wsoftmax(etaa, p12,     w, n);
    wsoftmax(etaj, p12+n,   w, n);
    wsoftmax(etak, p12+n*2, w, n);
    cvs[0] = nk_ran_multinomial1(n, p12);
    cvs[1] = nk_ran_multinomial1(n, p12+n);
    cvs[2] = nk_ran_multinomial1(n, p12+n*2);
}

int rcausal2(double* etaa, double* etaj, double* etak, double* w, int n, int* cvs){
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    wsoftmax(etaa, p12,     w, n);
    wsoftmax(etaj, p12+n,   w, n);
    wsoftmax(etak, p12+n*2, w, n);
    nk_ran_multinomial2(n, p12, cvs);
    nk_ran_multinomial2(n, p12+n, cvs+2);
    nk_ran_multinomial2(n, p12+n*2, cvs+4);
}


void pwhmfm(double* lbfj, double* lbfk, double* lbfmrj, double* lbfmrk, double* etaa, double* etaj, double* etak, double Pi1_j, double Pi1_k, double* w, int n, double* Psi, double* Zjall){
    int l;
    
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* zj; zj = (double*)calloc(n, sizeof(double));
    double* zk; zk = (double*)calloc(n, sizeof(double));
    double* zjk; zjk = (double*)calloc(n, sizeof(double));
    
    double Pi1_a = 0.07862747919984960920381;
    double Pi0_a = 1.0-Pi1_a;
    double Pi0_j = 1.0-Pi1_j;
    double Pi0_k = 1.0-Pi1_k;
    
    
    double geta = 0;
    for(l=0; l<n; l++){
        if(geta<lbfj[l]){geta=lbfj[l];}
        //if(geta<lbfk[l]){geta=lbfk[l];}
        //if(geta<lbfmrj[l]){geta=lbfmrj[l];}
        //if(geta<lbfmrk[l]){geta=lbfmrk[l];}
    }
    geta /= 2.0;
    //geta -= 10.0;
    
    wsoftmax(etaa, p12,     w, n);
    wsoftmax(etaj, p12+n,   w, n);
    wsoftmax(etak, p12+n*2, w, n);
    
    double offdiagprior=1.0;
    for(l=0; l<n; l++){
        offdiagprior -= p12[l+n]*p12[l+n*2];
    }
    
    double totj = 0.0;
    double totk = 0.0;
    // merginal for j and k
    for(l=0; l<n; l++){
        
        totj  += Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        totk  += Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        zj[l]  = Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        zk[l]  = Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        zjk[l]  = Psi[1] * Pi1_a * p12[l]     * exp(lbfj[l]+lbfk[l]  -geta*2.);
        zjk[l] += Psi[2] * Pi1_j * p12[l+n]   * exp(lbfj[l]+lbfmrj[l]-geta*2.);
        zjk[l] += Psi[3] * Pi1_k * p12[l+n*2] * exp(lbfk[l]+lbfmrk[l]-geta*2.);
        
    }
    
    double tmpzj;
    double tot = Psi[0] * Pi0_j*exp(-geta) * (Pi0_k*exp(-geta) + totk)      +      (Psi[1]*Pi0_a + Psi[2]*Pi0_j + Psi[3]*Pi0_k) * exp(-2.*geta);
    for(l=0; l<n; l++){
        tmpzj = zj[l];
        zj[l] = Psi[0] * zj[l]             * (Pi0_k*exp(-geta) + (totk - zk[l])/offdiagprior) + zjk[l];
        zk[l] = Psi[0] * zk[l]             * (Pi0_j*exp(-geta) + (totj - tmpzj)/offdiagprior) + zjk[l];
        
        tot += zj[l];
    }
    
    for(l=0; l<n; l++){
        Zjall[l] += zj[l];
        //zj[l] /= tot;
        //zk[l] /= tot;
    }
    
    free(p12);
    free(zjk);
    free(zj);
    free(zk);
}





void pwhm(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp13, double* pPhi0){
    int j, k;
    double* eta; eta = (double*)calloc(n*3, sizeof(double));
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    //double* bf12; bf12 = (double*)calloc(n, sizeof(double));
    //double* bf1; bf1 = (double*)calloc(n, sizeof(double));
    //double* bf2; bf2 = (double*)calloc(n, sizeof(double));
    //double* bfmr1; bfmr1 = (double*)calloc(n, sizeof(double));
    //double* bfmr2; bfmr2 = (double*)calloc(n, sizeof(double));
    
    double beta_vt[3];
    double beta_loc0[3];
    double beta_loc[3];
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double Phi0 = pPhi0[0];
    double Phi1 = 1.0 - pPhi0[0];
    
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1]; 
    beta_vt[2] = betas[2];
    
    beta_loc0[0] = 0.0;
    beta_loc0[1] = betas[4];
    beta_loc0[2] = betas[4];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[4];
    
    
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    //fprintf(stderr, "geta=%lf\n", geta);
    for(j=0; j<n; j++){
        eta[j] = eta[j+n] = eta[j+n*2] = beta_vt[vt[j]];
        eta[j]     += beta_loc0[loc1[j]];// null
        eta[j+n]   += beta_loc[loc1[j]]; // peak1 is causal
        eta[j+n*2] += beta_loc[loc2[j]]; // peak2 is causal
    }
    wsoftmax(eta,     p12,     w, n);
    wsoftmax(eta+n,   p12+n,   w, n);
    wsoftmax(eta+n*2, p12+n*2, w, n);
    
    double diagprior=0.0;
    double diaglikel=0.0;
    pp13[0] = Pi0*Pi0 * exp(-2.*geta); //*Phi0 
    pp13[4] = Pi0     * exp(-2.*geta); //*Phi1
    for(j=0; j<n; j++){
        pp13[1]  += Pi0 * Pi1 * p12[j+n]   * exp(lbf1[j]-geta); //*Phi0*Pmaster
        pp13[2]  += Pi0 * Pi1 * p12[j+n*2] * exp(lbf2[j]-geta); //*Phi0*Pmaster
      //pp13[3]  // linkage
      //pp13[4]  // null (pleio + causal)
        pp13[5] += Pi1 * p12[j]           * exp(lbf1[j]+lbf2[j]  -geta*2.); //*Phi1*Pdepend 
        pp13[6] += Pi1 * p12[j+n]         * exp(lbf1[j]+lbfmr1[j]-geta*2.); //*Phi1*Pmaster/2 
        pp13[7] += Pi1 * p12[j+n*2]       * exp(lbf2[j]+lbfmr2[j]-geta*2.); //*Phi1*Pmaster/2 
        
        diaglikel += Pi1*Pi1*p12[j+n]*p12[j+n*2]*exp(lbf1[j]+lbf2[j]  -geta*2.);
        diagprior += p12[j+n]*p12[j+n*2];
    }
    
    //pp13[3] = pp13[1]*pp13[2]/Pi0/Pi0; // linkage
    pp13[3] = (pp13[1]*pp13[2]/Pi0/Pi0 - diaglikel)/(1.-diagprior); // linkage
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    
    sumTo1(pp13, 8);
    
    free(eta);
    free(p12);
}





void pwhmOld(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp13, double* pPhi0){
    int j, k;
    double* eta; eta = (double*)calloc(n*3, sizeof(double));
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    //double* bf12; bf12 = (double*)calloc(n, sizeof(double));
    //double* bf1; bf1 = (double*)calloc(n, sizeof(double));
    //double* bf2; bf2 = (double*)calloc(n, sizeof(double));
    //double* bfmr1; bfmr1 = (double*)calloc(n, sizeof(double));
    //double* bfmr2; bfmr2 = (double*)calloc(n, sizeof(double));
    
    double beta_vt[3];
    double beta_loc0[3];
    double beta_loc[3];
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double Phi0 = pPhi0[0];
    double Phi1 = 1.0 - pPhi0[0];
    
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1]; 
    beta_vt[2] = betas[2];
    
    beta_loc0[0] = 0.0;
    beta_loc0[1] = betas[4];
    beta_loc0[2] = betas[4];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[4];
    
    
    
    double geta = 0;
    for(j=0; j<n; j++){
        if(geta<lbf1[j]){geta=lbf1[j];}
        if(geta<lbf2[j]){geta=lbf2[j];}
        if(geta<lbfmr1[j]){geta=lbfmr1[j];}
        if(geta<lbfmr2[j]){geta=lbfmr2[j];}
    }
    geta /= 2.0;
    geta -= 10.0;
    //fprintf(stderr, "geta=%lf\n", geta);
    for(j=0; j<n; j++){
        //bf1[j]   = exp(lbf1[j]-geta);
        //bf2[j]   = exp(lbf2[j]-geta);
        //bfmr1[j]   = exp(lbfmr1[j]-geta);
        //bfmr2[j]   = exp(lbfmr2[j]-geta);
        //bf12[j]  = bf1[j] * bf2[j];
        eta[j] = eta[j+n] = eta[j+n*2] = beta_vt[vt[j]];
        eta[j]     += beta_loc0[loc1[j]];// null
        eta[j+n]   += beta_loc[loc1[j]]; // peak1 is causal
        eta[j+n*2] += beta_loc[loc2[j]]; // peak2 is causal
    }
    wsoftmax(eta,     p12,     w, n);
    wsoftmax(eta+n,   p12+n,   w, n);
    wsoftmax(eta+n*2, p12+n*2, w, n);
    
    double pj, pk; 
    //pp2[1] = pj = pk = Pi0;
    pp13[0] = Pi0*Pi0 * exp(-2.*geta); //*Phi0 
    pp13[9] = Pi0     * exp(-2.*geta);     //*Phi1
    for(j=0; j<n; j++){
        pp13[1]  += Pi0 * Pi1 * p12[j+n]   * exp(lbf1[j]-geta); //*Phi0*Pmaster
        pp13[2]  += Pi0 * Pi1 * p12[j]     * exp(lbf1[j]-geta); //*Phi0*Pdepend
        
        pp13[3]  += Pi0 * Pi1 * p12[j+n*2] * exp(lbf2[j]-geta); //*Phi0*Pmaster
        pp13[4]  += Pi0 * Pi1 * p12[j]     * exp(lbf2[j]-geta); //*Phi0*Pdepend
        
        pp13[10] += Pi1 * p12[j+n]         * exp(lbf1[j]+lbfmr1[j]-geta*2.); //*Phi1*Pmaster/2 
        pp13[11] += Pi1 * p12[j+n*2]       * exp(lbf2[j]+lbfmr2[j]-geta*2.); //*Phi1*Pmaster/2 
        pp13[12] += Pi1 * p12[j]           * exp(lbf1[j]+lbf2[j]  -geta*2.); //*Phi1*Pdepend 
        
        //pj     += Pi1 * (p12[j]+p12[j+n]  ) /2.0 * bf1[j];
        //pk     += Pi1 * (p12[j]+p12[j+n*2]) /2.0 * bf2[j];
        //pp2[1] += Pi1 * (p12[j]+p12[j+n]+p12[j+n*2])/3.0*bf1[j]*bf2[j];//bf12[j];
    }
    //pp2[0] = pj*pk;
 
    pp13[5] = pp13[1]*pp13[3]/Pi0/Pi0; // *Pm*Pm*Phi0
    pp13[6] = pp13[1]*pp13[4]/Pi0/Pi0; // *Pm*Pd*Phi0
    pp13[7] = pp13[2]*pp13[3]/Pi0/Pi0; // *Pm*Pd*Phi0
    pp13[8] = pp13[2]*pp13[4]/Pi0/Pi0; // *Pd*Pd*Phi0
    
    pp13[1] *= exp(-geta);
    pp13[2] *= exp(-geta);
    pp13[3] *= exp(-geta);
    pp13[4] *= exp(-geta);
    
    sumTo1(pp2, 2);
    sumTo1(pp13, 13);
    
    free(eta);
    free(p12);
    //free(bf12);
    //free(bf1);
    //free(bf2);
    //free(bfmr1);
    //free(bfmr2);
}






// eqtl?
void pwhm1(double* lbf1, double* lbf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp5, double* pPhi0){
    int j, k;
    double* bf1; bf1 = (double*)calloc(n, sizeof(double));
    double* bf2; bf2 = (double*)calloc(n, sizeof(double));
    double* eta; eta = (double*)calloc(n, sizeof(double));
    double* p12; p12 = (double*)calloc(n, sizeof(double));
    double* bf12; bf12 = (double*)calloc(n, sizeof(double));
    
    double beta_vt[3];
    double beta_loc[3];
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double Phi0 = pPhi0[0];
    double Phi1 = 1.0 - pPhi0[0];
    
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1]; 
    beta_vt[2] = betas[2];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[3];
    
    for(j=0; j<n; j++){
        bf1[j] = exp(lbf1[j]); 
        bf2[j] = exp(lbf2[j]);
        bf12[j] = bf1[j] * bf2[j];
        eta[j]  = beta_vt[vt[j]];
        eta[j] += beta_loc[loc1[j]];
    }
    wsoftmax(eta,     p12,     w, n);
    
    double pj, pk; 
    pp2[1] = pj = pk = Pi0;
    pp5[0] = Phi0*Pi0*Pi0 + Phi1*Pi0;
    for(j=0; j<n; j++){
        pp5[1] += Phi0 * Pi0 * Pi1 * p12[j] * bf1[j];
        
        pp5[2] += Phi0 * Pi0 * Pi1 * p12[j] * bf2[j];
        
        pp5[3] += Phi1 * Pi1 * p12[j] * bf12[j];
        
        pj     += Pi1 * p12[j] * bf1[j];
        pk     += Pi1 * p12[j] * bf2[j];
        pp2[1] += Pi1 * p12[j] * bf12[j];
    }
    pp2[0] = pj*pk;
 
    pp5[4]  = pp5[1]*pp5[2]/Pi0/Pi0/Phi0;
    
    sumTo1(pp2, 2);
    sumTo1(pp5, 5);
    
    free(eta);
    free(p12);
    free(bf12);
    free(bf1);
    free(bf2);
}



// atac - eqtl
void pwhm12(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss){
    int j, k;
    double* eta; eta = (double*)calloc(n*3, sizeof(double));
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    
    double beta_vt[3];
    double beta_loc[3];
    
    double beta2_vt[3];
    double beta2_loc[3];
    double beta2_loc0[3];
    
    double Pi0 = 1.0 - betas[0];   // eqtl null
    double Pi1 = betas[0];
    
    double Pi0_2 = 1.0 - betas2[0]; // atac null
    double Pi1_2 = betas2[0];
    
    // eqtl
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1];
    beta_vt[2] = betas[2];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[3];
    
    //atac
    beta2_vt[0] = 0.0;
    beta2_vt[1] = betas2[1];
    beta2_vt[2] = betas2[2];
    
    beta2_loc[0] = 0.0;
    beta2_loc[1] = betas2[3];
    beta2_loc[2] = betas2[4];
    
    beta2_loc0[0] = 0.0;
    beta2_loc0[1] = betas2[4];
    beta2_loc0[2] = betas2[4];
    
    double geta = 0.0;
    
    for(j=0; j<n; j++){
        eta[j]     = beta_vt[vt[j]];  // eqtl
        eta[j+n]   = beta2_vt[vt[j]]; // atac dep
        eta[j+n*2] = beta2_vt[vt[j]]; // atac mas
        
        eta[j]     += beta_loc[loc1[j]];   // eqtl
        eta[j+n]   += beta2_loc0[loc2[j]]; // atac depend
        eta[j+n*2] += beta2_loc[loc2[j]];  // atac master
    }
    wsoftmax(eta,     p12,     w, n);
    wsoftmax(eta+n,   p12+n,   w, n);
    wsoftmax(eta+n*2, p12+n*2, w, n);
    
    pp12[0] = Pi0*Pi0_2*exp(-2.*geta);
    pp12[4] = Pi0_2    *exp(-2.*geta);
    pp12[6] = Pi0_2    *exp(-2.*geta);
    pp12[8] = Pi0      *exp(-2.*geta);
    double diagprior = 0.0;
    double diaglikel = 0.0;
    for(j=0; j<n; j++){
        pp12[1] += Pi0_2 * Pi1   * p12[j]      * exp(lbf1[j]-geta); // eqtl no atac
        pp12[2] += Pi0   * Pi1_2 * p12[j+n*2]  * exp(lbf2[j]-geta); // master
      //pp12[3] // linkage
      //pp12[4] // null pleiotropy
        pp12[5]+=  Pi1_2 * p12[j+n]    * exp(lbf1[j]+lbf2[j]  -2.*geta); // depend <-> eqtl
      //pp12[6] // null causal atac->eqtl
        pp12[7] += Pi1_2 * p12[j+n*2]  * exp(lbf2[j]+lbfmr2[j]-2.*geta); // master -> eqtl
      //pp12[8] // null causal eqtl->atac
        pp12[9] += Pi1   * p12[j]      * exp(lbf1[j]+lbfmr1[j]-2.*geta); // eqtl   -> atac
        
        diaglikel += Pi1 * Pi1_2 * p12[j] * p12[j+n*2] * exp(lbf1[j]+lbf2[j]  -2.*geta);
        diagprior += p12[j] * p12[j+n*2];
    }
    
    //pp12[3]  = pp12[1]*pp12[2]/Pi0/Pi0_2; // eqtl master linkage
    pp12[3]  = (pp12[1]*pp12[2]/Pi0/Pi0_2 - diaglikel)/(1.-diagprior); // eqtl master linkage
    
    pp12[1] *= exp(-geta);
    pp12[2] *= exp(-geta);
    
    sumTo1(pp12, 10);
    
    free(eta);
    free(p12);
}
void pwhm12old(double* lbf1, double* lbf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss){
    int j, k;
    double* eta; eta = (double*)calloc(n*3, sizeof(double));
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* bf1; bf1 = (double*)calloc(n, sizeof(double));
    double* bf2; bf2 = (double*)calloc(n, sizeof(double));
    double* bf12; bf12 = (double*)calloc(n, sizeof(double));
    
    double beta_vt[3];
    double beta_loc[3];
    
    double beta2_vt[3];
    double beta2_loc[3];
    double beta2_loc0[3];
    
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double Pi0_2 = 1.0 - betas2[0];
    double Pi1_2 = betas2[0];
    
    double Phi0 = pPhi0[0];
    double Phi1 = 1.0 - pPhi0[0];
    
    // eqtl
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1];
    beta_vt[2] = betas[2];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[3];
    
    //atac
    beta2_vt[0] = 0.0;
    beta2_vt[1] = betas2[1];
    beta2_vt[2] = betas2[2];
    
    beta2_loc[0] = 0.0;
    beta2_loc[1] = betas2[3];
    beta2_loc[2] = betas2[4];
    
    beta2_loc0[0] = 0.0;
    beta2_loc0[1] = betas2[4];
    beta2_loc0[2] = betas2[4];
    
    double coef = 1.0;
    
    for(j=0; j<n; j++){
        bf1[j] = exp(lbf1[j]-0.0);
        bf2[j] = exp(lbf2[j]-0.0);
        bf12[j] = bf1[j] * bf2[j];
        
        eta[j]     = beta_vt[vt[j]]; // eqtl
        eta[j+n]   = beta2_vt[vt[j]]; // atac dep
        eta[j+n*2] = beta2_vt[vt[j]]; // atac mas
        
        //if(lbf2[j]>50){fprintf(stderr, "%s %d %lf %lf\n", rss[j], loc2[j], beta2_loc[loc2[j]], bf12[j]);}
        
        eta[j]     += beta_loc[loc1[j]];   // eqtl
        eta[j+n]   += beta2_loc0[loc2[j]]; // atac depend
        eta[j+n*2] += beta2_loc[loc2[j]];  // atac master
    }
    wsoftmax(eta,     p12,     w, n);
    wsoftmax(eta+n,   p12+n,   w, n);
    wsoftmax(eta+n*2, p12+n*2, w, n);
    
    
    double pj, pk;
    pj = Pi0   / coef;
    pk = Pi0_2 / coef;
    pp2[1] = Pi0_2;
    pp12[0] = (Phi0*Pi0*Pi0_2 + Phi1*Pi0_2)/coef/coef;
    for(j=0; j<n; j++){
        pp12[1] += Phi0 * (Pi0_2/coef) * (Pi1   * p12[j] * bf1[j]); // eqtl no atac
        
        pp12[2] += Phi0 * (Pi0  /coef) * (Pi1_2 * p12[j+n*2] /2.0 * bf2[j]); // master
        pp12[3] += Phi0 * (Pi0  /coef) * (Pi1_2 * p12[j+n]   /2.0 * bf2[j]); // dep
        
        pp12[4] += Phi1 *                (Pi1_2 * p12[j+n*2] /2.0 * bf12[j]); // master
        pp12[5] += Phi1 *                (Pi1_2 * p12[j+n]   /2.0 * bf12[j]); // depend
        
        pj      += Pi1   *  p12[j]                    * bf1[j];
        pk      += Pi1_2 * (p12[j+n*2]+p12[j+n]) /2.0 * bf2[j];
        pp2[1]  += Pi1_2 * (p12[j+n*2]+p12[j+n]) /2.0 * bf12[j];
    }
    pp2[0] = pj*pk;
    
    pp12[6]  = pp12[1]*pp12[2]/(Pi0/coef)/(Pi0_2/coef)/Phi0;
    pp12[7]  = pp12[1]*pp12[3]/(Pi0/coef)/(Pi0_2/coef)/Phi0;
    
    sumTo1(pp2, 2);
    sumTo1(pp12, 8);
    
    
    free(eta);
    free(p12);
    free(bf1);
    free(bf2);
    free(bf12);
}





// atac - gwas
void pwhm13(double* lbf1, double* lbf2, double* lbfmr, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0){
    int j, k;
    double* eta; eta = (double*)calloc(n*2, sizeof(double));
    double* p12; p12 = (double*)calloc(n*2, sizeof(double));
    
    double beta_vt[3];
    double beta_loc0[3];
    double beta_loc[3];
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double del1 = 0.994987;  // prob being master caQTL estimated in multi-peak caQTL model 
    double del0 = 1.-del1;
    
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1]; 
    beta_vt[2] = betas[2];
    
    beta_loc0[0] = 0.0;
    beta_loc0[1] = betas[4];
    beta_loc0[2] = betas[4];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[4];
    
    double geta = 0.0;
    double pgwas = 0.0;  // flat prior for gwas
    
    for(j=0; j<n; j++){
        eta[j] = eta[j+n] = beta_vt[vt[j]];
        eta[j]     += beta_loc0[loc[j]];
        eta[j+n]   += beta_loc[loc[j]];
        pgwas += w[j];
    }
    //fprintf(stderr, "%lf %d\n", pgwas, n);
    wsoftmax(eta,     p12,     w, n); // null
    wsoftmax(eta+n,   p12+n,   w, n); // causal in peak
    pgwas = 1.0/pgwas;
    
    
    pp12[0] = Pi0;
    pp12[6] = 1.0;
    for(j=0; j<n; j++){
        pp12[1] += Pi1 * p12[j+n]     * exp(lbf1[j]) * del1; // master & no gwas
        pp12[2] += Pi1 * p12[j]       * exp(lbf1[j]) * del0; // depend & no gwas
        
        pp12[3] += Pi0 * pgwas * w[j] * exp(lbf2[j]);  // gwas only
        
        pp12[6] += Pi1 * p12[j+n]     * exp(lbf1[j]+lbfmr[j]-2.*geta) * del1; // master causal
        pp12[7] += Pi1 * p12[j]       * exp(lbf1[j]+lbfmr[j]-2.*geta) * del0; // depend causal
        pp12[6] += Pi1 * p12[j+n]     * exp(lbf1[j]+lbf2[j] -2.*geta) * del1; // master pleiotropy
        pp12[7] += Pi1 * p12[j]       * exp(lbf1[j]+lbf2[j] -2.*geta) * del0; // depend pleiotropy
    }
    
    pp12[4] = pp12[1]*pp12[3]/Pi0;  // gwas master linkage
    pp12[5] = pp12[2]*pp12[3]/Pi0;  // gwas depend linkage
    
    sumTo1(pp12, 8);
    
    free(eta);
    free(p12);
}



void pwhm13old(double* lbf1, double* lbf2, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0){
    int j, k;
    double* eta; eta = (double*)calloc(n*2, sizeof(double));
    double* p12; p12 = (double*)calloc(n*2, sizeof(double));
    double* bf12; bf12 = (double*)calloc(n, sizeof(double));
    double* bf1; bf1 = (double*)calloc(n, sizeof(double));
    double* bf2; bf2 = (double*)calloc(n, sizeof(double));
    
    double beta_vt[3];
    double beta_loc0[3];
    double beta_loc[3];
    double Pi0 = 1.0 - betas[0];
    double Pi1 = betas[0];
    
    double Phi0 = pPhi0[0];
    double Phi1 = 1.0 - pPhi0[0];
    
    double Del0 = pDel0[0];
    double Del1 = 1.0-Del0;
    
    beta_vt[0] = 0.0;
    beta_vt[1] = betas[1]; 
    beta_vt[2] = betas[2];
    
    beta_loc0[0] = 0.0;
    beta_loc0[1] = betas[4];
    beta_loc0[2] = betas[4];
    
    beta_loc[0] = 0.0;
    beta_loc[1] = betas[3];
    beta_loc[2] = betas[4];
    
    double pgwas = 0.0;  // flat prior for gwas
    
    for(j=0; j<n; j++){
        bf1[j]   = exp(lbf1[j]);
        bf2[j]   = exp(lbf2[j]);
        bf12[j]  = bf1[j] * bf2[j];
        eta[j] = eta[j+n] = beta_vt[vt[j]];
        eta[j]     += beta_loc0[loc[j]];
        eta[j+n]   += beta_loc[loc[j]];
        pgwas += w[j];
    }
    //fprintf(stderr, "%lf %d\n", pgwas, n);
    wsoftmax(eta,     p12,     w, n); // null
    wsoftmax(eta+n,   p12+n,   w, n); // causal in peak
    pgwas = 1.0/pgwas;
    
    double pj, pk, pjk;
    pj = pk = pjk = 0.0;
    pp12[0] = Phi0*Pi0*Del0 + Phi1*Pi0;
    for(j=0; j<n; j++){
        pp12[1] += Phi0 * Del0 * Pi1 * p12[j+n]  /2.0 * bf1[j]; // causal in peak
        pp12[2] += Phi0 * Del0 * Pi1 * p12[j]    /2.0 * bf1[j]; // null
        
        pp12[3] += Phi0 * Del1 * Pi0 * pgwas * w[j] * bf2[j];  // gwas only
        
        pp12[4] += Phi1 * Pi1 * p12[j+n]   /2.0 * bf12[j];
        pp12[5] += Phi1 * Pi1 * p12[j]     /2.0 * bf12[j];
        
        pj  += Pi1  * (p12[j]+p12[j+n])/2.0 * bf1[j];  // atac    
        pk  +=        pgwas         * w[j]  * bf2[j];  // gwas 
        pjk += Pi1  * (p12[j]+p12[j+n])/2.0 * bf12[j]; // atac-gwas
        
    }
    pp3[0] = Pi0     + pj;
    pp3[1] = Pi0*pk  + pj*pk;
    pp3[2] = Pi0     + pjk;
    
    pp12[6] = pp12[1]*pp12[3]/Pi0/Del0/Phi0;
    pp12[7] = pp12[2]*pp12[3]/Pi0/Del0/Phi0;
    
    sumTo1(pp3, 3);
    sumTo1(pp12, 8);
    
    free(eta);
    free(p12);
    free(bf12);
    free(bf1);
    free(bf2);
}





























