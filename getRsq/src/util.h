#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>

extern const gsl_rng_type * rngT;
extern gsl_rng * rng;

void init_gsl_rand();
double runif();

void colocWithGWAS(double* lbfj, double* lbfk, double* eta0, double Pi1_j, double Pi1_k, double* w, int n, double* pp13);

double rk1(double x, double z);
void rk(double* x0, double* x, double* xk, long n, long nk);

double expit(double x);
void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci);
void pwhmnewWithPia(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double Pi1_0, double* w, int n, double* pp2, double* pp13, int* lci);

void expand(int* pos, double* val, int n, int* pos2, double* val2, int n2, double* w);
int mini(int a, int b);
int maxi(int a, int b);
void clearAs0(double* x, int n);
void softmax(double* eta, double* y, int n);
void wsoftmax(double* eta, double* y, double* w, int n);
//void pwhm(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp12, double* pphi);
void pwhm(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp12, double* pPhi0);
void pwhm12(double* bf1, double* bf2, double* bfmr1, double* bfrm2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss);
void pwhm12old(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* betas2, double* pp2, double* pp12, double* pPhi0, char** rss);
void pwhm1(double* bf1, double* bf2, int* vt, int* loc1, int* loc2, double* w, int n, double* betas, double* pp2, double* pp5, double* pPhi0);
void pwhm13(double* lbf1, double* lbf2, double* lbfmr, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0);
void pwhm13old(double* lbf1, double* lbf2, int* vt, int* loc, double* w, int n, double* betas, double* pp3, double* pp12, double* pPhi0, double* pDel0);


void pwhmfm(double* lbfj, double* lbfk, double* lbfmrj, double* lbfmrk, double* etaa, double* etaj, double* etak, double Pi1_j, double Pi1_k, double* w, int n, double* Psi, double* Zjall);
void pwhmfm0(double* lbfj, double* etaj, double Pi1_j, double* w, int n, double* Zj0all);
void pwhmnew(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci);
//void pwhmNewAtacGwas(double* lbf1, double* lbf2, double Pi1atac, double* w, int n, double* pp13);
void pwhmnewataceqtl(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j);

void pwhmnewataceqtlAllParam(double* lbf1, double* lbf2, double* eta0, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp13, int* lci, char** rss, int fid2j);

int nk_ran_multinomial1(int n, double* p);
void nk_ran_multinomial2(int n, double* p, int* res);
int rcausal1(double* etaa, double* etaj, double* etak, double* w, int n, int* cvs);
int rcausal2(double* etaa, double* etaj, double* etak, double* w, int n, int* cvs);
int rnorm(double* x, int n, double a, double b, double s, double* y);
int rnorm2(double* x1, double* x2, int n, double a, double b1, double b2, double s, double* y);
void getIntrPost(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* Psi, double* w, int n, double* pp2, double* pp13, int* lci);

void getLeadVar(double* lbf1, double* lbf2, double* lbfmr1, double* lbfmr2, double* eta0, double* eta1, double* eta2, double Pi1_1, double Pi1_2, double* w, int n, double* pp2, double* pp13, int* lci, double* maxp, int* maxj);
