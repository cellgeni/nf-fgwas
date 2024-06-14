void pwhmfm(double* lbfj, double* lbfk, double* lbfmrj, double* lbfmrk, double* etaa, double* etaj, double* etak, double Pi1_j, double Pi1_k, double* w, int n, double* Zjall){
    int l;
    
    double* p12; p12 = (double*)calloc(n*3, sizeof(double));
    double* zj; zj = (double*)calloc(n, sizeof(double));
    double* zk; zk = (double*)calloc(n, sizeof(double));
    double* zjk; zjk = (double*)calloc(n, sizeof(double));
    
    
    double Phi[3] ={0.934350, 0.053918, 0.011732};
    
    double Pi1_a = 0.015649;
    double Pi0_a = 1.0-Pi1_a;
    double Pi0_j = 1.0-Pi1_j;
    double Pi0_k = 1.0-Pi1_k;
    
    
    double geta = 0;
    for(l=0; l<n; l++){
        if(geta<lbfj[l]){geta=lbfj[l];}
        if(geta<lbfk[l]){geta=lbfk[l];}
        if(geta<lbfmrj[l]){geta=lbfmrj[l];}
        if(geta<lbfmrk[l]){geta=lbfmrk[l];}
    }
    geta /= 2.0;
    geta -= 10.0;
    
    wsoftmax(etaa, p12,     w, n);
    wsoftmax(etaj, p12+n,   w, n);
    wsoftmax(etak, p12+n*2, w, n);
    
    double offdiagprior=1.0;
    for(l=0; l<n; l++){
        offdiagprior -= p12[l+n]*p12[l+n*2];
    }
    
    //double tot  = (Phi[0] * Pi0_j * Pi0_k + Phi[1] * Pi0_0 + Phi[2] * (Pi0_j+Pi0_k) / 2.0) * exp(-2.*geta);
    double totj = 0.0;
    double totk = 0.0;
    for(l=0; l<n; l++){
        
        totj  += Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        totk  += Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        zj[l]  = Pi1_j * p12[l+n]   * exp(lbfj[l]-geta);
        zk[l]  = Pi1_k * p12[l+n*2] * exp(lbfk[l]-geta);
        
        
        zjk[l]  = Phi[1]    * Pi1_a * p12[l]     * exp(lbfj[l]+lbfk[l]  -geta*2.);
        zjk[l] += Phi[2]/2. * Pi1_j * p12[l+n]   * exp(lbfj[l]+lbfmrj[l]-geta*2.);
        zjk[l] += Phi[2]/2. * Pi1_k * p12[l+n*2] * exp(lbfk[l]+lbfmrk[l]-geta*2.);
        
    }
    
    double tmpzj;
    double tot = Phi[0] * Pi0_j*exp(-geta) * (Pi0_k*exp(-geta) + totk)      +      (Phi[1]*Pi0_a + Phi[2]*(Pi0_j+Pi0_k)/2.0) * exp(-2.*geta);
    for(l=0; l<n; l++){
        tmpzj = zj[l];
        zj[l] = Phi[0] * zj[l]             * (Pi0_k*exp(-geta) + (totk - zk[l])/offdiagprior) + zjk[l];
        zk[l] = Phi[0] * zk[l]             * (Pi0_j*exp(-geta) + (totj - tmpzj)/offdiagprior) + zjk[l];
        
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

