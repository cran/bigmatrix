#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "R.h"

void tiger_clime_ladm(double * Sigma, double * SS, double * omg, double * x, int *dd, int * ite_cnt, double * lambda, int *nnlambda, double * gamma, int * max_ite, double * rho, int *col_cnz, int * row_idx, double * prec)
{
    int i,j,k,m,dim,dim_sq,tmp_size_a,size_a,size_a_pre,w_idx,nlambda;
    int ite_ext,gap_ext,max_ite1,cnz, tmp_m;
    double ilambda,eps1,alp_dif,mu_dif,max_dif,threshold;
    double omg_sum1, omg_dif_sum1, omg_out1;

    dim = *dd;
    dim_sq = dim*dim;
    nlambda = *nnlambda;
    threshold = 1/((*gamma)*(*rho));
    double *omg0 = (double*) malloc(dim*sizeof(double));
    double *omg1 = (double*) malloc(dim*sizeof(double));
    int *idx_a = (int*) malloc(dim*sizeof(int)); //sizes of active sets
    int *idx_i = (int*) malloc(dim*sizeof(int)); //sizes of inactive sets
    double *alp_tild = (double*) malloc(dim*sizeof(double));
    double *mu_grad = (double*) malloc(dim*sizeof(double));
    double *omg_pre = (double*) malloc(dim*sizeof(double));
    double *alp = (double*) malloc(dim*sizeof(double));
    double *mu = (double*) malloc(dim*sizeof(double));
    double *e_i = (double*) malloc(dim*sizeof(double));
    double *S_omg = (double*) malloc(dim*sizeof(double));
    double *omg_tild = (double*) malloc(dim*sizeof(double));
    double *y_i = (double*) malloc(dim*sizeof(double));
    double *Sy = (double*) malloc(dim*sizeof(double));

    max_ite1 = *max_ite;
    eps1 = * prec;
    cnz = 0;

    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++) {
            e_i[j] = 0;
            omg0[j] = 0;
            idx_i[j] = 1;
            alp[j] = 0;
            mu[j] = 0;
        }
        // idx_i[i] = 0;
        e_i[i] = 1;
        size_a = 0;
            
        for(m=0; m<nlambda; m++) {
            for(j=0; j<dim; j++) {
                //alp[j] = 0;
                //mu[j] = 0;
                //omg0[j] = 0;
                //idx_i[j] = 1;
            }
            //size_a = 0;
            ilambda = lambda[m];
            gap_ext = 1;
            ite_ext = 0;
            tmp_m = m*dim_sq+i*dim;
            max_dif = 1;
            while((gap_ext !=0 || max_dif > eps1) && ite_ext < max_ite1){ // && max_dif > eps
            //while(ite_ext < max_ite1){ // && max_dif > eps
                // update alpha
                for(j=0; j<dim; j++){
                    S_omg[j]=0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        S_omg[j]+=Sigma[w_idx*dim+j]*omg0[w_idx];
                    }
                    alp_tild[j]=S_omg[j]+mu[j]-e_i[j];
                }
                alp_dif = 0;
                for(j=0; j<dim; j++){
                    if (alp_tild[j]<=-ilambda){
                        alp[j]=-ilambda;
                        alp_dif = fabs(alp_tild[j]+ilambda)>alp_dif ? fabs(alp_tild[j]+ilambda) : alp_dif;
                    }
                    else {
                        if (alp_tild[j]>=ilambda){
                            alp[j] = ilambda;
                            alp_dif = fabs(alp_tild[j]-ilambda)>alp_dif ? fabs(alp_tild[j]-ilambda) : alp_dif;
                        }
                        else
                            alp[j] = alp_tild[j];
                    }
                }

                // update omega
                size_a_pre = size_a;
                for(j=0; j<dim; j++)
                    y_i[j] = e_i[j]+alp[j]-mu[j];
                for(j=0; j<dim; j++){
                    Sy[j] = 0;
                    omg_pre[j] = omg0[j];
                    for(k=0; k<dim; k++){
                        Sy[j] += Sigma[j*dim+k]*y_i[k];
                    }
                }
                for(j=0; j<dim; j++){
                    if(idx_i[j] == 1) {
                        omg_tild[j] = 0;
                        for(k=0; k<size_a; k++){
                            w_idx = idx_a[k];
                            omg_tild[j] += SS[w_idx*dim+j]*omg0[w_idx];
                        }
                        omg_tild[j] = omg0[j]+(Sy[j]-omg_tild[j])/(*gamma);
                        if(fabs(omg_tild[j])<=threshold) {
                            omg1[j] = 0;
                        }
                        else{
                            if(omg_tild[j]>threshold)
                                omg1[j] = omg_tild[j] - threshold;
                            else
                                omg1[j] = omg_tild[j] + threshold;
                            idx_a[size_a] = j;
                            size_a++;
                            idx_i[j] = 0;
                        }
                        omg0[j] = omg1[j];
                    }
                    else {
                        omg_tild[j] = 0;
                        for(k=0; k<size_a; k++){
                            w_idx = idx_a[k];
                            omg_tild[j] += SS[w_idx*dim+j]*omg0[w_idx];
                        }
                        omg_tild[j] = omg0[j]+(Sy[j]-omg_tild[j])/(*gamma);
                        if(fabs(omg_tild[j])<=threshold) {
                            omg1[j] = 0;
                        }
                        else{
                            if(omg_tild[j]>threshold)
                                omg1[j] = omg_tild[j] - threshold;
                            else
                                omg1[j] = omg_tild[j] + threshold;
                        }
                        omg0[j] = omg1[j];
                    }
                }
                
                tmp_size_a = 0;
                for(j=0; j<dim; j++){
                    if (omg0[j] == 0)
                        idx_i[j] = 1;
                    else {
                        idx_a[tmp_size_a] = j;
                        tmp_size_a++;
                        idx_i[j] = 0;
                    }
                }
                size_a = tmp_size_a;
                gap_ext = size_a - size_a_pre;

                omg_sum1 = 0;
                omg_dif_sum1 = 0;
                for(j=0; j<dim; j++){
                    omg_dif_sum1 += fabs(omg_pre[j] - omg0[j]);
                    omg_sum1 += fabs(omg0[j]);
                    S_omg[j]=0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        S_omg[j]+=Sigma[w_idx*dim+j]*omg0[w_idx];
                    }
                }
                omg_out1 = omg_dif_sum1/omg_sum1;

                // update mu
                mu_dif = 0;
                for(j=0; j<dim; j++){
                    mu_grad[j]=S_omg[j]-alp[j]-e_i[j];
                    mu[j] += mu_grad[j];
                    mu_dif = fabs(mu_grad[j])>mu_dif ? fabs(mu_grad[j]) : mu_dif;
                }
                max_dif = omg_out1>mu_dif ? omg_out1 : mu_dif;
                //max_dif = omg_out1;
                ite_ext++;
            }
            for(j=0;j<size_a;j++){
                w_idx = idx_a[j];
                omg[tmp_m+w_idx] = omg0[w_idx];
                if(w_idx != i) {
                    x[cnz] = omg0[w_idx];
                    row_idx[cnz] = m*dim+w_idx;
                    cnz++;
                }
            }
            ite_cnt[m*dim+i] = ite_ext;
        }
        col_cnz[i+1]=cnz;
    }
    free(omg0);
    free(omg1);
    free(idx_a);
    free(idx_i);
    free(alp_tild);
    free(mu_grad);
    free(omg_pre);
    free(alp);
    free(mu);
    free(e_i);
    free(S_omg);
    free(omg_tild);
    free(y_i);
    free(Sy);
}
