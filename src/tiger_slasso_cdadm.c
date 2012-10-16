#include "R.h"
#include "math.h"

void tiger_slasso_cdadm(double *Z, double *omg, double * x, int *dd, int *nn, int * ite_cnt_ext, int * ite_cnt_int, double *lambda, int *nnlambda, int *maxite, int *col_cnz, int *row_idx, double * prec)
{
    int i,j,k,m,ndata,dim,dim1,dim_sq,ileft,iright,junk_a,size_a,size_a_pre,w_idx,rs_idx,nlambda;
    int ite_ext,ite_int1,ite_int2,gap_ext,max_ite,max_ite1,max_ite2, cnz, tmp_m, y_col;
    double gap_int,ilambda,ilambda_tmp,tmp1,tmp2,omg_dif,eps,eps1,eps2,max_dif, tau0, tau1;

    dim = *dd;
    ndata = *nn;
    dim_sq = dim*dim;
    dim1 = dim-1;
    nlambda = *nnlambda;
    double *omg0 = (double*) malloc(dim1*sizeof(double));
    double *omg1 = (double*) malloc(dim1*sizeof(double));
    int *idx_a = (int*) malloc(dim1*sizeof(int)); //active set
    int *idx_i = (int*) malloc(dim1*sizeof(int)); //inactive set
    int *idx_col = (int*) malloc(dim1*sizeof(int)); //column set
    double *S_col = (double*) malloc(dim1*sizeof(double));
    double *Z_omg = (double*) malloc(ndata*sizeof(double));
    double *gamma_col = (double*) malloc(dim1*sizeof(double));
    double *omg_tild = (double*) malloc(dim1*sizeof(double));
    double *SS = (double*) malloc(dim1*dim1*sizeof(double));
    double *Sy = (double*) malloc(dim1*sizeof(double));

    max_ite = *maxite;
    max_ite1 = 1e2;
    max_ite2 = 1e2;
    eps = * prec;
    eps1 = 1e-3;
    eps2 = 1e-3;
    cnz = 0;

    for(i=0; i<dim; i++){
        y_col = 0;
        for(j=0; j<dim1; j++) {
            if(y_col == i) y_col++;
            idx_col[j] = y_col;
            y_col++;
        }
        for(j=0; j<dim1; j++) {
            omg0[j] = 0;
            idx_i[j] = 1;
            Sy[j] = 0;
            for(m=0; m<ndata; m++){
                Sy[j] += Z[i*ndata+m]*Z[idx_col[j]*ndata+m];
            }
            for(m=0; m<dim1; m++){
                SS[m*dim1+j]=0;
                for(k=0; k<ndata; k++){
                    SS[m*dim1+j] += Z[idx_col[j]*ndata+k] * Z[idx_col[m]*ndata+k];
                }
            }
            S_col[j] = SS[j*dim1+j];
        }
        // idx_i[i] = 0;
        size_a = 0;
            
        for(m=0; m<nlambda; m++) {
            ilambda = lambda[m];
            tau0 = 1;
            ite_ext = 0;
            tmp_m = m*dim_sq+i*dim;
            max_dif = 1;
            while(max_dif > eps && ite_ext < max_ite){ // && max_dif > eps
                // update omega
                ilambda_tmp = tau0*ilambda*ndata;
//if(i==0 && m==0)
//printf("ite_ext=%d,tau0=%f \n",ite_ext,tau0);
                for(j=0; j<dim1; j++) {
                    gamma_col[j] = ilambda_tmp/S_col[j];
                }
                gap_ext = 1;
                ite_int1 = 0;
                while(gap_ext !=0 && ite_int1<max_ite1){
                    size_a_pre = size_a;
                    for(j=0; j<dim1; j++){
                        if(idx_i[j] == 1){
                            omg_tild[j] = 0;
                            for(k=0; k<size_a; k++){
                                w_idx = idx_a[k];
                                omg_tild[j] += SS[w_idx*dim1+j]*omg0[w_idx];
                            }
                            omg_tild[j] = (Sy[j]-omg_tild[j]+S_col[j]*omg0[j])/S_col[j];
                            if(fabs(omg_tild[j])<=gamma_col[j]) {
                                omg1[j] = 0;
                            }
                            else{
                                if(omg_tild[j]>gamma_col[j])
                                    omg1[j] = omg_tild[j] - gamma_col[j];
                                else
                                    omg1[j] = omg_tild[j] + gamma_col[j];
                                idx_a[size_a] = j;
                                size_a++;
                                idx_i[j] = 0;
                            }
                            omg0[j] = omg1[j];
                        }
                    }
                    gap_ext = size_a - size_a_pre;
                    gap_int = 1;
                    ite_int2 = 0;
                    while(gap_int>eps2 && ite_int2<max_ite2){
                        tmp1 = 0;
                        tmp2 = 0;
                        for(j=0; j<size_a; j++){
                            w_idx = idx_a[j];
                            omg_tild[w_idx] = 0;
                            for(k=0; k<size_a; k++){
                                rs_idx = idx_a[k];
                                omg_tild[w_idx] += SS[rs_idx*dim1+w_idx]*omg0[rs_idx];
                            }
                            omg_tild[w_idx] = (Sy[w_idx]-omg_tild[w_idx]+S_col[w_idx]*omg0[w_idx])/S_col[w_idx];
                            if (fabs(omg_tild[w_idx]) <= gamma_col[w_idx]) {
                                omg1[w_idx] = 0;
                            }
                            else {
                                if (omg_tild[w_idx]>gamma_col[w_idx])
                                    omg1[w_idx] = omg_tild[w_idx] - gamma_col[w_idx];
                                else
                                    omg1[w_idx] = omg_tild[w_idx] + gamma_col[w_idx];
                                tmp2 = tmp2+fabs(omg1[w_idx]);
                            }
                            omg_dif = omg1[w_idx]-omg0[w_idx];
                            tmp1 = tmp1+fabs(omg_dif);
                            omg0[w_idx] = omg1[w_idx];
                        }
                        gap_int = tmp1/tmp2;
                        ite_int2++;  
                    }
                    ite_cnt_int[m*dim+i] += ite_int2;
                
                    junk_a = 0;
                    for(j=0; j<size_a; j++){
                        w_idx = idx_a[j];
                        if (omg1[w_idx] == 0){
                            junk_a++;
                            idx_i[w_idx] = 1;
                        }
                        else
                            idx_a[j-junk_a] = w_idx;
                    }
                    size_a = size_a - junk_a;
                    ite_int1++;
                }

                tau1 = 0;
                for(j=0; j<ndata; j++){
                    Z_omg[j]=0;
                    for(k=0; k<size_a; k++){
                        w_idx = idx_a[k];
                        Z_omg[j]+=Z[idx_col[w_idx]*ndata+j]*omg1[w_idx];
                    }
                    Z_omg[j] = Z_omg[j]-Z[i*ndata+j];
                    tau1 += Z_omg[j]*Z_omg[j];
                }
                tau1=sqrt(tau1/ndata);
                //if(ite_ext>000 && (ite_ext%20)==0) printf("col=%d, lambda=%f, ite_ext=%d, omg_out1=%f,omg_dif_sum1=%f,omg_sum1=%f \n",i,ilambda,ite_ext,omg_out1,omg_dif_sum1,omg_sum1);
                max_dif = fabs(tau1-tau0)/tau1;
                tau0 = tau1;
                ite_ext++;
            }
            ite_cnt_ext[m*dim+i] = ite_ext;
            
            omg[tmp_m+i] = 1/(tau0*tau0);
            for(j=0;j<size_a;j++){
                w_idx = idx_a[j];
                omg[tmp_m+idx_col[w_idx]] = -omg1[w_idx]*omg[tmp_m+i];
                x[cnz] = omg[tmp_m+idx_col[w_idx]];
                row_idx[cnz] = m*dim+idx_col[w_idx];
                cnz++;
            }
        }
        col_cnz[i+1]=cnz;
    }
    free(omg0);
    free(omg1);
    free(idx_a);
    free(idx_i);
    free(idx_col);
    free(S_col);
    free(gamma_col);
    free(Z_omg);
    free(omg_tild);
    free(SS);
    free(Sy);
}
