#include "dist_est_fn.h"
int main(int argc, char **argv)
{
  FILE *fp,*treefile,*seqfile,*ratefile,*Qfile;
  char aaRmat[80],*seqc,(*names)[11];
  int nsite,nrates,ntaxa,endsite,i,*count,*locs,nchar,model,smodel,j;
  int nchar2,nb,*seq,k,*s,status,which=1,de=1;
  double llikg,ub,*rates,*wt,*likmat,llik,alpha,kappa,*fr,*Qk;
  double one=1.0,*utreec,*utreecn,*M,*W,*Lam,bound,p,q,pl,*ratee,neg=-1.0;

  /* control file */
  fp=fopen(argv[--argc],"r");
  dist_est_paramf(fp,&treefile,&Qfile,&ratefile,&seqfile,aaRmat,&nchar,
		  &model,&nrates,&kappa,&ub,&de);
  fclose(fp);
  /* printf("nchar: %i\n",nchar); */
  /* printf("model: %i\n",model); */
  /* printf("nrate: %i\n",nrates); */
  /* printf("ub: %f\n",ub); */

  
  /* sequence */
  fscanf(seqfile,"%i %i", &ntaxa, &nsite);
  /* printf("%i %i\n", ntaxa, nsite); */
  seqc = (char *) malloc((size_t) ntaxa*nsite*sizeof(char));
  seq = (int *) malloc((size_t) ntaxa*nsite*sizeof(int));
  count = (int *) malloc((size_t) nsite*sizeof(int));
  locs = (int *) malloc((size_t) nsite*sizeof(int));
  names = (char (*)[11]) malloc((size_t) ntaxa*sizeof(*names));
  rinterleavef_nolim(seqfile,&ntaxa,&nsite,names,seqc);
  if (nchar==20)
    for(i = 0; i < ntaxa*nsite; i++) seq[i]=l2ip(seqc[i]);
  if (nchar==4)
    for(i = 0; i < ntaxa*nsite; i++) seq[i]=l2i(seqc[i]);
  endsite=patt_freq_locs(ntaxa,nsite,seq,count,locs);
  /* printf("%i %i\n",nsite,ntaxa); */

  /* substitution model */
  smodel=pmodel2smodel(model,aaRmat,nchar);
  /* frequencies */
  fr = (double *) malloc((size_t) nchar*sizeof(double));
  if ((*aaRmat=='\0' && smodel!=0 && smodel!=3) || (*aaRmat!='\0' && model==3)){ 
    char_freq_count(nchar,ntaxa,endsite,seq,count,fr);
    /* for(i = 0; i < nchar; i++) printf("%f\n",fr[i]); */
  }
  if (*aaRmat!='\0' && model==2) *fr=-1;
  /* K for F84, kappa for HKY */
  if(smodel==2 || smodel==8) Qk=&kappa;
  /* Q matrix for GTR */
  if (smodel == 3){
    Qk = (double *) malloc((size_t) nchar*nchar*sizeof(double));
    for(i = 0; i < nchar; i++) 
      for(j = 0; j < nchar; j++) fscanf(Qfile,"%lf",&Qk[i+j*nchar]);
  }
  subst_model_setparam(nchar,smodel,fr,Qk,&W,&Lam);
  /* for(i = 0; i < nchar; i++) printf("%.16e\n",fr[i]); */


  /* rates and initial weights */
  rates = (double *) malloc((size_t) nrates*sizeof(double));
  wt = (double *) malloc((size_t) nrates*sizeof(double));
  if (ratefile==NULL){
    for(i = 0; i < nrates; i++) rates[i] = i*ub/(nrates-1);
  }
  else{
    for(i = 0; i < nrates; i++) fscanf(ratefile, "%lf", &rates[i]);
  }
  for(i = 0; i < nrates; i++) wt[i] = 1.0/nrates;
  /* for(i = 0; i < nrates; i++) fprintf(stdout, "%f %f\n", rates[i], wt[i]); */
      
  /* Obtain likmat */
  likmat = (double *) malloc((size_t) nsite*nrates*sizeof(double));
  utreec=(double *) malloc((size_t) 4*(ntaxa-1)*sizeof(double));
  utreecn=(double *) malloc((size_t) 4*(ntaxa-1)*sizeof(double));
  s=(int *) malloc((size_t) ntaxa*sizeof(int));
  nchar2 = nchar*nchar;
  nb = 2*ntaxa-3;
  M = (double *) malloc((size_t) nchar2*nb*sizeof(double));
  /* printf("endsite: %i\n",endsite); */
  read_treecsnl(treefile,ntaxa,names,utreec,1);
  utreec[2+(ntaxa-2)*4] += utreec[3+(ntaxa-2)*4];
  utreec[3+(ntaxa-2)*4] = 0.0;
  /* pr_utreec(stdout,ntaxa,utreec); */
  /* for(i = 0; i < ntaxa; i++) printf("%s\n",names[i]); */
  memcpy(utreecn,utreec,(size_t) 4*(ntaxa-1)*sizeof(double));
  for(j = 0; j < nrates; j++){
/*     printf("%i\n",j); */
    if(rates[j] <= 0.0){
      for(k = 0; k < endsite; k++)
	likmat[k+j*endsite] = lik_zero_rate_site(nchar,ntaxa,&seq[k],endsite,fr);
    }
    else{
      for(i = 0; i < (ntaxa-1); i++) 
	for(k = 2; k < 4; k++) 
	  utreecn[k+i*4]=rates[j]*utreec[k+i*4];
      /* pr_utreec(stdout,ntaxa,utreecn); */
      transm_rate(smodel,nchar,ntaxa,utreecn,1,&one,neg,fr,W,Lam,0,M,&one,&one);
      for(k = 0; k < endsite; k++){
	for(i = 0; i < ntaxa; i++) s[i] = seq[k+i*endsite];
	/* for(i=0; i<ntaxa; i++) printf("%c",i2lp(s[i])); printf("\n"); */
	
	likmat[k+j*endsite] = lik_deriv_siteb(nchar,ntaxa,s,fr,utreecn,1,&one,&one,
					      M,&one,&one,0);
      }
    }
    /* for(k = 0; k < endsite; k++) printf("%.16e\n",likmat[k+j*endsite]); */
/*     if(rates[j]==1.0){ */
/*       p=0.0; */
/*       for(i = 0; i < endsite; i++) */
/* 	p += count[i]*log(likmat[i+j*endsite]); */
/*       printf("%.10e\n", p); */
/*     } */
  }

  if(de){
    optrates_harwell(endsite, nrates, count, likmat, wt, rates);
    /* optrates_lbfgsb(endsite, nrates, count, likmat, wt); */
    llik = lik_calc(endsite, nrates, count, likmat, wt); 
    printf("DE log-likelihood %.10e\n", llik);
    fp = fopen("DE.dat", "w");
    for(i = 0; i < nrates; i++) fprintf(fp, "%.16e %.16e\n", rates[i], wt[i]);
    fclose(fp);
  }

  llikg = gamma_mle(endsite, nrates, rates, count, likmat, &alpha);
  printf("DGPE log-likelihood %.10e (alpha=%.5f)\n", llikg, alpha);

  /* rate estimation */
  ratee=(double *) malloc(4*endsite*sizeof(double));
  /* DE rate estimates */
  mode_rate_estf(endsite,nrates,rates,wt,likmat,ratee);
  mean_rate_estf(endsite,nrates,rates,wt,likmat,&ratee[endsite]);
  /* gamma rate estimates */
  cdfgam(&which, &pl, &q, &rates[0], &alpha, &alpha, &status, &bound);
  for(i = 1; i < nrates; i++){
    cdfgam(&which, &p, &q, &rates[i], &alpha, &alpha, &status, &bound);
    wt[i-1] = p - pl;
    pl = p;
  }
  wt[nrates-1] = 1 - p;
  mode_rate_estf(endsite,nrates,rates,wt,likmat,&ratee[2*endsite]);
  mean_rate_estf(endsite,nrates,rates,wt,likmat,&ratee[3*endsite]);

  fp = fopen("rate_est.dat", "w");
  for(i = 0; i < nsite; i++){
    for(j = 0; j < 4; j++)
      fprintf(fp,"%f ", ratee[locs[i]+j*endsite]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  return(0);
}
