#include "dgpe.h"
#include "aa_empirical.h"
char get_param(int *k, char **arg, char *value){ 
  char opt='z';
  
  *value='\0';
  if(*k==0) return('z');
  if(arg[(*k)][0]!='-'){ /* not an option, rather a value */
    sprintf(value,"%s",arg[(*k)]);
    *k = (*k)-1;
    if(arg[(*k)][0]=='-') opt=arg[(*k)][1];
  }
  else{ /* is an option */
    opt=arg[(*k)][1];
    if(arg[(*k)][1]!='\0' && arg[(*k)][2]!='\0') sprintf(value,"%s",&arg[(*k)][2]);
  }
  *k = (*k)-1;
  return(opt);
}
void winterleavef(FILE *fp, int numsp, int nsite, char name[][11], char *seq)
{
  int nblock, i, j, k, l, kend; 

  nblock = nsite/60 + 1;
  
  fprintf(fp,"%i %i\n", numsp, nsite); 
  for(l = 0; l < nblock; l++){
    for(i = 0; i < numsp; i++){
      if (l == 0)
	for(j = 0; j < 10; j++) putc(name[i][j], fp);
      else
	for(j = 0; j < 10; j++) putc(' ', fp);
      if (l == nblock - 1)
	kend = nsite % 60;
      else
	kend = 60;
      for(k = 0; k < kend; k++){
	if (k % 10 == 0) putc(' ', fp);
	putc(seq[k + 60*l +  i*nsite], fp);
      }
      putc('\n', fp);
    }
    putc('\n', fp);
  }

}
void dorder(int n, double *x, int *o){
  int i,j,to,io;
  double xm; 

  for(i = 0; i < n; i++) o[i] = i;
  for(i = 0; i < n; i++){
    xm = x[o[i]];
    to = o[i];
    io = i;
    for(j = i+1; j < n; j++)
      if (x[o[j]] < xm){
	to = o[j];
	xm = x[o[j]];
	io = j;
      }
    o[io] = o[i];
    o[i] = to;
  }
}
double *dgpe_rate_est(FILE *seqfile, FILE *treefile, FILE *seqfileo,
		      double quant, int *nsites, double *llikg, double *alpha){
  /* IN: seqfile - FILE pointer to sequence data
   *         NOTE: sequence data assumes names are at most 10 characters long
   *     treefile - FILE pointer to newick tree
   *     seqfileo - If not NULL, FILE pointer for output sequence data
   *     quant - sites with rates >= qth quantile are output to seqfileo. 
   *       Only used if seqfileo is not NULL
   * OUT: nsite - number of sites
   *      llikg - maximized log likelihood
   *      alpha - optimal alpha 
   * RETURN: ratee - ratee[i] gives the rate estimate for the ith site */
  char aaRmat[]="lg.dat",*seqc,(*names)[11],*seqco;
  int nchar2,endsite,fso,*seq,nchar=20,nsite,ntaxa,smodel=9,nrates=101,*s,nb;
  int i,j,k,which=1,*o,start,*count,*locs,status;
  double one=1.0,*fr,*rates,*wt,ub=10.0,*likmat,*utreec,*utreecn,*M,neg=-1.0;
  double *ratee,*ratees,p,pl,q,bound,Qk,*W,*Lam;
  /* FILE *ofile; /\* testing *\/ */
  
  /* sequence */
  fso=fscanf(seqfile,"%i %i", &ntaxa, &nsite);
  *nsites=nsite;
  if(fso<2){
    printf("first two entries of sequence file should be ntaxa, nsite");
    exit(1);
  }
  seqc = (char *) malloc((size_t) ntaxa*nsite*sizeof(char));
  seq = (int *) malloc((size_t) ntaxa*nsite*sizeof(int));
  count = (int *) malloc((size_t) nsite*sizeof(int));
  names = (char (*)[11]) malloc((size_t) ntaxa*sizeof(*names));
  locs = (int *) malloc((size_t) nsite*sizeof(int));
  rinterleavef_nolim(seqfile,&ntaxa,&nsite,names,seqc);
  for(i = 0; i < ntaxa*nsite; i++) seq[i]=l2ip(seqc[i]);
  endsite=patt_freq_locs(ntaxa,nsite,seq,count,locs);

  /* frequencies and substitution model*/
  fr = (double *) malloc((size_t) nchar*sizeof(double));
  char_freq_count(nchar,ntaxa,endsite,seq,count,fr);
  subst_model_setparam(nchar,smodel,fr,&Qk,&W,&Lam);

  /* rates and initial weights */
  rates = (double *) malloc((size_t) nrates*sizeof(double));
  wt = (double *) malloc((size_t) nrates*sizeof(double));
  for(i = 0; i < nrates; i++) rates[i] = i*ub/(nrates-1);
  for(i = 0; i < nrates; i++) wt[i] = 1.0/nrates;
  
  /* Obtain likmat */
  likmat = (double *) malloc((size_t) nsite*nrates*sizeof(double));
  utreec=(double *) malloc((size_t) 4*(ntaxa-1)*sizeof(double));
  utreecn=(double *) malloc((size_t) 4*(ntaxa-1)*sizeof(double));
  s=(int *) malloc((size_t) ntaxa*sizeof(int));
  nchar2 = nchar*nchar;
  nb = 2*ntaxa-3;
  M = (double *) malloc((size_t) nchar2*nb*sizeof(double));
  read_treecsnl(treefile,ntaxa,names,utreec,1);
  utreec[2+(ntaxa-2)*4] += utreec[3+(ntaxa-2)*4];
  utreec[3+(ntaxa-2)*4] = 0.0;
  memcpy(utreecn,utreec,(size_t) 4*(ntaxa-1)*sizeof(double));
  for(j = 0; j < nrates; j++){
    if(rates[j] <= 0.0){
      for(k = 0; k < endsite; k++)
	likmat[k+j*endsite] = lik_zero_rate_site(nchar,ntaxa,&seq[k],endsite,
						 fr);
    }
    else{
      for(i = 0; i < (ntaxa-1); i++) 
	for(k = 2; k < 4; k++) 
	  utreecn[k+i*4]=rates[j]*utreec[k+i*4];
      transm_rate(smodel,nchar,ntaxa,utreecn,1,&one,neg,fr,W,Lam,0,M,&one,&one);
      for(k = 0; k < endsite; k++){
	for(i = 0; i < ntaxa; i++) s[i] = seq[k+i*endsite];
	likmat[k+j*endsite] = lik_deriv_siteb(nchar,ntaxa,s,fr,utreecn,1,&one,
					      &one, M,&one,&one,0);
      }
    }
  }

  /* alpha estimate and lnl */
  ratee=(double *) malloc((size_t) endsite*sizeof(double));
  *llikg = gamma_mle(endsite, nrates, rates, count, likmat, alpha);

  /* rate estimate*/
  for(i=0; i<endsite; i++) ratee[i]=0.0;
  cdfgam(&which, &pl, &q, &rates[0], alpha, alpha, &status, &bound);
  for(i = 1; i < nrates; i++){
    cdfgam(&which, &p, &q, &rates[i], alpha, alpha, &status, &bound);
    wt[i-1] = p - pl;
    pl = p;
  }
  wt[nrates-1] = 1 - p;
  mean_rate_estf(endsite,nrates,rates,wt,likmat,ratee);
  ratees=(double *) malloc((size_t) nsite*sizeof(double));
  for(i=0; i<nsite; i++){
    ratees[i]=ratee[locs[i]];
  }

  /* output sequence data */
  if(seqfileo!=NULL){
    o=(int *) malloc((size_t) nsite*sizeof(int));
    dorder(nsite,ratees,o);
    /* ofile=fopen("tmp.order","w"); */
    /* for(i=0; i<nsite; i++){ */
    /*   fprintf(ofile,"%i\n",o[i]); */
    /* } */
    /* fclose(ofile); */
    start=quant*nsite;
    seqco=(char *) malloc((size_t) (nsite-start)*ntaxa*sizeof(double));
    for(i=start; i<nsite; i++)
      for(j=0; j<ntaxa; j++)
	seqco[i-start+j*(nsite-start)]=seqc[o[i]+j*nsite];
    winterleavef(seqfileo,ntaxa,nsite-start,names,seqco);
    free(o); free(seqco); 
  }

  free(seqc); free(seq); free(names); free(locs); free(fr); free(rates);
  free(wt); free(W); free(Lam); free(likmat); free(utreec); free(utreecn);
  free(s); free(ratee);
  return(ratees);
}
void dist_est_paramf(FILE *fp, FILE **treefile, FILE **Qfile, FILE **ratefile, 
		     FILE **seqfile,  char *aaRmat, int *nchar, int *model, 
		     int *nrate, double *kappa, double *ub, int *de){
  char c,opt[100],value[100];
  int is_option,done=0;
  *treefile=NULL;
  *Qfile=NULL;
  *ratefile=NULL;
  *seqfile=NULL;
  *aaRmat='\0';
  *nchar=4;
  *nrate=101;
  *ub=10.0;
  *de=1; /* perform DE estimation */
  while(!done){
    c = ctl_paraml(fp,opt,value,&is_option);
    if(is_option){
/*       printf("option: %s value: %s\n",opt, value); */
      if(strcmp(opt,"treefile")==0) 
	*treefile=fopene(value,"r","treefile is not available");

      if(strcmp(opt,"nchar")==0) *nchar=atoi(value);
      if(strcmp(opt,"model")==0) *model=atoi(value);
      if(strcmp(opt,"aaRatefile")==0) sprintf(aaRmat,"%s",value);
      if(strcmp(opt,"Qfile")==0) 
	*Qfile=fopene(value,"r","Qfile is not available");
      if(strcmp(opt,"kappa")==0) *kappa=atof(value);
      if(strcmp(opt,"ttratio")==0) *kappa=-atof(value);
      
      if(strcmp(opt,"nrate")==0) *nrate=atoi(value);
      if(strcmp(opt,"ratefile")==0) 
	*ratefile=fopene(value,"r","ratefile is not available");
      if(strcmp(opt,"ub")==0) *ub=atof(value);

      if(strcmp(opt,"seqfile")==0) 
	*seqfile=fopene(value,"r","seqfile is not available");
      if(strcmp(opt,"DEest")==0) *de=atoi(value);
    }
    if(c==EOF) done=1;
  }
}
double lik_zero_rate_site(int nchar, int ntaxa, int *s, int incs, 
			  double *fr){
  int is_const=1,c,i,ncs;
  double l=0.0;
  
  ncs = ntaxa*incs;
  for(i=0; i<ncs; i+=incs)
    if(s[i]<nchar){
      c=s[i];
      break;
    }
  for(; i < ncs; i+=incs)
    if(s[i]<nchar && s[i] != c){
      is_const=0;
      break;
    }
      
  if(is_const) l=fr[c];
  return(l);
}
void numder(double x, double alpha, double *f)
{
#define H 1.0e-6
  double h = H;
  int which = 1, status;
  double p, ph, pnh, q, bound, anh, ah;
  
  anh = alpha-h;
  cdfgam(&which, &pnh, &q, &x, &anh, &anh, &status, &bound);
  cdfgam(&which, &p, &q, &x, &alpha, &alpha, &status, &bound);
  ah = alpha+h;
  cdfgam(&which, &ph, &q, &x, &ah, &ah, &status, &bound);
  f[0] = p;
  f[1] = (ph - pnh)/(2.0*h);
  f[2] = (ph - 2.0*p + pnh)/h;
}
void numderv(int nrates, double *wt, double *wtp, double *wtpp, double *br, 
	    double alpha)
{
  int i,k;
  double fl[3], fu[3];

  numder(br[0], alpha, fl);
  for(i = 1; i < nrates; i++){
    numder(br[i], alpha, fu);
    wt[i-1] = fu[0] - fl[0];
    wtp[i-1] = fu[1] - fl[1];
    wtpp[i-1] = fu[2] - fl[2];
    for(k = 0; k < 3; k++) fl[k] = fu[k];
  }
  wt[nrates-1] = 1-fl[0];
  wtp[nrates-1] = -fl[1];
  wtpp[nrates-1] = -fl[2];
/*   for(i = 1; i < nrates; i++) */
/*     printf("%.6e %.6e %.6e\n", wt[i], wtp[i], wtpp[i]); */
}
double gamma_mle(int endsite, int nrates, double *rates, int *count,
		 double *likmat, double *alpha_max)
{
#define ALPHAL 0.0005
#define ALPHAU 100.0
  int i, j;
  double lgi, lgip, lgipp, f[3], *br, *wt, *wtp, *wtpp;
  double alpha, alphau=ALPHAU, alphal=ALPHAL;

  br = (double *) malloc((size_t) nrates*sizeof(double));
  wt = (double *) malloc((size_t) nrates*sizeof(double));
  wtp = (double *) malloc((size_t) nrates*sizeof(double));
  wtpp = (double *) malloc((size_t) nrates*sizeof(double));
  
  br[0] = 0.0;
  for(j = 1; j < nrates; j++) br[j] = (rates[j] + rates[j-1])/2.0; 
  alpha = 1.0;
  do{
    numderv(nrates, wt, wtp, wtpp, br, alpha);
    f[0] = 0.0; f[1] = 0.0; f[2] = 0.0;
    for(i = 0; i < endsite; i++){
      lgi = 0.0; lgip = 0.0; lgipp = 0.0;
      for(j = 0; j < nrates; j++){
	lgi += wt[j]*likmat[i + j*endsite];
	lgip += wtp[j]*likmat[i + j*endsite];
	lgipp += wtpp[j]*likmat[i + j*endsite];
      }
      f[0] += count[i]*log(lgi);
      f[1] += count[i]*lgip/lgi;
      f[2] += count[i]*(lgipp/lgi - (lgip/lgi)*(lgip/lgi));
    }
  }while(bnewtr(&alphal, &alphau, &alpha, f));
  *alpha_max = alpha;
  printf("%.10e %.10e %.10e\n", alpha, f[0], f[1]); 

  free(br);
  free(wt);
  free(wtp);
  free(wtpp);
  return(f[0]);
}
double lik_calc(int endsite, int nrates, int *count, double *likmat, double *wt)
{
  int i,j;
  double llik, lg;
  llik = 0.0;
  for(i = 0; i < endsite; i++){
    lg = 0.0;
    for(j = 0; j < nrates; j++) lg += wt[j]*likmat[i+j*endsite];
    llik += count[i]*log(lg);
  }

  return(llik);
}
void mean_rate_estf(int endsite, int nrates, double *rate, double *wt, 
		    double *likmat, double *ratee){
  int i,j;
  double lg,joint;
  for(i = 0; i < endsite; i++){
    ratee[i] = 0.0;
    lg = 0.0;
    for(j = 0; j < nrates; j++){
      joint=wt[j] * likmat[i+j*endsite];
      ratee[i] += joint* rate[j];
      lg += joint;
    }
    ratee[i] /= lg;
  }
  return;
}
void makeprotfreqs(double *fr, double *W, double *Lam) 
{
  long i, mineig,nchar=20;
  mineig = 0;
  for (i = 0; i < nchar; i++)
    if (fabs(Lam[i]) < fabs(Lam[mineig]))
      mineig = i;
  for(i = 0; i < nchar; i++) fr[i] =  fabs(W[i+mineig*nchar]);
}
void chQ_freq(int nbin, double *W, double *Lam, double *fr, double *fro)
{
  /* binning/udistfbf.c 
   *     Q = Pi^(-1) W Lam W' 
   * decomposition used to obtain P(t) as
   *     P(t) = Pi^(-1) W exp(Lam t) W'
   * where W is determined from eigen-decomp
   *     Pi^(1/2) Q Pi^(-1/2) = U' Lam U,   W = (Pi)^(1/2) U'
   * To adjust for frequencies Pi*:
   * obtain Q_ij* = Q_ij pi_j* / pi_j
   *        Q_ii* = - sum_{j neq j} Q_ij*
   *        M_ij* =  sqrt(pi_i*) Q_ij* / sqrt(pi_j*)
   *        U and Lam from M = U' Lam U (eigen-decomp)
   *        W = (Pi*)^(1/2) U' */
  int i,j,l;
  double R[400], Offdiag[400], Rt, mu;

/*   for(i = 0; i < nbin; i++) */
/*     for(j = 0; j < nbin; j++){ */
/*       Rt = 0.0; */
/*       for(l = 0; l < nbin; l++) */
/* 	Rt += W[i+l*nbin]*W[j+l*nbin]*Lam[l]; */
/*       Rt /= fro[i]; */
/*       if (i != j) */
/* 	printf("%f\n", Rt); */
/*     } */
/*   exit(0); */

  for(i = 0; i < nbin; i++){ /* R_ij = Q_ij / (pi_i pi_j) */
    for(j = (i+1); j < nbin; j++){
      Rt = 0.0;
      for(l = 0; l < nbin; l++){
	Rt += W[i+l*nbin]*W[j+l*nbin]*Lam[l];
      }
      Rt /= (fro[i]*fro[j]);
      R[i+j*nbin] = Rt; R[j+i*nbin] = Rt;
    }
  }

  for(i = 0; i < nbin; i++){ /* R_ii = Q_ii* */
    Rt = 0.0;
    for(j = 0; j < nbin; j++){
      if (j != i)
	Rt -= R[i+j*nbin]*fr[j];
    }
    R[i+i*nbin] = Rt/fr[i]; /* will cancel with fr[i] below */
  }

/*   for(i = 0; i < nbin; i++) */
/*     for(j = 0; j < nbin; j++) printf("%f\n",R[i+j*nbin]*fr[j]); exit(0); */

  for(i = 0; i < nbin; i++){ /* R_ij = entry of M_ij* above */
    for(j = i; j < nbin; j++){
      R[i+j*nbin] *= sqrt(fr[i]*fr[j]);
      if (j != i) R[j+i*nbin] = R[i+j*nbin];
    }
  }

  eigenRealSym(R, nbin, Lam, Offdiag); /* eigen-decomp of M* */
  for(i = 0; i < nbin; i++){ /* W = (Pi*)^(1/2) U' */
    for(j = 0; j < nbin; j++){
      W[i+j*nbin] = R[j+i*nbin] * sqrt(fr[i]);
    }
  }

  /* rescale eigenvalues so that -sum_i pi_i Q_ii = 1 */
  mu = 0.0;
  for(i = 0; i < nbin; i++)
    for(l = 0; l < nbin; l++) mu -= Lam[l]*W[i+l*nbin]*W[i+l*nbin];
  for(i = 0; i < nbin; i++) Lam[i] /= mu;
}
double int_elamrt(double t, double lam, double alpha, int deriv, double *de){
  double te,dt;
  te = -lam*t/alpha;
  dt= exp(-alpha*log(1+te));
  if(deriv==0) return(dt);
  if(deriv==1 || deriv==2){
    de[0]=lam*exp(-(alpha+1)*log(1+te));
    if(deriv==2) de[1]=(alpha+1)*lam*lam*exp(-(alpha+2)*log(1+te))/alpha;
  }
  if(deriv== 3 || deriv==4){
    de[0] = dt*(te/(1+te)-log(1+te));
    if(deriv==4) 
      de[1]=de[0]*de[0]/dt + 
	dt*(te*te/(alpha*(1+te)*(1+te))-te/(alpha*(1+te)));
  }
  return(dt);
}
void pJCm(double t, int nchar, double *fr, int deriv, 
	    double *M, double *Mp, double *Mpp){
  int i,j;
  double m,emt;
  m = nchar/(nchar - 1.0);
  emt=exp(-t*m);
  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      if(i==j){
	M[i+j*nchar]=1.0/nchar+(1.0-1.0/nchar)*emt;
	if(deriv==1 || deriv==2) Mp[i+j*nchar]=-m*(1.0-1.0/nchar)*emt;
	if(deriv==2) Mpp[i+j*nchar]=m*m*(1.0-1.0/nchar)*emt;
      }
      if(i!=j){
	M[i+j*nchar]=1.0/nchar*(1.0-exp(-m*t));
	if(deriv==1 || deriv==2) Mp[i+j*nchar]=m*1.0/nchar*emt;
	if(deriv==2) Mpp[i+j*nchar]=-m*m*1.0/nchar*emt;
      }
    }
}
void pJCm_gen(double t, int nchar, double *fr, double *W, double *Lam, 
		int deriv, double *M, double *Mp, double *Mpp){
  pJCm(t,nchar,fr,deriv,M,Mp,Mpp);
}
void pJCma(double t, int nchar, double *fr, int deriv, double alpha,
	   double *M, double *Mp, double *Mpp){
  int i,j;
  double mu,emut,demu[2];
  mu = nchar/(nchar - 1.0);
  emut=int_elamrt(t,-1.0*mu,alpha,deriv,demu);

  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      if(i==j){
	M[i+j*nchar]=1.0/nchar+(1.0-1.0/nchar)*emut;
	if(deriv > 0) 
	  Mp[i+j*nchar]=(1.0-1.0/nchar)*demu[0];
	if(deriv==2 || deriv==4) Mpp[i+j*nchar]=(1.0-1.0/nchar)*demu[1];
      }
      if(i!=j){
	M[i+j*nchar]=1.0/nchar-emut/nchar;
	if(deriv > 0) 
	  Mp[i+j*nchar]=-demu[0]/nchar;
	if(deriv==2 || deriv==4) Mpp[i+j*nchar]=-demu[1]/nchar;
      }
    }
}
void pJCma_gen(double t, int nchar, double *fr, double *W, double *Lam, 
	       int deriv, double alpha, double *M, double *Mp, double *Mpp){
  pJCma(t,nchar,fr,deriv,alpha,M,Mp,Mpp);
}
double muF81f(int nchar, double *fr){
  int l;
  double mu = 0.0;
  for(l = 0; l < nchar; l++) mu += fr[l]*(1.0-fr[l]);
  return(1.0/mu);
}
void pF81m(double t, int nchar, double *fr, double mu, int deriv, 
	   double *M, double *Mp, double *Mpp){
  int i,j;
  double emut;
  emut=exp(-mu*t);
  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      if(i==j){
	M[i+j*nchar]=fr[j]+(1.0-fr[j])*emut;
	if(deriv==1 || deriv==2) Mp[i+j*nchar]=-mu*(1.0-fr[j])*emut;
	if(deriv==2) Mpp[i+j*nchar]=mu*mu*(1.0-fr[j])*emut;
      }
      if(i!=j){
	M[i+j*nchar]=fr[j]*(1.0-exp(-mu*t));
	if(deriv==1 || deriv==2) Mp[i+j*nchar]=mu*fr[j]*emut;
	if(deriv==2) Mpp[i+j*nchar]=-mu*mu*fr[j]*emut;
      }
    }
}
void pF81m_gen(double t, int nchar, double *fr, double *mu, double *Lam, 
		 int deriv, double *M, double *Mp, double *Mpp){
  pF81m(t,nchar,fr,*mu,deriv,M,Mp,Mpp);
}
void pF81ma(double t, int nchar, double *fr, double mu, int deriv, double alpha,
	   double *M, double *Mp, double *Mpp){
  int i,j;
  double emut,demu[2];
  emut=int_elamrt(t,-1.0*mu,alpha,deriv,demu);
  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      if(i==j){
	M[i+j*nchar]=fr[j]+(1.0-fr[j])*emut;
	if(deriv > 0) Mp[i+j*nchar]=(1-fr[j])*demu[0];
	if(deriv==2 || deriv==4) Mpp[i+j*nchar]=(1-fr[j])*demu[1];
      }
      if(i!=j){
	M[i+j*nchar]=fr[j]-fr[j]*emut;
	if(deriv > 0) Mp[i+j*nchar]=-fr[j]*demu[0];
	if(deriv==2 || deriv==4) Mpp[i+j*nchar]=-fr[j]*demu[1];
      }
    }
}
void pF81ma_gen(double t, int nchar, double *fr, double *mu, double *Lam, 
		int deriv, double alpha, double *M, double *Mp, double *Mpp){
  pF81ma(t,nchar,fr,*mu,deriv,alpha,M,Mp,Mpp);
}
int is_transitionf(int i, int j){
  return((i==0 && j==2) || (i==2 && j==0) || (i==1 && j==3) || (i==3 && j==1));
}
double muF84f(double K, double *fr){
  double mu = 0.0;
  int l;
/*   printf("%f\n",K); */
  for(l = 0; l < 4; l++) mu += fr[l]*(1.0-fr[l]);
  mu += 2*K*( fr[0]*fr[2]/(fr[0]+fr[2]) + fr[1]*fr[3]/(fr[1]+fr[3]) );
  return(1.0/mu);
}
void pF84m(double t, double mu, double K, double *fr, int deriv, 
	   double *M, double *Mp, double *Mpp){
  int i,j;
  double Pr,Py,Pj,muK,emuKt,emut;
  Pr=fr[0]+fr[2];
  Py=fr[1]+fr[3];
  muK=mu*(K+1); emuKt=exp(-muK*t); emut=exp(-mu*t);

  for(j = 0; j < 4; j++){
    if (j == 0 || j == 2) Pj=Pr;
    if (j == 1 || j == 3) Pj=Py;
    for(i = 0; i < 4; i++){
      if(i==j){
	M[i+j*4]=fr[j]+fr[j]*(1/Pj-1)*emut+(1-fr[j]/Pj)*emuKt;
	if(deriv==1 || deriv==2)
	  Mp[i+j*4]=-mu*fr[j]*(1/Pj-1)*emut-muK*(1-fr[j]/Pj)*emuKt;
	if(deriv==2)
	  Mpp[i+j*4]=mu*mu*fr[j]*(1/Pj-1)*emut + muK*muK*(1-fr[j]/Pj)*emuKt;
      }
      if(i!=j){
	if(is_transitionf(i,j)){
	  M[i+j*4]=fr[j]+fr[j]*(1/Pj-1)*emut-(fr[j]/Pj)*emuKt;
	  if(deriv==1 || deriv==2)
	    Mp[i+j*4]=-mu*fr[j]*(1/Pj-1)*emut+muK*(fr[j]/Pj)*emuKt;
	  if(deriv==2)
	    Mpp[i+j*4]=mu*mu*fr[j]*(1/Pj-1)*emut-muK*muK*(fr[j]/Pj)*emuKt;
	}
	else{
	  M[i+j*4]=fr[j]*(1.0-emut);
	  if(deriv==1 || deriv==2) Mp[i+j*4] = mu*fr[j]*emut;
	  if(deriv==2) Mpp[i+j*4]=-mu*mu*fr[j]*emut;
	}
      }
    }
  }
}
void pF84m_gen(double t, int nchar, double *fr, double *mu, double *K, 
	       int deriv, double *M, double *Mp, double *Mpp){
  pF84m(t,*mu,*K,fr,deriv,M,Mp,Mpp);
}
void pF84ma(double t, double mu, double K, double *fr, int deriv, double alpha,
	    double *M, double *Mp, double *Mpp){
  int i,j;
  double Pr,Py,Pj,emuKt,emut,demuK[2],demu[2];
  Pr=fr[0]+fr[2];
  Py=fr[1]+fr[3];
  emut=int_elamrt(t,-1.0*mu,alpha,deriv,demu);
  emuKt=int_elamrt(t,-1.0*mu*(K+1),alpha,deriv,demuK);
/*   if (deriv==4) printf("%f %f\n",emut,emuKt); */

  for(j = 0; j < 4; j++){
    if (j == 0 || j == 2) Pj=Pr;
    if (j == 1 || j == 3) Pj=Py;
    for(i = 0; i < 4; i++){
      if(i==j){
	M[i+j*4]=fr[j]+fr[j]*(1/Pj-1)*emut+(1-fr[j]/Pj)*emuKt;
	if(deriv > 0)
	  Mp[i+j*4]=fr[j]*(1/Pj-1)*demu[0]+(1-fr[j]/Pj)*demuK[0];
	if(deriv==2 || deriv==4)
	  Mpp[i+j*4]=fr[j]*(1/Pj-1)*demu[1]+(1-fr[j]/Pj)*demuK[1];
      }
      if(i!=j){
	if(is_transitionf(i,j)){
	  M[i+j*4]=fr[j]+fr[j]*(1/Pj-1)*emut-(fr[j]/Pj)*emuKt;
	  if(deriv > 0)
	    Mp[i+j*4]=fr[j]*(1/Pj-1)*demu[0]-(fr[j]/Pj)*demuK[0];
	  if(deriv==2 || deriv==4)
	    Mpp[i+j*4]=fr[j]*(1/Pj-1)*demu[1]-(fr[j]/Pj)*demuK[1];
	}
	else{
	  M[i+j*4]=fr[j]-fr[j]*emut;
	  if(deriv > 0) Mp[i+j*4] = -fr[j]*demu[0];
	  if(deriv==2 || deriv==4) Mpp[i+j*4]=-fr[j]*demu[1];
	}
      }
    }
  }
}
void pF84ma_gen(double t, int nchar, double *fr, double *mu, double *K, 
		int deriv, double alpha, double *M, double *Mp, double *Mpp){
  pF84ma(t,*mu,*K,fr,deriv,alpha,M,Mp,Mpp);
}
double tt2K(double R, double *fr){
  double pY,pR,K;
  pY=fr[1]+fr[3]; pR=fr[0]+fr[2];
  K=(R*pR*pY-fr[0]*fr[2]-fr[1]*fr[3])/(fr[0]*fr[2]/pR+fr[1]*fr[3]/pY);
  return(K);
}
double K2tt(double K, double *fr){
  double pY,pR,R;
  pY=fr[1]+fr[3]; pR=fr[0]+fr[2];
  R=K*(fr[0]*fr[2]/pR+fr[1]*fr[3]/pY)/(pR*pY) + 
    (fr[0]*fr[2]+fr[1]*fr[3])/(pR*pY);
  return(R);
}
double muHKYf(double kappa, double *fr){
  double mu;
  mu=2.0*(kappa*(fr[0]*fr[2]+fr[1]*fr[3])+(fr[0]+fr[2])*(fr[1]+fr[3]));
  return(1.0/mu);
}
void pHKYm(double t, int nchar, double *fr, double mu, double kappa,
	   int deriv, double *M, double *Mp, double *Mpp){
  int i,j;
  double Pr,Py,A,emuAt,emut,muA,Pj;

  Pr=fr[0]+fr[2];
  Py=fr[1]+fr[3];
  emut=exp(-mu*t);

  for(j = 0; j < 4; j++){

    if (j == 0 || j == 2) Pj=Pr;
    if (j == 1 || j == 3) Pj=Py;
    A = 1+Pj*(kappa-1); muA=mu*A;
    emuAt=exp(-muA*t);

    for(i = 0; i < 4; i++){
      if(i==j){
	M[i+j*4] = fr[j]+fr[j]*(1/Pj-1)*emut+(1-fr[j]/Pj)*emuAt;
	if(deriv==1 || deriv==2)
	  Mp[i+j*4] = -mu*fr[j]*(1/Pj-1)*emut-muA*(1-fr[j]/Pj)*emuAt;
	if(deriv==2)
	  Mpp[i+j*4] = mu*mu*fr[j]*(1/Pj-1)*emut+muA*muA*(1-fr[j]/Pj)*emuAt;
      }
      if(i!=j){
	if(is_transitionf(i,j)){
	  M[i+j*4] = fr[j]+fr[j]*(1/Pj-1)*emut-fr[j]*emuAt/Pj;
	  if(deriv==1 || deriv==2)
	    Mp[i+j*4] = -mu*fr[j]*(1/Pj-1)*emut+muA*fr[j]*emuAt/Pj;
	  if(deriv==2)
	    Mpp[i+j*4] = mu*mu*fr[j]*(1/Pj-1)*emut-muA*muA*fr[j]*emuAt/Pj;
	}
	else{
	  M[i+j*4] = fr[j]-fr[j]*emut;
	  if(deriv==1 || deriv==2) Mp[i+j*4] = mu*fr[j]*emut;
	  if(deriv==2) Mpp[i+j*4] = -mu*mu*fr[j]*emut;
	}
      }
    }}
}
void pHKYm_gen(double t, int nchar, double *fr, double *mu, double *kappa, 
	       int deriv, double *M, double *Mp, double *Mpp){
  pHKYm(t,nchar,fr,*mu,*kappa,deriv,M,Mp,Mpp);
}
void pHKYma(double t, int nchar, double *fr, double mu, double kappa,
	    int deriv, double alpha, double *M, double *Mp, double *Mpp){
  int i,j;
  double Pr,Py,A,emuAt,demuA[2],emut,demu[2],muA,Pj;

  Pr=fr[0]+fr[2];
  Py=fr[1]+fr[3];
  emut=int_elamrt(t,-1.0*mu,alpha,deriv,demu);

  for(j = 0; j < 4; j++){

    if (j == 0 || j == 2) Pj=Pr;
    if (j == 1 || j == 3) Pj=Py;
    A = 1+Pj*(kappa-1); muA=mu*A;
    emuAt=int_elamrt(t,-1.0*muA,alpha,deriv,demuA);

    for(i = 0; i < 4; i++){
      if(i==j){
	M[i+j*4] = fr[j]+fr[j]*(1/Pj-1)*emut+(1-fr[j]/Pj)*emuAt;
	if(deriv > 0)
	  Mp[i+j*4] = fr[j]*(1/Pj-1)*demu[0]+(1-fr[j]/Pj)*demuA[0];
	if(deriv==2 || deriv==4)
	  Mpp[i+j*4] = fr[j]*(1/Pj-1)*demu[1]+(1-fr[j]/Pj)*demuA[1];
      }
      if(i!=j){
	if(is_transitionf(i,j)){
	  M[i+j*4] = fr[j]+fr[j]*(1/Pj-1)*emut-fr[j]*emuAt/Pj;
	  if(deriv > 0)
	    Mp[i+j*4] = fr[j]*(1/Pj-1)*demu[0]-fr[j]*demuA[0]/Pj;
	  if(deriv==2 || deriv==4)
	    Mpp[i+j*4] = fr[j]*(1/Pj-1)*demu[1]-fr[j]*demuA[1]/Pj;
	}
	else{
	  M[i+j*4] = fr[j]-fr[j]*emut;
	  if(deriv > 0) Mp[i+j*4] = -fr[j]*demu[0];
	  if(deriv==2 || deriv==4) Mpp[i+j*4] = -fr[j]*demu[1];
	}
      }
    }}
}
void pHKYma_gen(double t, int nchar, double *fr, double *mu, double *kappa, 
	       int deriv, double alpha, double *M, double *Mp, double *Mpp){
  pHKYma(t,nchar,fr,*mu,*kappa,deriv,alpha,M,Mp,Mpp);
}
double tt2kappa(double R, double *fr){
  double pR,pY;
  pR=fr[0]+fr[2];
  pY=fr[1]+fr[3];
  return(R*pR*pY/(fr[0]*fr[2]+fr[1]*fr[3]));
}
double kappa2tt(double kappa,  double *fr){
  double pR,pY;
  pR=fr[0]+fr[2];
  pY=fr[1]+fr[3];
  return(kappa*(fr[0]*fr[2]+fr[1]*fr[3])/(pR*pY));
}
void rmnegpv(int m, double *pv, int incp){
  int i,mi;
  double sum=0.0;
  mi=m*incp;
  for(i=0; i<mi; i+=incp)
    if(pv[i] < 0.0)
      pv[i] = 0.0;
    else
      sum += pv[i];
  for(i=0; i<mi; i+=incp) pv[i] /= sum;
}
void rmnegtrans(int m, double *P){
  int i,j;

  for(i=0; i<m; i++)
    for(j=0; j<m; j++)
      if(P[i+j*m] < 0.0){
	rmnegpv(m,&P[i],m);
	break;
      }
  return;
}
void Q2frWLam(int nbin, double *Q, double *fr, double *W, double *Lam){
  /*     Q = Pi^(-1) W Lam W' 
   *  decomposition used to obtain P(t) as
   *     P(t) = Pi^(-1) W exp(Lam t) W'
   * where W is determined from eigen-decomp
   *     Pi^(1/2) Q Pi^(-1/2) = U' Lam U,   W = (Pi)^(1/2) U' 
   * Assumes Q corresponds to GTR model
   */
  int i,j;
  double R[400], Offdiag[400], mu;

  /* pi_j = [1 + sum_i neq j Q_ji/Q_ij]**(-1) */
  for(j = 0; j < nbin; j++){
    fr[j] = 1.0;
    for(i = 0; i < nbin; i++)
      if (j != i) fr[j] += Q[j+i*nbin]/Q[i+j*nbin];
    fr[j] = 1.0 / fr[j];
  }
  for(i = 0; i < nbin; i++)
    for(j = 0; j < nbin; j++) /* R =  Pi^(1/2) Q Pi^(-1/2) */
      R[i+j*nbin] = sqrt(fr[i])*Q[i+j*nbin]/sqrt(fr[j]);

  /* eigen-decomp of Pi^(1/2) Q Pi^(-1/2) */
  eigenRealSym(R, nbin, Lam, Offdiag); 

  for(i = 0; i < nbin; i++){ /* W = (Pi*)^(1/2) U' */
    for(j = 0; j < nbin; j++)
      W[i+j*nbin] = R[j+i*nbin] * sqrt(fr[i]);
  }

  /* rescale eigenvalues so that -sum_i pi_i Q_ii = 1 */
  mu = 0.0;
  for(i = 0; i < nbin; i++)
    for(j = 0; j < nbin; j++) mu -= Lam[j]*W[i+j*nbin]*W[i+j*nbin];
  for(i = 0; i < nbin; i++) Lam[i] /= mu;
}
void pGTRm(double t, int nchar, double *fr, double *W, double *Lam, 
	   int deriv, double *M, double *Mp, double *Mpp){
  int i,j,k;
  double w;
  for(i = 0; i < nchar*nchar; i++) M[i]=0.0;
  if(deriv==1 || deriv==2) for(i = 0; i < nchar*nchar; i++) Mp[i]=0.0;
  if(deriv==2) for(i = 0; i < nchar*nchar; i++) Mpp[i]=0.0;
  
  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      for(k = 0; k < nchar; k++){
	w=W[i + k*nchar] * W[j+ k*nchar] * exp(Lam[k]*t);
	M[i+j*nchar] += w;
	if(deriv==1 || deriv==2) Mp[i+j*nchar] += Lam[k]*w;
	if(deriv==2) Mpp[i+j*nchar] += Lam[k]*Lam[k]*w;
      }
      M[i+j*nchar] /= fr[i];
      if(deriv==1 || deriv==2) Mp[i+j*nchar] /= fr[i];
      if(deriv==2) Mpp[i+j*nchar] /= fr[i];
    }
  rmnegtrans(nchar,M);
}
void pGTRma(double t, int nchar, double *fr, double *W, double *Lam, 
	    int deriv, double alpha, double *M, double *Mp, double *Mpp){
  int i,j,k;
  double dt,de[2];
  for(i = 0; i < nchar*nchar; i++) M[i]=0.0;
  if(deriv > 0) for(i = 0; i < nchar*nchar; i++) Mp[i]=0.0;
  if(deriv==2 || deriv==4) for(i = 0; i < nchar*nchar; i++) Mpp[i]=0.0;
  
  for(i = 0; i < nchar; i++)
    for(j = 0; j < nchar; j++){
      for(k = 0; k < nchar; k++){
	dt=int_elamrt(t,Lam[k],alpha,deriv,de);
	M[i+j*nchar] += W[i + k*nchar] * W[j+ k*nchar] * dt;
	if(deriv > 0){
	  Mp[i+j*nchar] += W[i + k*nchar] * W[j+ k*nchar] * de[0];
	}
	if(deriv==2 || deriv==4)
	  Mpp[i+j*nchar] += W[i + k*nchar] * W[j+ k*nchar] * de[1];
      }
      M[i+j*nchar] /= fr[i];
      if(deriv > 0) Mp[i+j*nchar] /= fr[i];
      if(deriv==2 || deriv==4) Mpp[i+j*nchar] /= fr[i];
    }
  rmnegtrans(nchar,M);
}
int subst_model_setparam(int nchar, int smodel, double *fr, double *Qk,
			 double **W, double **Lam){
  /* IN:
   * nchar - number of character states
   * smodel - 0:JC, 1:F81, 2:F84, 3:GTR, 4:PAM, 5:JTT, 6:WAG, 7:mtREV24, 8:HKY,
   *          9:LG
   * Qk - rate matrix for GTR, ttratio, K or Kappa for F84, HKY, < 0 => ttratio; 
   * not used otherwise. 
   *
   * IN/OUT:
   * fr - frequencies; reset for JC, determined from Q for GTR
   *      if *fr < 0 for aa models, model-based frequencies used
   * OUT:
   * W,Lam - contain parameters for substitution model
   *
   * RETURN:
   * 2,1 or 0 according to whether both W and Lam were allocated, just W or
   * neither
   */
  int allocWLam=0,aamodel=0,i;
  double *fro,*Wo,*Lamo;
  if (smodel==4 || smodel==5 || smodel==6 || smodel==7 || smodel==9) aamodel=1;

  if (smodel==0){ /* JC */
    for(i = 0; i < nchar; i++) fr[i] = 1.0/nchar;
  }
  if (smodel==1){ /* F81 */
    *W=(double *) malloc((size_t) sizeof(double));
    **W=muF81f(nchar,fr);
/*     for(i = 0; i < nchar; i++) printf("%f\n",fr[i]); printf("\n"); */
/*     printf("%f\n",**W); */
    allocWLam=1;
  }
  if (smodel==2){ /* F84 */
    *W=(double *) malloc((size_t) sizeof(double));
    *Lam=(double *) malloc((size_t) sizeof(double));
/*     printf("%.16e\n",*Qk); */
/*     for(i = 0; i < nchar; i++) printf("%f\n",fr[i]); printf("\n"); */
    **Lam=*Qk;
    if(*Qk < 0) 
      **Lam=tt2K(-(*Qk),fr);
    **W= muF84f(**Lam,fr);
/*     printf("%f %f\n",**W,**Lam); */
    allocWLam=2;
  }
  if (smodel==3){ /* GTR */
    *W = (double *) malloc((size_t) nchar*nchar*sizeof(double));
    *Lam = (double *) malloc((size_t) nchar*sizeof(double));
    Q2frWLam(nchar,Qk,fr,*W,*Lam);
    allocWLam=2;
  }
  if (smodel==4){ /* PAM */
    Wo = pamprobmat; Lamo = pameigmat;
  }
  if (smodel==5){ /* JTT */
    Wo = jttprobmat; Lamo = jtteigmat;
  }
  if (smodel==6){ /* WAG */
    Wo = wagprobmat; Lamo = wageigmat;
  }
  if (smodel==7){ /* mtREV24 */
    Wo = mtREV24probmat; Lamo = mtREV24eigmat;
  }
  if (smodel==8){ /* HKY */
    *W=(double *) malloc((size_t) sizeof(double));
    *Lam=(double *) malloc((size_t) sizeof(double));
    **Lam=*Qk;
    if(*Qk < 0) 
      **Lam=tt2kappa(-(*Qk),fr);
    **W = muHKYf(**Lam,fr);
    allocWLam=2;
  }
  if(smodel==9){ /* LG */
    Wo = lgprobmat; Lamo = lgeigmat;
  }
  if (aamodel){  /* amino acid models */
    if (*fr > 0){ /* different frequencies than for model */
      *W = (double *) malloc((size_t) nchar*nchar*sizeof(double));
      *Lam = (double *) malloc((size_t) nchar*sizeof(double));
      fro = (double *) malloc((size_t) nchar*sizeof(double));
      makeprotfreqs(fro,Wo,Lamo);
      memcpy(*W,Wo,(size_t) nchar*nchar*sizeof(double));
      memcpy(*Lam,Lamo,(size_t) nchar*sizeof(double));
      chQ_freq(nchar,*W,*Lam,fr,fro);
      allocWLam=2;
    }
    else{
      *W=Wo; *Lam=Lamo;
      makeprotfreqs(fr,*W,*Lam);
    }
  }
  return(allocWLam);
}
int subst_model_setpm(int smodel,
		      void (*(*pm))(double t, int nchar, double *fr, 
				    double *W, double *Lam, int deriv, 
				    double *M, double *Mp, 
				    double *Mpp)){
  /* IN:
   * smodel - 0:JC, 1:F81, 2:F84, 3:GTR, 4:PAM, 5:JTT, 6:WAG, 7:mtREV24, 8:HKY,
   *          9:LG
   *
   * OUT:
   * pm - transition matrix function 
   */
  switch(smodel){
  case 0:
    *pm=&pJCm_gen; 
    break;
  case 1:
    *pm=&pF81m_gen; 
    break;
  case 2:
    *pm=&pF84m_gen; ;
    break;
  case 3: case 4: case 5: case 6: case 7: case 9:
    *pm=&pGTRm; 
    break;
  case 8:
    *pm=&pHKYm_gen;
    break;
  }
  return(0);
}
int subst_model_setpma(int smodel,
		       void (*(*pm))(double t, int nchar, double *fr, 
				     double *W, double *Lam, int deriv, 
				     double alpha,
				     double *M, double *Mp, double *Mpp)){
  /* IN:
   * smodel - 0:JC, 1:F81, 2:F84, 3:GTR, 4:PAM, 5:JTT, 6:WAG, 7:mtREV24, 8:HKY,
   *          9:LG
   *
   * OUT:
   * pm - transition matrix function 
   */
  switch(smodel){
  case 0:
    *pm=&pJCma_gen; 
    break;
  case 1:
    *pm=&pF81ma_gen; 
    break;
  case 2:
    *pm=&pF84ma_gen; ;
    break;
  case 3: case 4: case 5: case 6: case 7: case 9:
    *pm=&pGTRma; 
    break;
  case 8:
    *pm=&pHKYma_gen;
    break;
  }
  return(0);
}
void transm_rate(int smodel, int nchar, int ntaxa, double *utreec, 
		 int nrate, double *rate, double alpha,
		 double *fr, double *W, double *Lam, int deriv, 
		 double *M, double *Mp, double *Mpp){
  int nchar2,nb,j,join,kup,r,i,k;
  double tr,one,*Mpt,*Mppt;
  void (*pm)(double t, int nchar, double *fr, double *W, double *Lam, 
	     int deriv, double *M, double *Mp, double *Mpp);
  void (*pma)(double t, int nchar, double *fr, double *W, double *Lam, 
	      int deriv, double alpha, double *M, double *Mp, double *Mpp);
  if (alpha < 0)
    subst_model_setpm(smodel,&pm);
  if (alpha > 0) 
    subst_model_setpma(smodel,&pma);
  nchar2 = nchar*nchar;
  nb = 2*ntaxa-3;
  for(j = 0; j < nrate; j++){
    for(join = 0; join < ntaxa-1; join++){
      kup = (join < ntaxa-2)?2:1;
      for(k = 0; k < kup; k++){
	r = (int) utreec[k+join*4];
	tr = utreec[k+2+join*4]*rate[j];
	Mpt=&one; Mppt=&one;
	if(deriv > 0) Mpt=&Mp[r*nchar2+j*nchar2*nb];
	if(deriv==2||deriv==4) Mppt=&Mpp[r*nchar2+j*nchar2*nb];
	if(alpha < 0)
	  (*pm)(tr,nchar,fr,W,Lam,deriv,&M[r*nchar2+j*nchar2*nb],Mpt,Mppt);
	if(alpha > 0)
	  (*pma)(tr,nchar,fr,W,Lam,deriv,alpha,&M[r*nchar2+j*nchar2*nb],
		 Mpt,Mppt);	
	if(deriv==1 || deriv==2)
	  for(i = 0; i < nchar2; i++) 
	    Mp[i+r*nchar2+j*nchar2*nb] *= rate[j];
	if(deriv==2)
	  for(i = 0; i < nchar2; i++) 
	    Mpp[i+r*nchar2+j*nchar2*nb] *= rate[j]*rate[j];
      }}}
}
int l2i(char c){
  if (c == 'a' || c == 'A') return(0);
  if (c == 'c' || c == 'C') return(1);
  if (c == 'g' || c == 'G') return(2);
  if (c == 't' || c == 'T') return(3);
  if (c == '-') return(5);
  return(6);
}
int l2ip(char c){
  switch(c){
  case 'A':
    return(alanine);
    
  case 'R':
    return(arginine);
    
  case 'N':
    return(asparagine);
    
  case 'D':
    return(aspartic);

  case 'C':
    return(cysteine);

  case 'Q':
    return(glutamine);

  case 'E':
    return(glutamic);

  case 'G':
    return(glycine);

  case 'H':
    return(histidine);

  case 'I':
    return(isoleucine);

  case 'L':
    return(leucine);

  case 'K':
    return(lysine);

  case 'M':
    return(methionine);

  case 'F':
    return(phenylalanine);

  case 'P':
    return(proline);

  case 'S':
    return(serine);

  case 'T':
    return(threonine);

  case 'W':
    return(tryptophan);

  case 'Y':
    return(tyrosine);

  case 'V':
    return(valine);

  case '-':
    return(21);
    
  case 'X':
    return(22);

  case '?':
    return(23);
    /* future: cases X and ? represent unknown aa 
     * B - asparagine or aspartic
     * Z - glutamine or glutamic
     * * - stop codon
     */
  }
  return(100);
}
int rseqf(FILE *fp, char *seq){
  int rsite = 0;
  char c;
  
  while((c = getc(fp)) != '\n' &&  c != EOF){
    if (c != ' ' && c != '\r' && c!= '\t'){
      *seq = c;
      seq++;
      rsite++;
    }
  }

  return(rsite);
}
char ignore_whitespace_seqfile(FILE *infile){
  char c;
  while((c = getc(infile)) == ' ' || c == '\n' || c == '\r' || c == '\t') ;
  if(c == EOF){
    printf("seqfile: EOF before entire sequence read");
    exit(1);
  }
  return(c);
}
int rinterleave_block(FILE *infile, int numsp, int nsite, char name[][11], 
		       char *seq, int *initb){
  char c;
  int i,j,lsite,lsiteo;
  
  /* printf("%i\n",numsp); */
  for(i = 0; i < numsp; i++){
    c = ignore_whitespace_seqfile(infile);
    if((*initb) == 1){
      name[i][0] = c;
      for(j = 1; j < 10; j++)name[i][j] = getc(infile);
      name[i][10] = '\0';
      c = ignore_whitespace_seqfile(infile);      
      /* printf("%s\n",name[i]); */
    }
    seq[i*nsite] = c;
    lsite = rseqf(infile,&seq[1+i*nsite])+1;
    /* printf("%i\n",lsite); */
    if(i==0) lsiteo=lsite;

    if(i > 0 && lsite != lsiteo){
      printf("seqfile: number of character states should be same for all taxa in a block\n");
      exit(1);
    }
  }
  *initb = 0;
  
  return(lsite);
}
void rinterleavef_nolim(FILE *infile, const int *numsp, const int *nsite, 
			char name[][11], char *seq)
{
  int lsite = 0, csite = 0, initb = 1;

  while(csite < *nsite){
    lsite = rinterleave_block(infile,*numsp,*nsite,name,&seq[csite],&initb);
    csite += lsite;
  }
}
void char_freq_count(int nchar, int numsp, int endsite, int *seq, int *count, 
		     double *freq)
{
  int i,j,nc=0,k;

  for(j = 0; j < nchar; j++){
    freq[j] = 0.0;
  }

  for(i = 0; i < endsite; i++){
    for(k = 0; k < numsp; k++)
      if(seq[i+k*endsite] < nchar)
	for(j = 0; j < nchar; j++)
	  if (seq[i+k*endsite] == j){
	    freq[j] += count[i];
	    nc += count[i];
	  }
  }

  for(j = 0; j < nchar; j++) freq[j] /= nc;
}
int patt_freq_locs(int numsp, int nsite, int *seq, int *count, int *locs)
{
  /* locs[i] gives the location of ith site in compressed seq matrix */
  int *s;
  int endsite,i,j,k,match;

  s = (int *) malloc((size_t) nsite*numsp*sizeof(int));
  for(i = 0; i < nsite; i++) count[i]=0;
  endsite = 0;
  for(i = 0; i < nsite; i++){
    /* check for a match */
    match = 0;
    for(k = 0; k < endsite; k++){
      match = 1;
      for(j = 0; j < numsp; j++){
	if (seq[i+j*nsite] != s[j+k*numsp]){
	  match = 0;
	  break;
	}
      }
      if (match){
	count[k] += 1;
	locs[i]=k;
	break;
      }
    }
    
    /* if not a match, new element */
    if (!match){
      for(j = 0; j < numsp; j++) s[j+endsite*numsp] = seq[i+j*nsite];
      count[endsite] = 1;
      locs[i]=endsite;
      endsite++;
    }
  }
  

  for(j = 0; j < numsp; j++)
    for(i = 0; i < endsite; i++)
      seq[i+j*endsite] = s[j+i*numsp];

  free(s);
  return(endsite);
}
int read_tree(FILE *fp, char *tstring)
{
  /* assumes 
   *   '[' and ']' delimit comments only
   *    '][' won't be used to end a comment and start a new one
   *   white space can be removed */
  int ltstring;
  char c;

  ltstring = 0;
  while((c = getc(fp)) != ';' && c!=EOF){ 
    if (c != ' ' && c != '\t' && c != '\r' && c != '\n'){ /* no white space */
      if (c == '[') /* start of comment; ignored */
	while((c = getc(fp)) != ']') ;
      else
	tstring[ltstring++]=c;
    }
  } 
  if(c==EOF){return(-1);}
  tstring[ltstring++]=';';
  tstring[ltstring++]='\0';
  return(ltstring);
}
double get_dlist(char *tstr, int ss, int de, int *ssn, int *ls, int *le)
{
  /* IN:
   * tstr - tree
   * ss - char position to start search
   * de - end position of dlist to be searched
   * OUT:
   * ssn - char position to start new search
   * ls - start pos of label; (ls-1) gives end of dlist found
   * le - end pos of label; 
   * RETURN: edge length (0.1 if not present) */

  char c,els[100]; /* note the limit on edge-length string */
  int lb=0,rb=0,lrb,i,j=0;
  double el=0.1;

  /* printf("%s\n",tstr); */
  /* printf("%i %i\n",ss,de); */
  /* for(i=ss; i<de; i++) putchar(tstr[i]); printf("\n"); */

  i=ss;
  lrb=ss-1;
  while(i <= de-1){
    if ((c=tstr[i]) == '(') lb++;
    if (c == ')'){rb++; lrb=i;}
    if (c == ',' && lb==rb) break;
    i++;
  }
  *ssn=i+1;
  *ls=lrb+1;
  i=lrb+1; 
  while((c=tstr[i]) != ':' && c != ','&& c != ')') i++;
  *le=i-1;
  if(c == ':'){
    while((c=tstr[++i]) != ',' && c != ')') els[j++]=c;
    els[j]='\0';
    el=atof(els);
  }
  return(el);
}
void underscore2blank(char *s){
  while(*s++ != '\0')
    if(*s=='_') *s=' ';
}
void treecsnl(int ntaxa, char *tstr, char names[][11], double *utreec, 
	      int has_names){
  char **ilabels; 
  int j;
  ilabels=treecs(ntaxa, tstr, names, utreec, has_names);
  for(j=0; j<ntaxa-2; j++)free(ilabels[j]);
  free(ilabels);
}
char **treecs(int ntaxa, char *tstr, char names[][11], double *utreec, 
	      int has_names){
  /* IN:
   * tstr - Newick tree w/o comments, whitespace and ; (see read_tree)
   * has_names = 1 if names are available
   * IN/OUT:
   * names - if has_names, used to determine ordering 0,...,ntaxa-1 of 
   *         terminal labels
   *         o/w gives names of the taxa labeled 0,...,ntaxa-1 in utreec
   * OUT:
   * utreec 
   * RETURN:
   * ilabels - internal nodel labels; (ntaxa-2); might be uninitialized
   *    not allocated as a contiguous block
   *
   * COMMENTS:
   * Based on observation that a root subtree is either
   *    (a) a leaf and has no `(' or ')' 
   * or (b) a dlist + label and el, and has equal numbers of `(' and ')' 
   *
   * violations of Newick standard
   * - quoted labels not allowed
   * - '[' and ']' only used to delimit comments
   * - limits of 10 non-null chars on leaf names
   * - underscores are not converted to blanks
   */

  int *ds; /* starting char positions of dlists */
  int *de; /* ending char positions of dlists */
  int cd=0,td=1; /* current and total dlists */
  int cj; /* current join */
  int cil; /* current internal label */
  int ctl=0; /* current terminal label */
  int *ul; /* utreec labels */
  int *ls,*le,*is_leaf,k,i,ss,K,j,ssn,ir,il;
  char **ilabels,tname[11];
  double *el; 

  el = (double *) malloc((size_t) ntaxa*sizeof(double));
  ds = (int *) malloc((size_t) 2*ntaxa*sizeof(int));
  de = (int *) malloc((size_t) 2*ntaxa*sizeof(int));
  ls = (int *) malloc((size_t) ntaxa*sizeof(int));
  le = (int *) malloc((size_t) ntaxa*sizeof(int));
  is_leaf = (int *) malloc((size_t) ntaxa*sizeof(int));
  ul = (int *) malloc((size_t) ntaxa*sizeof(int));
  ilabels = (char **)malloc((size_t) (ntaxa-2)*sizeof(char *));

  /* if(has_names) */
  /*   for(i = 0; i < ntaxa; i++) underscore2blank(names[i]); */

  for(i = 0; i < ntaxa; i++) el[i]=0;
  
  *ds=0; *de=0; *is_leaf=0;
  while(tstr[*de] != '\0')(*de)++;
  while(tstr[*de] != ')') (*de)--; /* root label and edge ignored */

  cj=ntaxa-2;
  cil=2*ntaxa-3;
  while (cj >= 0){/* Main loop */

    /* obtain all root subtrees of the current dlist and update dlists */
    ss=ds[cd]+1;
    K=0;
    while(ss <= de[cd]-1){
      el[K]=get_dlist(tstr,ss,de[cd],&ssn,&ls[K],&le[K]);
      is_leaf[K]=(ls[K] == ss)?1:0;
      if (!is_leaf[K]){
	ds[td]=ss;
	de[td]=ls[K]-1;
	td++;
      }
      ss=ssn;
      K++;
    }

    /* attach ulabels to subtrees */
    for(i = 0; i < K-2; i++){
      ilabels[cil-ntaxa]=(char *)malloc(sizeof(char));
      ilabels[cil-ntaxa][0]='\0';
      cil--;
    }
    for(i=0; i < K; i++){
      if (!is_leaf[i]){
	ilabels[cil-ntaxa]=(char *)malloc((le[i]-ls[i]+2)*sizeof(char));
	memcpy(ilabels[cil-ntaxa],&tstr[ls[i]],(le[i]-ls[i]+1)*sizeof(char));
	ilabels[cil-ntaxa][le[i]+1-ls[i]]='\0';
	/* underscore2blank(ilabels[cil-ntaxa]); */
	ul[i]=cil;
	cil--;
      }
      if (is_leaf[i] && !has_names){
	memcpy(names[ctl],&tstr[ls[i]],(le[i]-ls[i]+1)*sizeof(char));
	for(k=le[i]+1-ls[i]; k < 10; k++) names[ctl][k]=' ';
	names[ctl][10]='\0';
	/* underscore2blank(names[ctl]); */
	ul[i]=ctl;
	ctl++;
      }
      if (is_leaf[i] && has_names){
	memcpy(tname,&tstr[ls[i]],(le[i]-ls[i]+1)*sizeof(char));
	for(k=le[i]+1-ls[i]; k < 10; k++) tname[k]=' ';
	tname[10]='\0';
	/* underscore2blank(tname); */
	for(j = 0; j < ntaxa; j++){
	  if (strncmp(tname,names[j],11)==0) break;
	}
	if(j >= ntaxa){
	  printf("treecs: unable to find %s among taxa\n",tname); exit(0);}
	ul[i]=j;
      } 
    }/* end of attach ulabels ... */
    
    /* printf("current subtree %i\n",ntaxa+cj); */
    /* for(j = ds[cd]; j <= de[cd]; j++) putchar(tstr[j]); printf("\n"); */
    /* printf("number of subtrees: %i\n",K); */
    /* for(i=0; i < K; i++){ */
    /*   printf("%i %f ",ul[i],el[i]); */
    /*   if(!is_leaf[i]){ */
    /* 	for(j = ds[tdo]; j <= de[tdo]; j++) putchar(tstr[j]); */
    /* 	tdo++; */
    /*   } */
    /*   printf("\n"); */
    /* } */
    /* tdo=td; */

    /* adjust labels of existing internal edges for multifurcation */
    if(K>2)
      for(i=cj+1; i<ntaxa-1; i++)
	for(k=0; k<2; k++) 
	  if(utreec[k+i*4]>ntaxa && utreec[k+i*4]<ntaxa+cj) utreec[k+i*4] -= K-2;

    /* add root subtrees to utreec */
    for(i = 0; i < K-2; i++){
      utreec[cj*4]=ul[i]; utreec[1+cj*4]=cj-1+ntaxa; 
      utreec[2+cj*4]=el[i]; utreec[3+cj*4]=0.0;
      cj--;
    }
    if (ul[K-2] < ul[K-1]){
      ir=K-2;il=K-1;
    }
    else{
      ir=K-1;il=K-2;
    }
    utreec[cj*4]=ul[ir]; utreec[2+cj*4]=el[ir];
    utreec[1+cj*4]=ul[il]; utreec[3+cj*4]=el[il];
    cj--;
    cd++;
    
    /* for(i=cj+1; i<ntaxa-1; i++) */
    /*   printf("%i %i %i %f %f\n",i+ntaxa,(int) utreec[i*4],(int) utreec[1+i*4], */
    /* 	     utreec[2+i*4],utreec[3+i*4]); */

  } /* end of Main loop */
  /* pr_utreec(stdout,ntaxa,utreec); */
  /* for(i = 0; i < ntaxa; i++) */
  /*   printf("%s\n",names[i]); */
  
  free(is_leaf); 
  free(ul); 
  free(le); 
  free(ls); 
  free(de); 
  free(ds); 
  free(el); 
  return(ilabels);
}
void read_treecsnl(FILE *treefile, int ntaxa, char names[][11], double *utreec, 
		   int has_names){
  char *tstring;
  /* char tstring[10000]; */
  tstring = (char *) malloc(100000*sizeof(char));
  read_tree(treefile,tstring);
  /* printf("%s\n",tstring); */
  treecsnl(ntaxa,tstring,names,utreec,has_names);
  free(tstring);
}
void backward_probb(int nchar, int r, double *M, double *f, double *bk,
		   double *b){
  int i,j;
  for(i = 0; i < nchar; i++){
    b[i] = 0.0;
    for(j = 0; j < nchar; j++) b[i] += M[i+j*nchar]*f[j]*bk[j];
  }
}
void backward_probsb(int nchar, int ntaxa, double *utreec, double *M, double *f, 
		    double *b){
  /* b[i+(r-ntaxa)*nchar] gives prob, starting from i, of data "backward" 
   * from r, including t_r */ 
  int nchar2,k,r,l,j;
  nchar2 = nchar*nchar;
  r = utreec[(ntaxa-2)*4]; /* root */
  for(j = 0; j < nchar; j++) b[j+(ntaxa-2)*nchar]=1;
  memcpy(&b[(ntaxa-3)*nchar],&f[r*nchar],nchar*sizeof(double));

  for(k = ntaxa-2; k >= 0; k--)
    for(j=0; j<2; j++){
      r = utreec[j+k*4]; l = utreec[1-j+k*4];
      if(r >= ntaxa && r < 2*ntaxa-3)
	backward_probb(nchar,r,&M[r*nchar2],&f[l*nchar],&b[k*nchar],
		       &b[(r-ntaxa)*nchar]);
    }
}
void forward_probsb(int nchar, int ntaxa, int *s, double *utreec, double *M, 
		   double *f){
  /* f[i+r*nchar] gives prob, starting from i, of data "forward" from r, 
   * including t_r */ 
  int nchar2,i,k,r,l,j;
  nchar2 = nchar*nchar;

  /* for(k=0; k<2*ntaxa-3; k++){ */
  /*   printf("%i\n",k); */
  /*   wmat(stdout,nchar,nchar,&M[k*nchar2]); */
  /* } */
  /* exit(0); */
      
  
  for(k = 0; k < ntaxa; k++)
    if(s[k] < nchar) 
      memcpy(&f[k*nchar],&M[s[k]*nchar+k*nchar2],nchar*sizeof(double));
    else
      for(i = 0; i < nchar; i++) f[i+k*nchar] = 1.0;

  for(k = ntaxa; k < 2*ntaxa-3; k++){
    r=utreec[(k-ntaxa)*4]; l=utreec[1+(k-ntaxa)*4];
    for(i = 0; i < nchar; i++){
      f[i+k*nchar] = 0;
      for(j = 0; j < nchar; j++) 
	f[i+k*nchar] += M[i+j*nchar+k*nchar2]*f[j+l*nchar]*f[j+r*nchar];
    }}
  
  r=utreec[(ntaxa-3)*4]; l=utreec[1+(ntaxa-3)*4];
  for(i = 0; i < nchar; i++) 
    f[i+(2*ntaxa-3)*nchar] = f[i+r*nchar]*f[i+l*nchar];
}
double lik_deriv_tedge(int nchar, int s, double *fr, double *bk, double *fl, 
		       double *Mp){
  int i;
  double lpr=0.0;

  if(s > nchar) return(0.0);
  for(i = 0; i < nchar; i++) 
    lpr += fr[i]*bk[i]*fl[i]*Mp[i+s*nchar];
  return(lpr);
}
double lik_deriv_iedge(int nchar, double *fr, double *bk, double *fl, 
		     double *fc1, double *fc2, double *Mp){
  int i,j;
  double oi,lpr=0.0;

  for(i = 0; i < nchar; i++){
    oi = 0.0;
    for(j = 0; j < nchar; j++)
      oi += fc1[j]*fc2[j]*Mp[i+j*nchar];
    lpr += fr[i]*bk[i]*fl[i]*oi;
  }
  return(lpr);
}
void lik_deriv_calcb(int nchar, int ntaxa, int *s, double *fr, double *utreec, 
		     double *Mp, double *f, double *b, double *lp){
  int nchar2,k,j,c1,c2,r,l,jup;

  nchar2 = nchar*nchar;
  for(k = 0; k < ntaxa-1; k++){
    jup = (k < ntaxa-2)?2:1;
    for(j = 0; j < jup; j++){
      r=utreec[j+k*4]; l=utreec[1-j+k*4];
      if(r < ntaxa) 
	lp[r] = lik_deriv_tedge(nchar,s[r],fr,&b[k*nchar],&f[l*nchar],
				&Mp[r*nchar2]);
      else{
	c1 = utreec[(r-ntaxa)*4]; c2 = utreec[1+(r-ntaxa)*4];
	lp[r] = lik_deriv_iedge(nchar,fr,&b[k*nchar],&f[l*nchar],&f[c1*nchar],
				&f[c2*nchar],&Mp[r*nchar2]);
      }
    }}
}
double lik_deriv_site_1rateb(int nchar, int ntaxa, int *s, double *fr, 
			     double *utreec, 
			     double *M, double *Mp, double *lp, int deriv){
  int nchar2,i,r;
  double lik=0.0,*f,*b;
  nchar2 = nchar*nchar;
  f = (double *) malloc((size_t) nchar*(2*ntaxa-2)*sizeof(double));
  forward_probsb(nchar,ntaxa,s,utreec,M,f);
  r = utreec[(ntaxa-2)*4];
  for(i = 0; i < nchar; i++)
    lik += fr[i]*f[i+(2*ntaxa-3)*nchar]*f[i+r*nchar];
  if (deriv){
    b = (double *) malloc((size_t) nchar*(ntaxa-1)*sizeof(double));
    backward_probsb(nchar,ntaxa,utreec,M,f,b);
    lik_deriv_calcb(nchar,ntaxa,s,fr,utreec,Mp,f,b,lp);
    free(b);
  }
  free(f); 
  return(lik);
}
double lik_deriv_siteb(int nchar, int ntaxa, int *s, double *fr, double *utreec, 
		      int nrate, double *rate, double *wt, double *M, 
		      double *Mp, double *lp, int deriv){
  int k,ir,nchar2,nb;
  double lik,*lpt;
  /* Main routine for likelihood calculations at a site
   * Input 
   *   M and MP:  nchar x nchar x nbranch x nrate giving transition probs 
   *              and derivs for all combinations of taxa and edges
   *              MP is not used if deriv = 0
   *   s: The character at the site
   *   rate, wt: the rate distribution
   *   fr: frequencies of character states for site
   *   utreec: the tree
   *   deriv = 1: calculate derivatives
   * Output
   *   lp: derivatives of the likelihood at the site
   * Return: site likelihood
   */  
  nchar2 = nchar*nchar;
  nb = 2*ntaxa-3;
  lik = 0.0;
  lpt = (double *) malloc((size_t) (2*ntaxa-3)*sizeof(double));
  if (deriv)
    for(k = 0; k < 2*ntaxa-3; k++) lp[k] = 0.0;
  for(ir = 0; ir < nrate; ir++){
    lik += wt[ir]*lik_deriv_site_1rateb(nchar,ntaxa,s,fr,utreec,&M[ir*nchar2*nb], 
				       &Mp[ir*nchar2*nb],lpt,deriv);
    if (deriv)
      for(k = 0; k < 2*ntaxa-3; k++) lp[k] += wt[ir]*lpt[k];
  }
  free(lpt);
  return(lik);
}
FILE *fopene(char *filename, char *rw, char *err){
  FILE *fp;
  if((fp=fopen(filename,rw))==NULL){
    printf("%s\n",err);
    exit(1);
  }
  return(fp);
}
char ctl_paraml(FILE *fp, char *opt, char *value, int *is_option){
  /* return last char on line */
  char c;
  int i;

  *is_option=0;
  while(((c = getc(fp)) == ' ' || c == '\t' || c=='\r') && c != EOF && c != '\n');
  if(c == EOF || c == '\n') return(c);
  if(c == '*'){ /* comment */
    while((c = getc(fp)) != '\n' && c != EOF) ;
    return(c);
  } 

  /* option name */
  *is_option=1;
  for(i=0; c !=' ' && c != '\t' && c != '\r'; i++){
    opt[i]=c;
    c=getc(fp);
  }
  opt[i] = '\0';
  /* skip '=' */
  while((c=getc(fp))==' ' || c=='\t' || c == '=') ;
  for(i=0; c != ' ' && c != '\t' && c != '\r' && c != '\n' && c != EOF; i++){
    value[i]=c;
    c=getc(fp);
  }
  value[i]= '\0';
  while(c != '\n' && c != EOF) c=getc(fp);
  return(c);
}
int pmodel2smodel(int model, char *aaRmat, int nchar){
  if (nchar==4){
    if(model==0) return(0); /* JC */
    if(model==2) return(1); /* F81 */
    if(model==3) return(2); /* F84 */
    if(model==4) return(8); /* HKY */
    if(model==7) return(3); /* GTR */
  }
  if (nchar != 4){
    if(model==0) return(0); /* JC */
    if(model==1) return(1); /* F81 */
    if(model==8) return(3); /* GTR */
  }
  if (*aaRmat!='\0'){
    if(strcmp(aaRmat,"dayhoff.dat")==0) return(4); /* PAM */
    if(strcmp(aaRmat,"jones.dat")==0) return(5); /* JTT */
    if(strcmp(aaRmat,"wag.dat")==0) return(6); /* WAG */
    if(strcmp(aaRmat,"mtREV24.dat")==0) return(7); /* mtREV24 */
    if(strcmp(aaRmat,"lg.dat")==0) return(9); /* LG */
  }
  printf("No substitution model available\n");
  exit(1);
}
int bnewtr(double *a, double *b, double *x, double *f)
{
  const int itmax = 100;
/*   const double tol1 = 1.0e-15, tol2 = 1.0e-9; */
  const double tol1 = 1.0e-6, tol2 = 1.0e-6;
  static int iter = 1;
  static double s;

  if ((*b - *a) < tol2){
    iter = 1;
    return(0);
  }
  s = f[1];
  if (fabs(s) < tol1){
    iter = 1;
    return(0);
  }
  if (f[2] < 0.0)  s = -s/f[2];
  if (s > 0.0){
    *a = *x;
    s = (*x + s > *b) ? 0.5*(*b-*a) : s;
  }
  else{
    *b = *x;
    s = (*x + s < *a) ? 0.5*(*a-*b) : s;
  }
  *x += s;
  if (iter == itmax){
    iter = 1;
    return(0);
  }
  iter++;
  return(1);
}
int eigenRealSym(double A[], int n, double Root[], double Offdiag[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return, 
   A has the right vectors and Root has the eigenvalues. work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(), 
   and then using the QL algorithm with implicit shifts.  

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, Offdiag);
   status=EigenTridagQLImplicit(Root, Offdiag, n, A);
   EigenSort(Root, A, n);

   return(status);
}
void EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0;i<n-1;i++) {
      p=d[k=i];
      for (j=i+1;j<n;j++)
         if (d[j] >= p) p=d[k=j];
      if (k != i) {
         d[k]=d[i];
         d[i]=p;
         for (j=0;j<n;j++) {
            p=U[j*n+i];
            U[j*n+i]=U[j*n+k];
            U[j*n+k]=p;
         }
      }
   }
}
void HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix 
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      } 
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}
int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.  
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the 
   tridiagonal matrix, or the output from HouseholderRealSym() to get the 
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+SIGN(r,g));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r=aa*sqrt(1+bb*bb); }
               else if(bb==0)             r=0;
               else            { aa/=bb;  r=bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}
