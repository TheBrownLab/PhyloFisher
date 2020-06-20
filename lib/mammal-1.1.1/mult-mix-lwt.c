#include "mult-mix-lwt.h"
int main(int argc, char **argv)
{
  FILE *seqfile,*freqfile,*outfile,*wtfile,*errfile,*lwtfile;
  char value[1000],opt;
  int i,fso,nclass,nclass_var,invar=0,nchar=20,c,j,no_opt,itmax=10000,pr=0;
  int ntaxa,iter,plusF=0;
  double *fr,*w,lnl,lnlo,dp,*lwt;
  double dlnl_tol=1.0e-16,dp_tol=1.0e-6,C=-1;

  freqfile=NULL; outfile=NULL; wtfile=NULL; nclass=-1; errfile=stdout;
  seqfile=NULL; lwtfile=NULL;
  argc--;
  while((opt=get_param(&argc,argv,value)) != 'z'){
    no_opt=1;
    if(opt=='c'){
      nclass=atoi(value); no_opt=0;
    }
    if(opt=='d'){plusF=1; no_opt=0;}
    if(opt=='e'){
      errfile=fopene(value,"w","errfile not available"); no_opt=0;
    }
    if(opt=='f'){
      freqfile=fopene(value,"r","freqfile not available"); no_opt=0;
    }
    /* if(opt=='0'){invar=1; no_opt=0;} */
    if(opt=='o'){
      outfile=fopene(value,"w","outfile not available"); no_opt=0;
    }
    if(opt=='p'){
      pr=1; no_opt=0;
    }
    if(opt=='w'){
      wtfile=fopene(value,"w","wtfile not available"); no_opt=0;
    }
    if(opt=='i'){
      seqfile=fopene(value,"r","sequence file not available"); no_opt=0;
    }
    if(opt=='m'){ itmax=atoi(value); no_opt=0; }
    if(opt=='l'){
      lwtfile=fopene(value,"r","likelihood weight file not available"); no_opt=0;
    }
    if(opt=='C'){ C=atof(value); no_opt=0; }
    if(opt=='D'){ nchar=4; no_opt=0; }
    if(no_opt) {printf("Option %c not available\n",opt); exit(1);}  
  }
  if(freqfile==NULL){printf("freqfile needed: -f freqfile\n"); exit(1);}
  if(outfile==NULL){outfile=fopene("tmp.freq","w","tmp.freq NA\n"); }
  if(wtfile==NULL){wtfile=fopene("tmp.wt","w","tmp.wt NA\n"); }
  if(seqfile==NULL){printf("seqfile needed: -i seqfile\n"); exit(1);}
  if(lwtfile==NULL){
    printf("likelihood weight file needed: -l lwtfile\n"); exit(1);
  }

  fso=fscanf(seqfile,"%i",&ntaxa); rewind(seqfile);
  lwt=(double *) malloc((size_t) ntaxa*sizeof(double)); /* likelihood weights */
  for(i=0; i<ntaxa; i++)
    if((fso=fscanf(lwtfile,"%lf",&lwt[i]))<=0){
      printf("lwtfile should have %i entries\n",ntaxa); exit(1);
    }

  fr=(double *) malloc(nclass*nchar*sizeof(double));  /* start frequencies */
  nclass_var=(invar)?(nclass-1):nclass;
  for(c=0; c<nclass_var; c++) 
    for(j=0; j<nchar; j++){
      fso=fscanf(freqfile,"%lf",&fr[j+c*nchar]);
      if(!fso){ printf("freqfile should have %ix%i entries\n",nchar,nclass_var);
	exit(1);
      }
    }
  if(invar)
    for(j=0; j<nchar; j++) fr[j+c*nchar] = 1.0/nchar;

  w=(double *) malloc(nclass*sizeof(double));
  
  iter=mult_mix_lwt(seqfile,nclass,lwt,C,plusF,itmax,&dlnl_tol,&dp_tol,pr,fr,w,
		    &lnl,&lnlo,&dp);

  if(iter>itmax) iter=itmax;
  fprintf(errfile,"%i %.16e %.3e %.3e\n",iter,lnl,dp,lnl-lnlo);
  /* for(c=0; c<nclass; c++) fprintf(outfile,"%.3e\n",w[c]); */
  for(c=0; c<nclass; c++){
    for(j=0; j<nchar; j++) fprintf(outfile,"%.16e ",fr[j+c*nchar]);
    fprintf(outfile,"\n");
  }
  for(c=0; c<nclass; c++) fprintf(wtfile,"%.16e\n",w[c]); 
  return(0);
}
