#include "dgpe.h"
int main(int argc, char **argv)
{
  FILE *fp,*treefile,*seqfile,*seqfileo;
  char value[1000],opt,aaRmat[]="lg.dat";
  int no_opt,nsite,i,j;
  int nchar2,nb,k,status,de=1,fso;
  double llikg,alpha,*ratee,q=0.75;

  seqfileo=NULL; seqfile=NULL; treefile=NULL;
  argc--;
  while((opt=get_param(&argc,argv,value)) != 'z'){
    if(opt=='i'){
      seqfile=fopene(value,"r","sequence file not available"); no_opt=0;
    }
    if(opt=='t'){
      treefile=fopene(value,"r","treefile file not available"); no_opt=0;
    }
    if(opt=='o'){
      seqfileo=fopene(value,"w","output seqfile file not available"); no_opt=0;
    }
    if(opt=='q'){
      q=atof(value); no_opt=0;
    }
    if(no_opt) {printf("Option %c not available\n",opt); exit(1);}  
  }
  if(seqfile==NULL){printf("seqfile needed: -i seqfile\n"); exit(1);}
  if(treefile==NULL){printf("treefile needed: -t treefile\n"); exit(1);}

  ratee=dgpe_rate_est(seqfile,treefile,seqfileo,q,&nsite,&llikg,&alpha);

  printf("DGPE log-likelihood %.10e (alpha=%.5f)\n", llikg, alpha);
  fp = fopen("rate_est.dat", "w");
  for(i = 0; i < nsite; i++) fprintf(fp,"%f\n", ratee[i]);
  fclose(fp);
  return(0);
}
