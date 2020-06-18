#include "mult-data.h"
int main(int argc, char **argv)
{
  FILE *seqfile;
  char value[1000],opt,(*names)[11];
  int no_opt,i,j,*n,*ntot,*count,nchar=20,ntaxa,nsite,*seq,nd,single=0;
  
  seqfile=NULL;
  argc--;
  while((opt=get_param(&argc,argv,value)) != 'z'){
    no_opt=1;
    if(opt=='i'){
      seqfile=fopene(value,"r","sequence file not available"); no_opt=0;
    }
    if(opt=='s'){ single=1; no_opt=0; }
    if(no_opt) {printf("Option %c not available\n",opt); exit(1);}  
  }
  if(seqfile==NULL){printf("sequence file needed: -i seqfile\n"); exit(1);}  

  seq=read_seq(seqfile,nchar,0,&ntaxa,&nsite,&names);
  n=(int *) malloc((size_t) nchar*nsite*sizeof(int));
  ntot=(int *) malloc((size_t) nsite*sizeof(int));
  count=(int *) malloc((size_t) nsite*sizeof(int));
  if(!single)
    mult_suff(ntaxa,nsite,seq,&nd,n,ntot,count);
  if(single){
    nd = nsite;
    for(i=0; i<nd*nchar; i++) n[i]=0;
    for(i=0; i<nd; i++){ ntot[i]=0; count[i]=1; }
    for(i=0; i<nd; i++){
      for(j=0; j<ntaxa; j++)
	if(seq[i+j*nsite]<nchar){
	  ntot[i]++;
	  n[seq[i+j*nsite]+i*nchar]++;
	}}
  }
  
  for(i=0; i<nd; i++){
    for(j=0; j<nchar; j++) printf("%2i ",n[j+i*nchar]);
    printf("%2i ",ntot[i]);
    printf("%i\n",count[i]);
  }
}
