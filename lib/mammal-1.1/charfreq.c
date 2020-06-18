#include "charfreq.h"
int main(int argc, char **argv)
{
  char (*names)[11],*seqc;
  int *seq,ntaxa,nsite,i,nchar,fso;
  double freq[20];

  if(argc<2){
    printf("The number of character states must be indicated:\n charfreq nchar < infile\n");
    exit(0);
  }
  nchar=atoi(argv[--argc]);
  fso=fscanf(stdin, "%i %i",&ntaxa,&nsite); 
  names=(char (*)[11]) malloc((size_t) ntaxa*sizeof(*names));
  seqc = (char *) malloc((size_t) ntaxa*nsite*sizeof(char));
  seq = (int *) malloc((size_t) ntaxa*nsite*sizeof(int));
  rinterleavef_nolim(stdin,&ntaxa,&nsite,names,seqc);
  if(nchar==4) for(i = 0; i < ntaxa*nsite; i++) seq[i]=l2i(seqc[i]);
  if(nchar==20) for(i = 0; i < ntaxa*nsite; i++) seq[i]=l2ip(seqc[i]);
  char_freq(nchar,ntaxa,nsite,seq,freq);
  for(i = 0; i < nchar; i++) printf("%.16e\n",freq[i]);
  return(0);
}
