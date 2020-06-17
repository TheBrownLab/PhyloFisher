#include "mammal-sigma.h"
int main(int argc, char **argv)
{
  int ntaxa,s,t;
  double *Sigma;

  Sigma=estimate_sigma(stdin,&ntaxa);
  for(s=0; s<ntaxa; s++){
    for(t=0; t<ntaxa; t++) printf("%.8e ",Sigma[s+t*ntaxa]);
    printf("\n");
  }
  return(0);
}
