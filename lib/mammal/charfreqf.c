#include "charfreq.h"
void rinterleavef_nolim(FILE *infile, const int *numsp, const int *nsite, 
			char name[][11], char *seq)
{
  int lsite = 0, csite = 0, initb = 1;

  while(csite < *nsite){
    lsite = rinterleave_block(infile,*numsp,*nsite,name,&seq[csite],&initb);
    csite += lsite;
  }
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
char ignore_whitespace_seqfile(FILE *infile){
  char c;
  while((c = getc(infile)) == ' ' || c == '\n' || c == '\r' || c == '\t') ;
  if(c == EOF){
    printf("seqfile: EOF before entire sequence read");
    exit(1);
  }
  return(c);
}
void char_freq(int nchar, int numsp, int nsite, int *seq, double *freq)
{
  int i,j,nc=0;

  for(j = 0; j < nchar; j++){
    freq[j] = 0.0;
  }
  for(i = 0; i < nsite*numsp; i++){
    if(seq[i] < nchar){
      for(j = 0; j < nchar; j++){
	if (seq[i] == j) freq[j]++;
      }
      nc++;
    }
  }

  for(j = 0; j < nchar; j++){
    freq[j] /= nc;
  }
}
void letter(int numsp, int nsite, char *seq)
{
  int i, j;
  for(i = 1; i < numsp; i++)
    for(j = 0; j < nsite; j++)
      if (seq[j + i*nsite] == '.') seq[j + i*nsite]=seq[j];
}
int l2i(char c){
  if (c == 'a' || c == 'A') return(0);
  if (c == 'c' || c == 'C') return(1);
  if (c == 'g' || c == 'G') return(2);
  if (c == 't' || c == 'T') return(3);
  if (c == '-') return(5);
  return(6);
}
void rmgap(int numsp, int *nsite, char *seqc){
  int i,j,k,sp,newnsite;

  newnsite = *nsite;
  for(i = *nsite-1; i >= 0; i--)
    for(j = 0; j < numsp; j++)
      if (seqc[i + j*(*nsite)] == '-'){
	for(k = i; k < newnsite-1; k++)
	  for(sp = 0; sp < numsp; sp++)
	    seqc[k + sp*(*nsite)] = seqc[k+1 + sp*(*nsite)];
	newnsite--;
	break;
      }

  for(j = 1; j < numsp; j++)
    for(i = 0; i < newnsite; i++)
      seqc[i + j*newnsite] = seqc[i + j*(*nsite)];

  *nsite = newnsite;
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
