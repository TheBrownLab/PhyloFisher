#include "mult-data.h"
int *read_seq(FILE *fp, int nchar, int rm_gap, int *ntaxa, int *nsite, 
	      char (**name)[11]){
  int i,*seq;
  char *seqc;
  
  seqc = read_seqc(fp,ntaxa,nsite,name);
  letter(*ntaxa,*nsite,seqc);
  if (rm_gap) rmgap(*ntaxa,nsite,seqc);
  seq = (int *) malloc((size_t) (*ntaxa)*(*nsite)*sizeof(int));
  if(nchar==20)
    for(i = 0; i < (*nsite)*(*ntaxa); i++) seq[i] = l2ip(seqc[i]);
  else
    for(i = 0; i < (*nsite)*(*ntaxa); i++) seq[i] = l2i(seqc[i]);
  free(seqc);
  return(seq);
}
char *read_seqc(FILE *fp, int *ntaxa, int *nsite, char (**name)[11]){
  char *seqc;
  int fso;

  fso=fscanf(fp, "%i %i",ntaxa,nsite);
  if(fso<=0){printf("Error reading sequence file"); exit(1);}
  *name=(char (*)[11]) malloc((size_t) (*ntaxa)*sizeof(**name));
  seqc = (char *) malloc((size_t) (*ntaxa)*(*nsite)*sizeof(char));
  rinterleavef_nolim(fp,ntaxa,nsite,*name,seqc);
  return(seqc);
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
void mult_suff(int ntaxa, int nsite, int *seq, int *nd, int *n, int *ntot,
	       int *count){
  
  int i,j,k,*nt,nchar=20,isdiff,found;
  
  nt=(int *) malloc(nchar*sizeof(int));
  *nd=0;
  for(i=0; i<nsite; i++){
    for(j=0; j<nchar; j++) nt[j]=0;
    for(j=0; j<ntaxa; j++) 
      if(seq[i+j*nsite]<nchar) nt[seq[i+j*nsite]]++;
    found=0;
    for(j=0; j<(*nd); j++){
      isdiff=0;
      for(k=0; k<nchar; k++)
	if(nt[k]!=n[k+j*nchar]){
	  isdiff=1;
	  break;
	}
      if(!isdiff){
	count[j]++;
	found=1;
	break;
      }
    }
    if(!found){
      memcpy(&n[(*nd)*nchar],nt,nchar*sizeof(int));
      ntot[(*nd)]=0;
      for(j=0; j<nchar; j++){
	ntot[(*nd)] += nt[j];
      }
      count[(*nd)]=1;
      (*nd)++;
    }
  }
  
}
FILE *fopene(char *filename, char *rw, char *err){
  FILE *fp;
  if((fp=fopen(filename,rw))==NULL){
    printf("%s\n",err);
    exit(1);
  }
  return(fp);
}
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
