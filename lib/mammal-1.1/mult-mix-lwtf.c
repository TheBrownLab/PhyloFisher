#include "mult-mix-lwt.h"
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
double xpownb(double x, double n){
  if(n<=0) return(1.0);
  if(x<=0) return(0.0);
  return(exp(n*log(x)));
}
int mult_mix_lwt(FILE *seqfile, int nclass, double *lwt, double C, int plusF,
		 int itmax, double *dlnl_tol, double *dp_tol, int pr,
		 double *fr,
		 double *w, double *lnl, double *lnlo, double *dp){
  /* IN: seqfile - FILE pointer to sequence data
   * NOTE: sequence data assumes names are at most 10 characters long
   *   nclass - the number of classes in the mixture
   *   lwt - double vector of length ntaxa. Gives the likelihood weights.
   *      TODO: Need to update to calculate with a separate routine
   *   C - penalty parameter. C=5 in paper
   *   plusF - plusF=0 => all classes are 
   *     optimized over. plusF=1 => included an additional +F component as the 
   *     last frequency vector in fr (see IN/OUT below). This parameter 
   *     should be set to 0.
   *   itmax, dlnl_tol, dp_tol - convergence criteria 
   *      itmax - the maximum number of iterations. itmax=10000 used in paper
   *      dlnl_tol - exit if log like - log like previous iteration < dlnl_tol.
   *         dlnl_tol=1.0e-16 used in paper
   *      dp_tol - exit if L1 distance between parameters, old-new, < dp_tol.
   *         dp_tol=1.0e-6 used in paper
   *   pr - 1 if the lnl for each iteration should be printed. 0 otherwise
   *IN/OUT: fr - double vector of length 20*nclassi. fr[j+c*20] gives the 
   *         frequency of the jth amino acid for the cth frequency class.
   *         On entry, the starting frequencies
   *         On exit, the estimated frequencies. 
   *OUT: w - double vector of length nclass. The estimated weights for the 
   *       mixture 
   *     lnl, lnlo, dp - used to monitor convergenece
   *       lnl - multinomial log likelihood, on exit
   *       lnlo - multinomial log likelihood, on previous iteration
   *       dp - L1 distance between parameters on exit and at previous iteration
   *RETURN: number of iterations*/
  char (*names)[11];
  int iter,i,x,j,nchar=20,ntaxa,*is_const,*aa,nclass_var,nsite,*seq,invar=0;
  int c,nclass_change;
  double marg,sum,*dn,*dntot,*wo,*post,*fro;
  
  seq=read_seq(seqfile,nchar,0,&ntaxa,&nsite,&names);
  dn=(double *) malloc((size_t) nsite*nchar*sizeof(double));
  dntot=(double *) malloc((size_t) nsite*sizeof(double));
  aa=(int *) malloc((size_t) nsite*sizeof(int));
  is_const=(int *) malloc((size_t) nsite*sizeof(int));
  for(i=0; i<nsite; i++){
    for(x=0; x<nchar; x++) dn[x+i*nchar]=0;
    is_const[i]=1; aa[i]=seq[i]; 
    for(j=0; j<ntaxa; j++)
      if((x=seq[i+j*nsite])<nchar){
	dn[x+i*nchar] += lwt[j];
	if(x!=aa[i]) is_const[i]=0;
      }
    dntot[i]=0.0;
    for(x=0; x<nchar; x++) dntot[i] += dn[x+i*nchar];
  }

  if(invar) nclass++; /* number of classes */
  nclass_var=(invar)?(nclass-1):nclass;
  nclass_change=(plusF)?(nclass_var-1):nclass_var;
  
  for(c=0; c<nclass; c++) w[c] = 1.0/nclass;
  wo=(double *) malloc(nclass*sizeof(double));
  post=(double *) malloc(nclass*sizeof(double));
  for(c=0; c<nclass; c++) post[c]=0.0;
  fro=(double *) malloc(nclass*nchar*sizeof(double));

  for(iter=1; iter<=itmax; iter++){

    memcpy(fro,fr,nchar*nclass*sizeof(double));
    memcpy(wo,w,nclass*sizeof(double));
    if(iter>1) *lnlo=*lnl;

    for(c=0; c<nclass; c++) w[c]=0.0;
    for(i=0; i<nclass_change*nchar; i++) fr[i]=0.0;
    if(invar) for(i=0; i<nchar; i++) fr[i+(nclass-1)*nchar]=0.0;

    *lnl=0.0;
    if(C>0) /* penalty adjustment */
      for(j=0; j<nclass_change*nchar; j++) *lnl += C*log(fro[j]);
    
    for(i=0; i<nsite; i++){

      if(!is_const[i]){
	for(marg=0.0, c=0; c<nclass_var; c++){
	  for(post[c]=wo[c], j=0; j<nchar; j++)
	    post[c] *= xpownb(fro[j+c*nchar],dn[j+i*nchar]);
	  marg += post[c];
	}
	if(invar) post[nclass-1]=0.0;
      }
      
      if(is_const[i]){
	for(marg=0.0, c=0; c<nclass_var; c++){
	  post[c] = wo[c]*xpownb(fro[aa[i]+c*nchar],dntot[i]);
	  marg += post[c];
	}
	if(invar){
	  post[nclass-1]=wo[nclass-1]*fro[aa[i]+(nclass-1)*nchar];
	  marg += post[nclass-1];
	}
      }

      *lnl += log(marg);
      for(c=0; c<nclass; c++){
	post[c] /= marg;
	w[c] += post[c];
      }

      if(!is_const[i])
	for(c=0; c<nclass_change; c++)
	  for(j=0; j<nchar; j++) 
	    fr[j+c*nchar] += dn[j+i*nchar]*post[c];
      if(is_const[i]){
	for(c=0; c<nclass_change; c++)
	  fr[aa[i]+c*nchar] += dntot[i]*post[c];
	if(invar) fr[aa[i]+(nclass-1)*nchar] += post[nclass-1];
      }
    }
    if(pr) printf("%i %f\n",iter,*lnl);
    
    if(C>0) /* penalty adjustment */
      for(c=0; c<nclass_change; c++)
	for(j=0; j<nchar; j++) fr[j+c*nchar] +=  C; 
	  /* if(w[c]>0.0) fr[j+c*nchar] +=  C; */
    /* if(w[c]>0) might prevent fr[j+c*nchar] from becoming equal-freq class
     * and causing w[c] to incr, slowing converg, but isnan fr is possible */

    for(sum=0, c=0; c<nclass; c++) sum += w[c];
    for(c=0; c<nclass; c++) w[c] /= sum;
    for(c=0; c<nclass_change; c++){
      for(sum=0.0, j=0; j<nchar; j++) sum += fr[j+c*nchar];
      for(j=0; j<nchar; j++) fr[j+c*nchar] /= sum;
    }
    if(invar){
      for(sum=0.0, j=0; j<nchar; j++) sum += fr[j+(nclass-1)*nchar];
      for(j=0; j<nchar; j++) fr[j+(nclass-1)*nchar] /= sum;
    }

    if(iter>1){
      for(*dp=0, c=0; c<nclass; c++) *dp += fabs(w[c]-wo[c]);
      for(c=0; c<nclass; c++)
	for(j=0; j<nclass; j++) *dp += fabs(fr[j+c*nchar]-fro[j+c*nchar]);
      if((*dp)<(*dp_tol) || ((*lnl)-(*lnlo))<(*dlnl_tol)) break;
    }
  }

  return(iter);
}
