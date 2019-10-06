
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>


#define EXTERN
#include "utils.h"
#undef EXTERN



char line[MAX_LINE_LEN];

double calc_dist(char *aln1,char *aln2,int len);
void protdist(char *inFile,double **mat);
char getNextBase(FILE **seqFile);
int getSeqLen(FILE *seqFile);
char  *readNextSeq(char *inFile,int *LEN,FILE *fp1,FILE *fp2,int *pflag,char **name);
char *removeGaps(char *seq,int len,int *nlen);

int Nseq;
char ** zorro_raw_seq;
char **zorro_align;
char **zorro_sequence;
char **zorro_names;
int *lens;
int alen;
double TOT_DIST;
int sample;
double **dists;
int uguide; // User provided guide tree
char guidetree[200];

/*
** A handy error function to prints  an error message and die
*/
void error(char *fmt, ... ) {
  va_list args;
  va_start(args, fmt);
  fprintf(stderr,"\n\nerror: ");
  vfprintf(stderr, fmt, args);
  fprintf(stderr,"  exiting...\n\n");
  va_end(args);
  exit(EXIT_FAILURE);
}

int readSeq(char *inFile){
  char *seq;
  int len,i;
  int tmp = 0, full_seq;
  FILE *fp1,*fp2;
  int pflag;
  if((fp1 = fopen(inFile,"r")) == NULL){
    error((char *)"readSeq: can't open %s for read, exiting\n",inFile);
  }
  // Calculate the number of sequences
  full_seq = tmp = 0;
  while(fgets(line,MAX_LINE_LEN,fp1) != NULL){
    if(line[0] == '>'){
      tmp++;
    }
  }
  full_seq = tmp;
  fclose(fp1);
  
  fp1 = fopen(inFile,"r");
  fp2 = fopen(inFile,"r");
  Nseq = 0;
  alen = -1;
  zorro_raw_seq = (char **)(malloc(tmp * sizeof(char *)));
  zorro_align = (char **)(malloc(tmp * sizeof(char *)));
  zorro_sequence = (char **)(malloc(tmp * sizeof(char *)));
  zorro_names = (char **)(malloc(tmp * sizeof(char *)));
  lens = (int *)(malloc(tmp * sizeof(int)));
  while(1){
    seq = readNextSeq(inFile,&len,fp1,fp2,&pflag,zorro_names+Nseq);
    if(len < 0){
      break;
    }
    if(alen == -1){
      alen = len;
      for(int k = 0 ; k < full_seq ; k++) {
    	  zorro_align[k] = (char*) ( malloc(alen * sizeof(char*)) );
    	  zorro_raw_seq[k] = (char*) ( malloc((alen + 1) * sizeof(char*)) );
      }
    }
    else if(len != alen){
      error("Wrong alignment\n");
    }
    //fprintf(stderr,">seq%d\n%s\n",Nseq,seq);
    // SW hacked so input is cleaner
    for(int k = 0; k < alen; k++) {
    	zorro_raw_seq[Nseq][k] = seq[k];
    	zorro_align[Nseq][k] = pep2num(seq[k]);
    }
    Nseq++;
  }
  if(tmp != Nseq){
    error("Non-MFA Sequence File\n");
  }
  fclose(fp1);
  fclose(fp2);
//  fprintf(stderr,"%d sequences, Length %d\n",Nseq,alen);

  for(i=0;i<Nseq;i++){
    zorro_sequence[i] = removeGaps(zorro_align[i],alen,&lens[i]);
    //fprintf(stderr,"%d\n",lens[i]);
  }

#ifdef NOWEIGHTING
  pairWeights = (double **)malloc(Nseq*sizeof(double *));

  for(i=0;i<Nseq;i++){
    pairWeights[i] = (double *)malloc(Nseq*sizeof(double)); 
  }
  protdist(inFile,pairWeights);
  
  dists = (double **)malloc(Nseq*sizeof(double *));

  for(i=0;i<Nseq;i++){
    dists[i] = (double *)malloc(Nseq*sizeof(double)); 
  }
  protdist(inFile,dists);
  
#endif
  

  return Nseq;
}


void protdist(char *inFile,double **mat){
  char command[200];
  char distfile[200];
  FILE *fp;
  int i,j,off1,off2;
  
  TOT_DIST = 0.0;

  sprintf(command,"/home/souravc/masking/calc_dist.pl %s %d",inFile,Nseq);
  system(command);
  fprintf(stderr,"Reading Distances\n");
  sprintf(distfile,"%s.trim.mat",inFile);
  if((fp = fopen(distfile,"r")) == NULL){
    error("Cannot open distance matrix %s\n",distfile);
  }
  for(i=0;i<Nseq;i++){
    if(fgets(line,MAX_LINE_LEN,fp) == NULL){
      error("%s : matrix in wrong format line %d\n",distfile,i+1);
    }
    off1 = 0;
    for(j=0;j<Nseq;j++){
      if(1 != sscanf(line+off1,"%lf%n",&mat[i][j],&off2)){
        error((char *) "bad format in distance matrix %s: %d %d %f %d\n%s\n",distfile,i,j,mat[i][j],off2,line+off1);
      }
      off1 += off2;
    }

  }
}

double calc_dist(char *aln1,char *aln2,int len){
  double f = 0.0;
  int i;
  for(i=0;i<len;i++){
    if(aln1[i] == aln2[i]){
      f += 1.0;
    }
  }
  f /= (double)(len); 
  return f;
}


char *removeGaps(char *seq,int len,int *nlen){
  char *tseq;
  int i,j;
  tseq = (char *)malloc(len*sizeof(char));
  for(i=0,j=0;i<len;i++){
    if(seq[i] != 20){
      tseq[j] = seq[i];
      j++;
    }
  }
  *nlen = j;
  return tseq;
}


char  *readNextSeq(char *inFile,int *LEN,FILE *fp1,FILE *fp2,int *pflag,char **name){

    char *seq;
    int i,len,counter;
    char c;
    
    len = getSeqLen(fp1);    
    *LEN = len;
    
    if(len == 0){
      fgets(line,MAX_LINE_LEN,fp2);
      fprintf(stderr,"Zero length sequence %s",line);
      //getSeqLen(fp2);
    }
    
    if(len <= 0){    
      return NULL;
    }
    
   
    
    seq = (char *)malloc(len*sizeof(char));
    
    
    c = getNextBase(&fp2);
    if(c != '>'){
      error((char *)"readNextSeq: %s not in FASTA format\nwrong character at start if seq %c\n",inFile,c);
    }  
    

    if(NULL == fgets(line,MAX_LINE_LEN,fp2))
      error((char *)"readNextSeq: problem reading first line of %s \n",inFile);
   
    
    
    for(i=0;i<strlen(line);i++){
      if(line[i]=='\0' || isspace(line[i])){
    	  if(isspace(line[i])) { line[i] = '\0'; }
    	  break;
      }
    }
    i--;
    
    *name = (char *)malloc((i+2)*sizeof(char));
    strncpy(*name,line+1,i+1);
 
    counter = 0;
    while(1){
      c = getNextBase(&fp2);
      if(c == EOF || c == '>' || c == '='){
	if(counter != len){
	  error((char *)"readNextSeq: %s not in FASTA format\nIncompatible length %d %d\n",inFile,counter,len);
	}
	if(c == '='){
	  *pflag = 1;
	}
	break;
      }
      
      if(c!= '>'){
    	  seq[counter] = c;
//	seq[counter] = pep2num(c);
      }
      else{
	error((char *)"readNextSeq: %s not in FASTA format\n",inFile);
      }
      counter++;
      if(counter > len){
	error((char *)"readNextSeq: %s not in FASTA format\n",inFile);
      }
    }

    return seq;
}




/*
** return the next character that is supposed to be a DNA base in
** the file.  This could still be an N or a * or something, checking
** isn't done here.  EOF returned when EOF or the end of the sequence
** is reached.
*/
char getNextBase(FILE **seqFile) {
  char c=-1;
  c = fgetc(*seqFile);
  for( ; ; c = fgetc(*seqFile) ) {    
    if(c == EOF) {
      break;
    } else if( isdigit(c) || isspace(c) || (c=='*') ) {
      continue;
    } else if(c == '='){
       fgetc(*seqFile);
       return c;
    }else if(c == '>') {
      /* ignore HEADER */   
      ungetc(c,*seqFile);
      return c;
    } else {
      break;
    }
  }
  return (c);
}


/*
"ARNDCQEGHILKMFPSTWYV"
*/
int pep2num(char c) {
  int r;
  
  switch(c) {
    case 'A':
    case 'a':
      return 0;
    break;
    case 'R':
    case 'r' :
      return 1;
    break;
    case 'N':
    case 'n' :
      return 2;
    break; 
    case 'D':
    case 'd' :
      return 3;
    break;
    case 'C':
    case 'c' :
      return 4;
    break;
    case 'Q':
    case 'q' :
      return 5;
    break;
    case 'E':
    case 'e' :
      return 6;
    break;
    case 'G':
    case 'g' :
      return 7;
    break;
    case 'H':
    case 'h' :
      return 8;
    break;
    case 'I':
    case 'i' :
      return 9;
    break;
    case 'L':
    case 'l':
      return 10;
    break;
    case 'K':
    case 'k':
      return 11;
    break;
    case 'M':
    case 'm':
      return 12;
    break;
    case 'F':
    case 'f':
      return 13;
    break;
    case 'P':
    case 'p' :
      return 14;
    break;
    case 'S':
    case 's' :
      return 15;
    break;
    case 'T':
    case 't' :
      return 16;
    break;
    case 'W':
    case 'w' :
      return 17;
    break;
    case 'Y':
    case 'y' :
      return 18;
    break;
    case 'V':
    case 'v' :
      return 19;
    break;
    case 'X':
    case 'x':
    case 'U':
    case 'u':
      r = (int)(20*(((double)rand())/((double)RAND_MAX)));
      return r;
    break;  
    case 'Z':
    case 'z':
      r = (int)(2*(((double)rand())/((double)RAND_MAX)));
      return (r+5);
    break;  
    case 'B':
    case 'b':
      r = (int)(2*(((double)rand())/((double)RAND_MAX)));
      return (r+2);
      break;
    case 'J':
    case 'j':
      r = (int)(2*(((double)rand())/((double)RAND_MAX)));
      return (r+9);
      break;
    case '-':
    case '.':
      return 20;
    break;
    default:
      error("Illegal character %c\n",c);
      return -1;
    break;
  }
}




/*
** Returns the length of the sequence in inFile.
*/
int getSeqLen(FILE *seqFile) {
  int seqLen;
  char base;

  
  
  if(NULL == fgets(line,MAX_LINE_LEN,seqFile)){
    return -1;
  }
    

  seqLen = 0;
  base = getNextBase(&seqFile);
  while (base != EOF) {
    if(base == '>'){  
      // Ignore header
      //fgets(line,MAX_LINE_LEN,seqFile);
      //base = getNextBase(&seqFile);
      break;
    }
    if(base == '='){
      break;
    }
    seqLen++;
    base = getNextBase(&seqFile);
  }

  return(seqLen);
}
