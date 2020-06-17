#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#ifndef PI
   #define PI 3.14159265358979323846
#endif
enum {
  alanine = 0, arginine, asparagine, aspartic, cysteine, 
  glutamine, glutamic, glycine, histidine, isoleucine,
  leucine, lysine, methionine, phenylalanine, proline,
  serine, threonine, tryptophan, tyrosine, valine
};
int *read_seq(FILE *fp, int nchar, int rm_gap, int *ntaxa, int *nsite, 
	      char (**name)[11]);
char *read_seqc(FILE *fp, int *ntaxa, int *nsite, char (**name)[11]);
void letter(int numsp, int nsite, char *seq);
int l2i(char c);
void rmgap(int numsp, int *nsite, char *seqc);
int l2ip(char c);
void rinterleavef_nolim(FILE *infile, const int *numsp, const int *nsite, 
			char name[][11], char *seq);
int rinterleave_block(FILE *infile, int numsp, int nsite, char name[][11], 
		       char *seq, int *initb);
int rseqf(FILE *fp, char *seq);
char ignore_whitespace_seqfile(FILE *infile);
double *estimate_sigma(FILE *seqfile, int *ntaxa);
