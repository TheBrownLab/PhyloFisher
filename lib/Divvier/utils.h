/////////////////////////////////////////////////////////////////
// utils.h
//
// Default constants for use in AMAP.  The emission
// probabilities were computed using the program used to build
// the BLOSUM62 matrix from the BLOCKS 5.0 dataset.  Transition
// parameters were obtained via unsupervised EM training on the
// BALIBASE 2.0 benchmark alignment database.
/////////////////////////////////////////////////////////////////



int readSeq(char *inFile);
void error(char *fmt, ... );
int pep2num(char c);
char *removeGaps(char *seq,int len,int *nlen);

#define addLogProb(x,y) ( (y) > (x) ) ? ((y) += log(1 + exp((x)-(y)))) : ((y)=(x)+log(1 + exp((y)-(x))))

#define MAX_LINE_LEN 2000
#define MAX_FILENAME_LEN 80
#define N_BASES 4
#define N_PEPT 20
#define MAX_LABEL_LEN 20

EXTERN int DO_HMM_APPROX;

EXTERN int Nseq;
EXTERN char **zorro_raw_seq;
EXTERN char **zorro_align;
EXTERN char **zorro_sequence;
EXTERN char **zorro_names;
EXTERN int *lens;
EXTERN int alen;
EXTERN double TOT_DIST;
EXTERN double **dists;
EXTERN int uguide; // User provided guide tree
EXTERN char guidetree[200];
#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))
