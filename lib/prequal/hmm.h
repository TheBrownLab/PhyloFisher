

#define ADD
#define LONG
#define ENDTRANS
#define NCSTATE

#include <stdbool.h>

// Stuff originally in trim.h
#define addLogProb(x,y) ( (y) > (x) ) ? ((y) += log(1 + exp((x)-(y)))) : ((y)=(x)+log(1 + exp((y)-(x))))

#define MAX_LINE_LEN 2000
#define MAX_FILENAME_LEN 80
#define N_BASES 4
#define N_PEPT 20
#define MAX_LABEL_LEN 20

void error(char *fmt, ... );

EXTERN int Nseq;
EXTERN char **align;
EXTERN char **sequence;
EXTERN char **names;
EXTERN int *lens;
EXTERN int alen;
EXTERN double TOT_DIST;
EXTERN double **dists;
EXTERN int uguide; // User provided guide tree
EXTERN char guidetree[200];

int pep2num(char c);

// Original content
EXTERN double initDistribDefault[];
EXTERN int NUMSTATE;

#ifndef LONG 
EXTERN double gapOpenDefault;
EXTERN double gapExtendDefault;
#else
EXTERN double gapOpenDefault;
EXTERN double gapOpenDefaultL;
EXTERN double gapExtendDefault;
EXTERN double gapExtendDefaultL;
#endif



EXTERN double selfM;
EXTERN double MtoXY;
EXTERN double selfXY;
EXTERN double XYtoM;
#ifdef LONG
EXTERN double MtoXYL;
EXTERN double selfXYL;
EXTERN double XYLtoM;
#endif
#ifdef NCSTATE
EXTERN double MtoC;
EXTERN double NtoM;
EXTERN double selfN;
EXTERN double selfC;
#endif

EXTERN double emitSingleDefault[20];
EXTERN double emitPairsDefault[20][20];
void initHMM(int len);
void calc_posterior(int len);


// Simon's functions
bool SimonGetPosteriors(int len, double **retPosteriors, bool doApprox, int ApproxSize);
bool SimonGetSubsetPosteriors(int len, double **retPosteriors, int **runList, int runListLength, bool doApprox, int ApproxSize);
bool MakePosteriors(int X, int Y, bool doApprox);
