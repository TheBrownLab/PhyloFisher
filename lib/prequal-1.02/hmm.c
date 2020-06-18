#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>


#define EXTERN
#include "hmm.h"
#undef EXTERN


#ifndef NCSTATE
#ifndef LONG
double initDistribDefault[] = { 0.6080327034f, 0.1959836632f, 0.1959836632f };
#else
double initDistribDefault[] = { 0.6814756989f, 8.615339902e-05f, 8.615339902e-05f, 0.1591759622f, 0.1591759622 };
#endif
#else
#ifndef LONG
double initDistribDefault[] = { 0.6080327034f, 0.0, 0.0, 0.1959836632f, 0.1959836632f };
#else
double initDistribDefault[] = { 0.6814756989f, 0.0, 0.0, 0.0, 0.0, 0.15926215055f, 0.15926215055f };
#endif
#endif


#ifdef NCSTATE
#ifndef LONG
int NUMSTATE = 5;
#else
int NUMSTATE = 7;
#endif
#endif

#ifndef LONG 
double gapOpenDefault = 0.01993141696f ;
double gapExtendDefault =  0.7943345308f;
#else
double gapOpenDefault = 0.0119511066f;
double gapOpenDefaultL = 0.008008334786f;
double gapExtendDefault = 0.3965826333f;
double gapExtendDefaultL = 0.8988758326;
#endif

double selfM;
double MtoXY;
double selfXY;
double XYtoM;
#ifdef LONG
double MtoXYL;
double selfXYL;
double XYLtoM;
#endif
#ifdef NCSTATE
double MtoC = 0.01;
double NtoM = 0.05;
double selfN;
double selfC;
#endif

//char *alphabetDefault = (char *)("ARNDCQEGHILKMFPSTWYV");

double emitSingleDefault[20] = {
  0.07831005f, 0.05246024f, 0.04433257f, 0.05130349f, 0.02189704f,
  0.03585766f, 0.05615771f, 0.07783433f, 0.02601093f, 0.06511648f,
  0.09716489f, 0.05877077f, 0.02438117f, 0.04463228f, 0.03940142f,
  0.05849916f, 0.05115306f, 0.01203523f, 0.03124726f, 0.07343426f
};

double emitPairsDefault[20][20] = {
  {0.02373072f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00244502f, 0.01775118f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00210228f, 0.00207782f, 0.01281864f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00223549f, 0.00161657f, 0.00353540f, 0.01911178f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00145515f, 0.00044701f, 0.00042479f, 0.00036798f, 0.01013470f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00219102f, 0.00253532f, 0.00158223f, 0.00176784f, 0.00032102f, 0.00756604f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00332218f, 0.00268865f, 0.00224738f, 0.00496800f, 0.00037956f, 0.00345128f, 0.01676565f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00597898f, 0.00194865f, 0.00288882f, 0.00235249f, 0.00071206f, 0.00142432f, 0.00214860f, 0.04062876f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00114353f, 0.00132105f, 0.00141205f, 0.00097077f, 0.00026421f, 0.00113901f, 0.00131767f, 0.00103704f, 0.00867996f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00318853f, 0.00138145f, 0.00104273f, 0.00105355f, 0.00094040f, 0.00100883f, 0.00124207f, 0.00142520f, 0.00059716f, 0.01778263f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00449576f, 0.00246811f, 0.00160275f, 0.00161966f, 0.00138494f, 0.00180553f, 0.00222063f, 0.00212853f, 0.00111754f, 0.01071834f, 0.03583921f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00331693f, 0.00595650f, 0.00257310f, 0.00252518f, 0.00046951f, 0.00312308f, 0.00428420f, 0.00259311f, 0.00121376f, 0.00157852f, 0.00259626f, 0.01612228f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00148878f, 0.00076734f, 0.00063401f, 0.00047808f, 0.00037421f, 0.00075546f, 0.00076105f, 0.00066504f, 0.00042237f, 0.00224097f, 0.00461939f, 0.00096120f, 0.00409522f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00165004f, 0.00090768f, 0.00084658f, 0.00069041f, 0.00052274f, 0.00059248f, 0.00078814f, 0.00115204f, 0.00072545f, 0.00279948f, 0.00533369f, 0.00087222f, 0.00116111f, 0.01661038f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00230618f, 0.00106268f, 0.00100282f, 0.00125381f, 0.00034766f, 0.00090111f, 0.00151550f, 0.00155601f, 0.00049078f, 0.00103767f, 0.00157310f, 0.00154836f, 0.00046718f, 0.00060701f, 0.01846071f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00631752f, 0.00224540f, 0.00301397f, 0.00285226f, 0.00094867f, 0.00191155f, 0.00293898f, 0.00381962f, 0.00116422f, 0.00173565f, 0.00250962f, 0.00312633f, 0.00087787f, 0.00119036f, 0.00180037f, 0.01346609f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00389995f, 0.00186053f, 0.00220144f, 0.00180488f, 0.00073798f, 0.00154526f, 0.00216760f, 0.00214841f, 0.00077747f, 0.00248968f, 0.00302273f, 0.00250862f, 0.00093371f, 0.00107595f, 0.00147982f, 0.00487295f, 0.01299436f, 0.0f, 0.0f, 0.0f},
  {0.00039119f, 0.00029139f, 0.00021006f, 0.00016015f, 0.00010666f, 0.00020592f, 0.00023815f, 0.00038786f, 0.00019097f, 0.00039549f, 0.00076736f, 0.00028448f, 0.00016253f, 0.00085751f, 0.00015674f, 0.00026525f, 0.00024961f, 0.00563625f, 0.0f, 0.0f},
  {0.00131840f, 0.00099430f, 0.00074960f, 0.00066005f, 0.00036626f, 0.00070192f, 0.00092548f, 0.00089301f, 0.00131038f, 0.00127857f, 0.00219713f, 0.00100817f, 0.00054105f, 0.00368739f, 0.00047608f, 0.00102648f, 0.00094759f, 0.00069226f, 0.00999315f, 0.0f},
  {0.00533241f, 0.00169359f, 0.00136609f, 0.00127915f, 0.00119152f, 0.00132844f, 0.00178697f, 0.00194579f, 0.00071553f, 0.01117956f, 0.00914460f, 0.00210897f, 0.00197461f, 0.00256159f, 0.00135781f, 0.00241601f, 0.00343452f, 0.00038538f, 0.00148001f, 0.02075171f}
};



double **Xfmatrix;
double **Yfmatrix;
double **Mfmatrix;
double **Xbmatrix;
double **Ybmatrix;
double **Mbmatrix;
#ifdef LONG
double **XLfmatrix;
double **YLfmatrix;
double **XLbmatrix;
double **YLbmatrix;
#endif
#ifdef NCSTATE
double *XNfmatrix;
double *YNfmatrix;
double *XNbmatrix;
double *YNbmatrix;
double *XCfmatrix;
double *YCfmatrix;
double *XCbmatrix;
double *YCbmatrix;
#endif

#define my_min(a,b) ((a)<(b)?(a):(b))	// Duplicate required
#define my_max(a,b) ((a)>(b)?(a):(b))

double pxy;

double *posterior;
int *sample;
int *Xpos;
int *Ypos;
int *states;


void forward(char *seqX,char *seqY,int lenX,int lenY);
void backward(char *seqX,char *seqY,int lenX,int lenY);
void addPosterior(int X,int Y);


/////////////////////////////////////////// New code /////////////////////////////////////
double ** allPosteriors = NULL;
// Stuff for the pairHMM approximation
int *leftBounds,*rightBounds; // The upper and lower bounds on the DP matrix for each run;
double *checkSpaceX, *checkSpaceY;
void BuildBounds(int X, int Y); // Fills left and right bounds with appropriate values, assuming the longest sequence is going accross (index j) in the matrix
void approx_forward(char *seqX,char *seqY,int lenX,int lenY);
void approx_backward(char *seqX,char *seqY,int lenX,int lenY);
void RunMakePosteriors(int X, int Y, bool DoApprox, int ApproxSize);
int BOUND_SIZE = 50;

// Dumb function for progress spinner

void ProgressSpinner(int suffix, int total) {
	static int count = 0;
	static char progress_spinner [] = "/-\\|";
	printf("\r%c",progress_spinner[count++]);
	if(suffix >= 0) { printf(" %d",suffix); }
	if(total >= 0) { printf(" / %d", total); }
	fflush(stdout);
	if(count == 4) { count = 0; }
}



// Function that Simon Whelan has created to calculate the posterior probabilities
// ---
// Runs similar to calc_posterior, but
bool SimonGetPosteriors(int len, double **retPosteriors, bool doApprox, int approxSize) {
	int i, j;
	assert(allPosteriors == NULL);
	allPosteriors = retPosteriors;
	// Memory allocations
	posterior = (double *) (malloc((len) * sizeof(double)));
	sample = (int *) (malloc((len) * sizeof(int)));
	Xpos = (int *) (malloc((len) * sizeof(int)));
	Ypos = (int *) (malloc((len) * sizeof(int)));
	states = (int *) (malloc((len) * sizeof(int)));
	checkSpaceX = (double *) malloc(len * sizeof(double));
	checkSpaceY = (double *) malloc(len * sizeof(double));
	if (doApprox) {
		printf("\nDeclaring left/right %d",len+3);
		leftBounds = (int *) malloc((len + 3) * sizeof(int));
		rightBounds = (int *) malloc((len + 3) * sizeof(int));
	}
	for (i = 0; i < len; i++) {
		posterior[i] = 0.0;
		sample[i] = 0;
	}
	// Main loop over sequences here
//	printf("\nCollecting all posterior probabilities\n");
	for (i = 0; i < Nseq; i++) {
		for (j = i + 1; j < Nseq; j++) {
			ProgressSpinner(i + 1,Nseq);
			RunMakePosteriors(i,j,doApprox, approxSize);
		}
	}

	printf("\r ... Done Get Posteriors");

	allPosteriors = NULL;
	return true;
}

bool SimonGetSubsetPosteriors(int len, double **retPosteriors, int **runList, int runListLength, bool doApprox, int ApproxSize) {
	int i, j, k;
	assert(allPosteriors == NULL);
	allPosteriors = retPosteriors;

	// Memory allocations
	posterior = (double *) (malloc((len) * sizeof(double)));
	sample = (int *) (malloc((len) * sizeof(int)));
	Xpos = (int *) (malloc((len) * sizeof(int)));
	Ypos = (int *) (malloc((len) * sizeof(int)));
	states = (int *) (malloc((len) * sizeof(int)));
	checkSpaceX = (double *) malloc(len * sizeof(double));
	checkSpaceY = (double *) malloc(len * sizeof(double));
	if (doApprox) {
		leftBounds = (int *) malloc((len + 3) * sizeof(int));
		rightBounds = (int *) malloc((len + 3) * sizeof(int));
	}
	for (i = 0; i < len; i++) {
		posterior[i] = 0.0;
		sample[i] = 0;
	}
//	printf("\nCollecting subset of posterior probabilities\n");
	for (i = 0; i < Nseq; i++) {
		ProgressSpinner(i + 1, Nseq);
		for (k = 0; k < my_min(runListLength,Nseq); k++) {
			int j = runList[i][k];
			assert(j >= 0 && j < Nseq);
			if (j == i) {
				continue;
			}
			RunMakePosteriors(i, j, doApprox, ApproxSize);
		}
	}
	printf("\r ... Done Get Posteriors");

	allPosteriors = NULL;
	return true;
}

// Function that manages approximation during posterior estimation for X and Y
void RunMakePosteriors(int X, int Y, bool DoApprox, int ApproxSize) {
	int i;
	BOUND_SIZE = ApproxSize;
	if (lens[X] > lens[Y]) {
		i = Y;
		Y = X;
		X = i;
	}
	for(i = 0; i < 5; i++) {
		if(MakePosteriors(X,Y,DoApprox)) { break; }
//		printf("\nBound failed on %d - %d comparison", X, Y); fflush(stdout);
		BOUND_SIZE *= 2;		// Expand the bound on fail
	}
	if(i == 5) {
		printf("\nApproximate pairHMM bounding failed on %d - %d comparison. Enforcing complete computation.\n", X, Y);
		MakePosteriors(X,Y,false);
	}
	// Transfer to PP
	for(i = 0 ; i < lens[X] ; i ++) {
		if(checkSpaceX[i] > allPosteriors[X][i]) { allPosteriors[X][i] = checkSpaceX[i]; }
	}
	for(i = 0 ; i < lens[Y] ; i ++) {
		if(checkSpaceY[i] > allPosteriors[Y][i]) { allPosteriors[Y][i] = checkSpaceY[i]; }
	}
}
// Function that calculates the posteriors for sequence X and Y
bool MakePosteriors(int X, int Y, bool DoApprox) {
	int i, j, Xp, Yp;
	double f;
	double maxPP = 0;

	//	printf("\nDoing posterior %d %d\n", X, Y); fflush(stdout);

	// Calculate the forward and backward matrices; note that pxy (the forward likelihood) is calculated in forward(...)
	if (DoApprox) {
		BuildBounds(X, Y);
		approx_forward(sequence[X], sequence[Y], lens[X], lens[Y]);
		approx_backward(sequence[X], sequence[Y], lens[X], lens[Y]);
	} else {
		forward(sequence[X], sequence[Y], lens[X], lens[Y]);
		backward(sequence[X], sequence[Y], lens[X], lens[Y]);
	}

	if (DoApprox) {
		int maxIndex = 0;
		bool BoundEdge = false;
//		printf("\nComparing %d : %d\n", X, Y);
		// Do sequence X
		for (i = 1; i <= lens[X]; i++) { // Loop through the sequence
			maxPP = 0;
			BoundEdge = false;
			maxIndex = -1;

			double BoundCheck = 0.0;
			if(leftBounds[i - 1] !=0) { BoundCheck = exp(Mfmatrix[i][leftBounds[i - 1]] + Mbmatrix[i][leftBounds[i - 1]] - pxy); }
			if(rightBounds[i - 1] != lens[Y]) { BoundCheck = my_max(exp(Mfmatrix[i][rightBounds[i - 1]] + Mbmatrix[i][rightBounds[i - 1]] - pxy),BoundCheck); }
			if(BoundCheck > 0.0001) {
//				printf("\n[X=%d,Y=%d]: Bound problem: %f for bound %d (%d,%d)",X,Y,BoundCheck,BOUND_SIZE,leftBounds[i-1],rightBounds[i-1]);
				return false;
			}

			for (j = 1; j <= lens[Y]; j++) {
				// The approximate bit
				if (j < leftBounds[i - 1] || j > rightBounds[i - 1]) {
					f = 0;
				} else {
#ifdef ADD
					f = exp(Mfmatrix[i][j] + Mbmatrix[i][j] - pxy);
#else
					f = Mfmatrix[i][j]+Mbmatrix[i][j]-pxy;
#endif
				}
				if (f > maxPP) {
					if(j < leftBounds[i-1] + 5 || j - 5 > rightBounds[i - 1]) { BoundEdge = true; }
					else { BoundEdge = false; }
					maxIndex = j;
					maxPP = f;
				}
			}
			//			printf(" : %f", maxPP);
			if (BoundEdge) {
				if(maxPP > 0.01 && i > BOUND_SIZE && i < lens[X] - BOUND_SIZE) {
//					printf("\nWARNING: Bounds on posteriors might be too tight for sequence %d (cf %d). Position %d has max %f\n", X, Y, maxIndex , maxPP);
					return false;
				}
			}
			int m = -1, n;
			if(i == m) {
				printf("Site[%d]",i);
				for(n = 1; n < lens[Y]; n++) {
					if(n % 5 == 0) { printf("\n[[%d]]",n); }
					printf(" %f",exp(Mfmatrix[i][n] + Mbmatrix[i][n] - pxy)); }
				printf("\nBoundcheck: %f",exp(Mfmatrix[i][leftBounds[i - 1]] + Mbmatrix[i][leftBounds[i - 1]] - pxy));
				printf("\nDone with checkSpace %f and bound %d\n", maxPP,BOUND_SIZE); exit(-1);
			}
			checkSpaceX[i-1] = maxPP;
		}
		// Do sequence Y
		for (i = 1; i <= lens[Y]; i++) { // Loop through the sequence
			maxPP = 0;
			BoundEdge = false;
			maxIndex = -1;
			for (j = 1; j <= lens[X]; j++) {
				if (i < leftBounds[j-1] || i > rightBounds[j-1]) {
					f = 0;
				} else {
#ifdef ADD
					f = exp(Mfmatrix[j][i] + Mbmatrix[j][i] - pxy);
#else
					f = Mfmatrix[j][i]+Mbmatrix[j][i]-pxy;
#endif
				}
				if (f > maxPP) {
					if(i < leftBounds[j-1] + 5 || i - 5 > rightBounds[j - 1]) { BoundEdge = true; }
					else { BoundEdge = false; }
					maxIndex = j;
					maxPP = f;
				}
			}
			if (BoundEdge) {
				if(maxPP > 0.01 && i > BOUND_SIZE && i < lens[Y] - BOUND_SIZE) {
//					printf("\nWARNING: Bounds on posteriors might be too tight for sequence %d (cf %d). Position %d has max %f\n", Y, X, maxIndex , maxPP);
					return false;
				}
			}
			checkSpaceY[i-1] = maxPP;
		}
	} else {
		// Do sequence X
		for (i = 1; i <= lens[X]; i++) { // Loop through the sequence
			maxPP = 0;
			for (j = 1; j <= lens[Y]; j++) {
#ifdef ADD
				f = exp(Mfmatrix[i][j] + Mbmatrix[i][j] - pxy);
#else
				f = Mfmatrix[i][j]+Mbmatrix[i][j]-pxy;
#endif
				if (f > maxPP) { maxPP = f; }
			}
			checkSpaceX[i-1] = maxPP;
		}

		// Do sequence Y
		for (i = 1; i <= lens[Y]; i++) { // Loop through the sequence
			maxPP = 0;
			for (j = 1; j <= lens[X]; j++) {
#ifdef ADD
				f = exp(Mfmatrix[j][i] + Mbmatrix[j][i] - pxy);
#else
				f = Mfmatrix[j][i]+Mbmatrix[j][i]-pxy;
#endif
				if (f > maxPP) { maxPP = f; }

			}
			checkSpaceY[i-1] = maxPP;
		}
	}
	return true;
}

/////////////////////////////////////// End my new code /////////////////////////////////

void initHMM(int len) {
	int i, j;
	double f;

	// Forward Algorithm Probability matrix

	Xfmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	Yfmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	Mfmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	Xfmatrix++;
	Yfmatrix++;
	Mfmatrix++;

#ifdef LONG
	XLfmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	YLfmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	XLfmatrix++;
	YLfmatrix++;
#endif
#ifdef NCSTATE
	XNfmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	YNfmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	XCfmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	YCfmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	XNfmatrix++;
	YNfmatrix++;
	XCfmatrix++;
	YCfmatrix++;
#endif

	for (i = -1; i <= len + 1; i++) {
		Xfmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		Yfmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		Mfmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		Xfmatrix[i]++;
		Yfmatrix[i]++;
		Mfmatrix[i]++;
#ifdef LONG
		XLfmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		YLfmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		XLfmatrix[i]++;
		YLfmatrix[i]++;
#endif
	}

	// Backward Algorithm Probability Matrix

	Xbmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	Ybmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	Mbmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
#ifdef LONG
	XLbmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
	YLbmatrix = (double **) (malloc((len + 3) * sizeof(double *)));
#endif
#ifdef NCSTATE
	XNbmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	YNbmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	XCbmatrix = (double *) (malloc((len + 3) * sizeof(double)));
	YCbmatrix = (double *) (malloc((len + 3) * sizeof(double)));
#endif

	for (i = -1; i <= len + 1; i++) {
		Xbmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		Ybmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		Mbmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
#ifdef LONG
		XLbmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
		YLbmatrix[i] = (double *) (malloc((len + 3) * sizeof(double)));
#endif
	}
	// Normalize Parameters

	for (i = 0; i < N_PEPT; i++) {
		for (j = i + 1; j < N_PEPT; j++) {
			emitPairsDefault[i][j] = emitPairsDefault[j][i];
		}
	}

	f = 0.0;
	for (i = 0; i < N_PEPT; i++) {
		for (j = 0; j < N_PEPT; j++) {
			f = f + emitPairsDefault[i][j];
		}
	}

	for (i = 0; i < N_PEPT; i++) {
		for (j = 0; j <= i; j++) {
			emitPairsDefault[i][j] /= f;
		}
	}

	f = 0.0;

	for (i = 0; i < N_PEPT; i++) {
		f += emitSingleDefault[i];
	}

	for (i = 0; i < N_PEPT; i++) {
		emitSingleDefault[i] = emitSingleDefault[i] / f;
	}

	// Take LOGS
#ifndef LONG
	for(i=0;i<3;i++) {
		initDistribDefault[i] = log(initDistribDefault[i]);
	}
#else
	for (i = 0; i < 5; i++) {
		initDistribDefault[i] = log(initDistribDefault[i]);
	}
#endif

	for (i = 0; i < N_PEPT; i++) {
		emitSingleDefault[i] = log(emitSingleDefault[i]);
		for (j = 0; j < N_PEPT; j++) {
			emitPairsDefault[i][j] = log(emitPairsDefault[i][j]);
		}
	}

#ifndef NCSTATE
#ifndef LONG
	MtoXY = log(gapOpenDefault);
	selfM = log(1 - 2*gapOpenDefault);
	selfXY = log(gapExtendDefault);
	XYtoM = log(1 - gapExtendDefault);
#else
	MtoXY = log(gapOpenDefault);
	selfM = log(1 - 2*gapOpenDefault-2*gapOpenDefaultL);
	selfXY = log(gapExtendDefault);
	XYtoM = log(1 - gapExtendDefault);
	MtoXYL = log(gapOpenDefaultL);
	selfXYL = log(gapExtendDefaultL);
	XYLtoM = log(1 - gapExtendDefaultL);
#endif  
#else
#ifndef LONG
	MtoXY = log(gapOpenDefault);
	selfM = log(1 - 2*gapOpenDefault-2*MtoC);
	selfXY = log(gapExtendDefault);
	XYtoM = log(1 - gapExtendDefault);
#else
	MtoXY = log(gapOpenDefault);
	selfM = log(1 - 2 * gapOpenDefault - 2 * gapOpenDefaultL - 2 * MtoC);
	selfXY = log(gapExtendDefault);
	XYtoM = log(1 - gapExtendDefault);
	MtoXYL = log(gapOpenDefaultL);
	selfXYL = log(gapExtendDefaultL);
	XYLtoM = log(1 - gapExtendDefaultL);
#endif  
	selfN = log(1 - NtoM);
	selfC = log(1 - NtoM);
	NtoM = log(NtoM);
	MtoC = log(MtoC);
#endif

	//fprintf(stderr,"Initialized HMM\n");
	//fprintf(stderr,"selfM %f MtoXY %f MtoXYL %f\n",selfM,MtoXY,MtoXYL);
	//fprintf(stderr,"selfXY %f XYtoM %f\nselfXYL %f XYLtoM %f\n",selfXY,XYtoM,
	//  selfXYL,XYLtoM);
	//fprintf(stderr,"selfN %f selfC %f NtoM %f MtoC %f\n",selfN,selfC,NtoM,MtoC);

}

void forward(char *seqX, char *seqY, int lenX, int lenY) {
	int i, j;
	double tmp;

	// Initialize Forward Probabilities

	Xfmatrix[0][0] = -FLT_MAX;
	Yfmatrix[0][0] = -FLT_MAX;
#ifndef NCSTATE
	Mfmatrix[0][0] = 0;
#else
	Mfmatrix[0][0] = -FLT_MAX;
	XNfmatrix[0] = 0;
	XCfmatrix[0] = 0;
#endif  

#ifdef LONG
	XLfmatrix[0][0] = -FLT_MAX;
	YLfmatrix[0][0] = -FLT_MAX;
#endif  
	for (i = -1; i < alen + 1; i++) {
		Xfmatrix[-1][i] = -FLT_MAX;
		Yfmatrix[-1][i] = -FLT_MAX;
		Mfmatrix[-1][i] = -FLT_MAX;
		Xfmatrix[i][-1] = -FLT_MAX;
		Yfmatrix[i][-1] = -FLT_MAX;
		Mfmatrix[i][-1] = -FLT_MAX;
#ifdef LONG
		XLfmatrix[-1][i] = -FLT_MAX;
		YLfmatrix[-1][i] = -FLT_MAX;
		XLfmatrix[i][-1] = -FLT_MAX;
		YLfmatrix[i][-1] = -FLT_MAX;
#endif
	}

	Mfmatrix[1][1] = initDistribDefault[0]
			+ emitPairsDefault[(int) (seqX[0])][(int) (seqY[0])];
	Xfmatrix[1][0] = initDistribDefault[1] + emitSingleDefault[(int) (seqX[0])];
	Yfmatrix[0][1] = initDistribDefault[2] + emitSingleDefault[(int) (seqY[0])];
#ifdef LONG
	XLfmatrix[1][0] = initDistribDefault[3]
			+ emitSingleDefault[(int) (seqX[0])];
	YLfmatrix[0][1] = initDistribDefault[4]
			+ emitSingleDefault[(int) (seqY[0])];
#endif

#ifdef NCSTATE
	XNfmatrix[1] = initDistribDefault[NUMSTATE - 2]
			+ emitSingleDefault[(int) (seqX[0])];
	YNfmatrix[1] = initDistribDefault[NUMSTATE - 1]
			+ emitSingleDefault[(int) (seqY[0])];

	// Calculate XNfmatrix[i] and YNfmatrix[j]
	for (i = 2; i <= lenX; i++) {
		XNfmatrix[i] = XNfmatrix[i - 1] + selfN
				+ emitSingleDefault[(int) (seqX[i - 1])];
	}
	for (j = 2; j <= lenY; j++) {
		YNfmatrix[j] = YNfmatrix[j - 1] + selfN
				+ emitSingleDefault[(int) (seqY[j - 1])];
	}
#endif  

	for (i = 0; i <= lenX; i++) {
		for (j = 0; j <= lenY; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			// Calculate Mfmatrix[i][j]
			if (i != 0 && j != 0) {
				if (i != 1 || j != 1) {
					tmp = selfM + Mfmatrix[i - 1][j - 1];
					addLogProb(Yfmatrix[i-1][j-1]+XYtoM,tmp);
					addLogProb(Xfmatrix[i-1][j-1]+XYtoM,tmp);
#ifdef LONG
					addLogProb(YLfmatrix[i-1][j-1]+XYLtoM,tmp);
					addLogProb(XLfmatrix[i-1][j-1]+XYLtoM,tmp);
#endif	  
#ifdef NCSTATE
					if (i == 1 && j > 1) {
						addLogProb(YNfmatrix[j-1]+NtoM,tmp);
					}
					if (j == 1 && i > 1) {
						addLogProb(XNfmatrix[i-1]+NtoM,tmp);
					}
#endif
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					Mfmatrix[i][j]
							= tmp
									+ emitPairsDefault[(int) (seqX[i - 1])][(int) (seqY[j
											- 1])];
				}
			} else {
				Mfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate Xfmatrix[i][j]
			if (i != 0) {
				if (i != 1 || j != 0) {
					tmp = MtoXY + Mfmatrix[i - 1][j];
					addLogProb(selfXY + Xfmatrix[i-1][j],tmp);
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					Xfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqX[i - 1])];
				}
			} else {
				Xfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate Yfmatrix[i][j]
			if (j != 0) {
				if (i != 0 || j != 1) {
					tmp = MtoXY + Mfmatrix[i][j - 1];
					addLogProb(selfXY + Yfmatrix[i][j-1],tmp);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					Yfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqY[j - 1])];
				}
			} else {
				Yfmatrix[i][j] = -FLT_MAX;
			}
#ifdef LONG
			// Calculate XLfmatrix[i][j]
			if (i != 0) {
				if (i != 1 || j != 0) {
					tmp = MtoXYL + Mfmatrix[i - 1][j];
					addLogProb(selfXYL + XLfmatrix[i-1][j],tmp);
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					XLfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqX[i - 1])];
				}
			} else {
				XLfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate YLfmatrix[i][j]
			if (j != 0) {
				if (i != 0 || j != 1) {
					tmp = MtoXYL + Mfmatrix[i][j - 1];
					addLogProb(selfXYL + YLfmatrix[i][j-1],tmp);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					YLfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqY[j - 1])];
				}
			} else {
				YLfmatrix[i][j] = -FLT_MAX;
			}
#endif

			//fprintf(stderr,"%d %d : %f %f %f %f %f\n",i,j,Mfmatrix[i][j],Xfmatrix[i][j],Yfmatrix[i][j],XLfmatrix[i][j],YLfmatrix[i][j]);

		}
	}

#ifdef NCSTATE
	// Calculate XCfmatrix[i] and YCfmatrix[j]
	XCfmatrix[0] = -FLT_MAX;
	for (i = 1; i <= lenX; i++) {
		tmp = selfC + XCfmatrix[i - 1];
		addLogProb(MtoC+Mfmatrix[i-1][lenY],tmp);
		//fprintf(stderr,"Xf %d :%f\n",i,tmp);
		XCfmatrix[i] = tmp + emitSingleDefault[(int) (seqX[i - 1])];
		//fprintf(stderr,"Xf %d :%f\n",i,XCfmatrix[i]);
	}
	YCfmatrix[0] = -FLT_MAX;
	for (j = 1; j <= lenY; j++) {
		tmp = selfC + YCfmatrix[j - 1];
		addLogProb(MtoC+Mfmatrix[lenX][j-1],tmp);
		//fprintf(stderr,"Yf %d :%f\n",j,tmp);
		YCfmatrix[j] = tmp + emitSingleDefault[(int) (seqY[j - 1])];
		//fprintf(stderr,"Yf %d :%f\n",j,YCfmatrix[j]);
	}
#endif

#ifdef ENDTRANS
	pxy = Mfmatrix[lenX][lenY] + initDistribDefault[0];
	addLogProb(Yfmatrix[lenX][lenY]+initDistribDefault[1],pxy);
	addLogProb(Xfmatrix[lenX][lenY]+initDistribDefault[2],pxy);
#ifdef LONG
	addLogProb(YLfmatrix[lenX][lenY]+initDistribDefault[3],pxy);
	addLogProb(XLfmatrix[lenX][lenY]+initDistribDefault[4],pxy);
#endif
	//fprintf(stderr,"pxy = %f\n",pxy);
#ifdef NCSTATE
	addLogProb(YCfmatrix[lenY]+initDistribDefault[NUMSTATE-2],pxy);
	addLogProb(XCfmatrix[lenX]+initDistribDefault[NUMSTATE-1],pxy);
#endif
#else
	pxy = Mfmatrix[lenX][lenY];
	addLogProb(Yfmatrix[lenX][lenY],pxy);
	addLogProb(Xfmatrix[lenX][lenY],pxy);
#ifdef LONG
	addLogProb(YLfmatrix[lenX][lenY],pxy);
	addLogProb(XLfmatrix[lenX][lenY],pxy);
#endif
	//fprintf(stderr,"pxy = %f\n",pxy);
#ifdef NCSTATE
	addLogProb(YCfmatrix[lenY],pxy);
	addLogProb(XCfmatrix[lenX],pxy);
#endif
#endif

	//fprintf(stderr,"pxy = %f\n",pxy);

	//for(i=0;i<=lenY;i++){
	//fprintf(stderr,"%d : %f %f %f\n",i,Mfmatrix[2][i],Xfmatrix[2][i],Yfmatrix[2][i]);
	//}
}

void backward(char *seqX, char *seqY, int lenX, int lenY) {
	int i, j;
	double tmp;

#ifdef ENDTRANS  
	Mbmatrix[lenX][lenY] = initDistribDefault[0];
	Xbmatrix[lenX][lenY] = initDistribDefault[1];
	Ybmatrix[lenX][lenY] = initDistribDefault[2];
#ifdef LONG
	XLbmatrix[lenX][lenY] = initDistribDefault[3];
	YLbmatrix[lenX][lenY] = initDistribDefault[4];
#endif
#ifdef NCSTATE
	YCbmatrix[lenY] = initDistribDefault[NUMSTATE - 2];
	XCbmatrix[lenX] = initDistribDefault[NUMSTATE - 1];
#endif  
#else
	Mbmatrix[lenX][lenY] = 0;
	Xbmatrix[lenX][lenY] = 0;
	Ybmatrix[lenX][lenY] = 0;
#ifdef LONG
	XLbmatrix[lenX][lenY] = 0;
	YLbmatrix[lenX][lenY] = 0;
#endif
#ifdef NCSTATE
	YCbmatrix[lenY] = 0;
	XCbmatrix[lenX] = 0;
#endif 
#endif  

	for (i = lenX; i >= 0; i--) {
		Mbmatrix[i][lenY + 1] = -FLT_MAX;
		Xbmatrix[i][lenY + 1] = -FLT_MAX;
		Ybmatrix[i][lenY + 1] = -FLT_MAX;
#ifdef LONG
		XLbmatrix[i][lenY + 1] = -FLT_MAX;
		YLbmatrix[i][lenY + 1] = -FLT_MAX;
#endif
	}

	for (j = lenY; j >= 0; j--) {
		Mbmatrix[lenX + 1][j] = -FLT_MAX;
		Xbmatrix[lenX + 1][j] = -FLT_MAX;
		Ybmatrix[lenX + 1][j] = -FLT_MAX;
#ifdef LONG
		XLbmatrix[lenX + 1][j] = -FLT_MAX;
		YLbmatrix[lenX + 1][j] = -FLT_MAX;
#endif

	}

#ifdef NCSTATE
	// Calculate XCbmatrix[i] and YCbmatrix[j]
	for (i = lenX - 1; i >= 0; i--) {
		XCbmatrix[i] = selfC + XCbmatrix[i + 1]
				+ emitSingleDefault[(int) (seqX[i])];

		//fprintf(stderr,"Xb %d :%f\n",i,XCbmatrix[i]);
	}
	for (j = lenY - 1; j >= 0; j--) {
		YCbmatrix[j] = selfC + YCbmatrix[j + 1]
				+ emitSingleDefault[(int) (seqY[j])];
		//fprintf(stderr,"Xb %d :%f\n",j,YCbmatrix[j]);
	}
#endif

	for (i = lenX; i >= 0; i--) {
		for (j = lenY; j >= 0; j--) {
			if (i == lenX && j == lenY) {
				continue;
			}

			// Calculate Mbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = selfM
						+ emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])]
						+ Mbmatrix[i + 1][j + 1];
			}

			if (i != lenX) {
				addLogProb(Xbmatrix[i+1][j]+emitSingleDefault[(int)(seqX[i])]+MtoXY,tmp);
#ifdef LONG
				addLogProb(XLbmatrix[i+1][j]+emitSingleDefault[(int)(seqX[i])]+MtoXYL,tmp);
#endif
			}

			if (j != lenY) {
				addLogProb(Ybmatrix[i][j+1]+emitSingleDefault[(int)(seqY[j])]+MtoXY,tmp);
#ifdef LONG
				addLogProb(YLbmatrix[i][j+1]+emitSingleDefault[(int)(seqY[j])]+MtoXYL,tmp);
#endif	
			}

#ifdef NCSTATE
			if (i == lenX && j < lenY) {
				addLogProb(YCbmatrix[j+1]+emitSingleDefault[(int)(seqY[j])]+MtoC,tmp);
			}
			if (j == lenY && i < lenX) {
				addLogProb(XCbmatrix[i+1]+emitSingleDefault[(int)(seqX[i])]+MtoC,tmp);
			}
#endif

			Mbmatrix[i][j] = tmp;

			// Calculate Xbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYtoM
						+ emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (i != lenX) {
				addLogProb(selfXY + Xbmatrix[i+1][j] + emitSingleDefault[(int)seqX[i]],tmp);
			}

			Xbmatrix[i][j] = tmp;
			// Calculate Ybmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYtoM + emitPairsDefault[(int) seqX[i]][(int) seqY[j]]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (j != lenY) {
				addLogProb(selfXY + Ybmatrix[i][j+1] + emitSingleDefault[(int)seqY[j]],tmp);
			}

			Ybmatrix[i][j] = tmp;

#ifdef LONG
			// Calculate XLbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYLtoM
						+ emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (i != lenX) {
				addLogProb(selfXYL + XLbmatrix[i+1][j] + emitSingleDefault[(int)seqX[i]],tmp);
			}

			XLbmatrix[i][j] = tmp;
			// Calculate YLbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYLtoM + emitPairsDefault[(int) seqX[i]][(int) seqY[j]]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (j != lenY) {
				addLogProb(selfXYL + YLbmatrix[i][j+1] + emitSingleDefault[(int)seqY[j]],tmp);
			}

			YLbmatrix[i][j] = tmp;
#endif
			//fprintf(stderr,"%d %d: %f %f %f %f %f\n",i,j,Mbmatrix[i][j],Xbmatrix[i][j],Ybmatrix[i][j],XLbmatrix[i][j],YLbmatrix[i][j]);
		}
	}

#ifdef NCSTATE
	// Calculate XNbmatrix[i] and YNbmatrix[j]
	XNbmatrix[lenX] = -FLT_MAX;
	for (i = lenX - 1; i >= 1; i--) {
		tmp = selfN + XNbmatrix[i + 1] + emitSingleDefault[(int) (seqX[i])];
		//fprintf(stderr,"Xb %d :%f\n",i,tmp);
		addLogProb(NtoM+Mbmatrix[i+1][0]+emitSingleDefault[(int)(seqX[i])],tmp);
		XNbmatrix[i] = tmp;
		//fprintf(stderr,"Xb %d :%f\n",i,XNbmatrix[i]);
	}
	YNbmatrix[lenY] = -FLT_MAX;
	for (j = lenY - 1; j >= 1; j--) {
		tmp = selfN + YNbmatrix[j + 1] + emitSingleDefault[(int) (seqY[j])];
		//fprintf(stderr,"Yb %d :%f\n",j,tmp);
		addLogProb(NtoM+Mbmatrix[0][j+1] + emitSingleDefault[(int)(seqY[j])],tmp);
		YNbmatrix[j] = tmp;
		//fprintf(stderr,"Yb %d :%f\n",j,YNbmatrix[j]);
	}
#endif  

}

////////////////////////////////// Approximate PP routines //////////////////////////////////////

void BuildBounds(int X, int Y) {
	int i;
	assert(lens[X] <= lens[Y]);
	int diff = lens[Y] - lens[X];
//	printf("\nBuilding bounds [%d,%d] : %d", X,Y, BOUND_SIZE + diff);
	double ratio = (double) lens[Y] / (double) lens[X];
	for (i = 0; i < lens[X] + 1; i++) {
		int start_pos = (int) ((double) i * ratio);
		leftBounds[i] = my_max(0, my_min(start_pos - BOUND_SIZE - diff, lens[Y]) );
		rightBounds[i] = my_min(lens[Y], my_max(0,start_pos + BOUND_SIZE + diff));
		//		printf("\n[%d] : %d - %d", i,leftBounds[i],rightBounds[i]);
	}
}

void approx_forward(char *seqX, char *seqY, int lenX, int lenY) {
	int i, j;
	double tmp;

	int skip_count = 0, total_count = 0;

	assert(lenX <= lenY);

	// Initialize Forward Probabilities

	Xfmatrix[0][0] = -FLT_MAX;
	Yfmatrix[0][0] = -FLT_MAX;
#ifndef NCSTATE
	Mfmatrix[0][0] = 0;
#else
	Mfmatrix[0][0] = -FLT_MAX;
	XNfmatrix[0] = 0;
	XCfmatrix[0] = 0;
#endif

#ifdef LONG
	XLfmatrix[0][0] = -FLT_MAX;
	YLfmatrix[0][0] = -FLT_MAX;
#endif
	for (i = -1; i < alen + 1; i++) {
		Xfmatrix[-1][i] = -FLT_MAX;
		Yfmatrix[-1][i] = -FLT_MAX;
		Mfmatrix[-1][i] = -FLT_MAX;
		Xfmatrix[i][-1] = -FLT_MAX;
		Yfmatrix[i][-1] = -FLT_MAX;
		Mfmatrix[i][-1] = -FLT_MAX;
#ifdef LONG
		XLfmatrix[-1][i] = -FLT_MAX;
		YLfmatrix[-1][i] = -FLT_MAX;
		XLfmatrix[i][-1] = -FLT_MAX;
		YLfmatrix[i][-1] = -FLT_MAX;
#endif
	}

	Mfmatrix[1][1] = initDistribDefault[0]
			+ emitPairsDefault[(int) (seqX[0])][(int) (seqY[0])];
	Xfmatrix[1][0] = initDistribDefault[1] + emitSingleDefault[(int) (seqX[0])];
	Yfmatrix[0][1] = initDistribDefault[2] + emitSingleDefault[(int) (seqY[0])];
#ifdef LONG
	XLfmatrix[1][0] = initDistribDefault[3]
			+ emitSingleDefault[(int) (seqX[0])];
	YLfmatrix[0][1] = initDistribDefault[4]
			+ emitSingleDefault[(int) (seqY[0])];
#endif

#ifdef NCSTATE
	XNfmatrix[1] = initDistribDefault[NUMSTATE - 2]
			+ emitSingleDefault[(int) (seqX[0])];
	YNfmatrix[1] = initDistribDefault[NUMSTATE - 1]
			+ emitSingleDefault[(int) (seqY[0])];

	// Calculate XNfmatrix[i] and YNfmatrix[j]
	for (i = 2; i <= lenX; i++) {
		XNfmatrix[i] = XNfmatrix[i - 1] + selfN
				+ emitSingleDefault[(int) (seqX[i - 1])];
	}
	for (j = 2; j <= lenY; j++) {
		YNfmatrix[j] = YNfmatrix[j - 1] + selfN
				+ emitSingleDefault[(int) (seqY[j - 1])];
	}
#endif

	for (i = 0; i <= lenX; i++) {
		for (j = 0; j <= lenY; j++) {
			if (i == 0 && j == 0) {
				continue;
			}
			total_count++;
			// Calculate Mfmatrix[i][j]
			if (i != 0 && j != 0) {
				if (j < leftBounds[i - 1] || j > rightBounds[i - 1]) {
					Mfmatrix[i][j] = Xfmatrix[i][j] = Yfmatrix[i][j] = -FLT_MAX;
#ifdef LONG
					XLfmatrix[i][j] = YLfmatrix[i][j] = -FLT_MAX;
#endif
					skip_count++;
					continue;
				}
				if (i != 1 || j != 1) {
					tmp = selfM + Mfmatrix[i - 1][j - 1];
					addLogProb(Yfmatrix[i - 1][j - 1] + XYtoM, tmp);
					addLogProb(Xfmatrix[i - 1][j - 1] + XYtoM, tmp);
#ifdef LONG
					addLogProb(YLfmatrix[i-1][j-1]+XYLtoM,tmp);
					addLogProb(XLfmatrix[i-1][j-1]+XYLtoM,tmp);
#endif
#ifdef NCSTATE
					if (i == 1 && j > 1) {
						addLogProb(YNfmatrix[j-1]+NtoM,tmp);
					}
					if (j == 1 && i > 1) {
						addLogProb(XNfmatrix[i-1]+NtoM,tmp);
					}
#endif
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					Mfmatrix[i][j] = tmp + emitPairsDefault[(int) (seqX[i - 1])][(int) (seqY[j - 1])];
				}
			} else {
				Mfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate Xfmatrix[i][j]
			if (i != 0) {
				if (i != 1 || j != 0) {
					tmp = MtoXY + Mfmatrix[i - 1][j];
					addLogProb(selfXY + Xfmatrix[i-1][j],tmp);
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					Xfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqX[i - 1])];
				}
			} else {
				Xfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate Yfmatrix[i][j]
			if (j != 0) {
				if (i != 0 || j != 1) {
					tmp = MtoXY + Mfmatrix[i][j - 1];
					addLogProb(selfXY + Yfmatrix[i][j-1],tmp);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					Yfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqY[j - 1])];
				}
			} else {
				Yfmatrix[i][j] = -FLT_MAX;
			}
#ifdef LONG
			// Calculate XLfmatrix[i][j]
			if (i != 0) {
				if (i != 1 || j != 0) {
					tmp = MtoXYL + Mfmatrix[i - 1][j];
					addLogProb(selfXYL + XLfmatrix[i-1][j],tmp);
					assert((int) (seqX[i - 1]) >= 0 && (int) (seqX[i - 1]) < 20);
					XLfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqX[i - 1])];
				}
			} else {
				XLfmatrix[i][j] = -FLT_MAX;
			}
			// Calculate YLfmatrix[i][j]
			if (j != 0) {
				if (i != 0 || j != 1) {
					tmp = MtoXYL + Mfmatrix[i][j - 1];
					addLogProb(selfXYL + YLfmatrix[i][j-1],tmp);
					assert((int) (seqY[j - 1]) >= 0 && (int) (seqY[j - 1]) < 20);
					YLfmatrix[i][j] = tmp
							+ emitSingleDefault[(int) (seqY[j - 1])];
				}
			} else {
				YLfmatrix[i][j] = -FLT_MAX;
			}
#endif

			//fprintf(stderr,"%d %d : %f %f %f %f %f\n",i,j,Mfmatrix[i][j],Xfmatrix[i][j],Yfmatrix[i][j],XLfmatrix[i][j],YLfmatrix[i][j]);

		}
	}

#ifdef NCSTATE
	// Calculate XCfmatrix[i] and YCfmatrix[j]
	XCfmatrix[0] = -FLT_MAX;
	for (i = 1; i <= lenX; i++) {
		tmp = selfC + XCfmatrix[i - 1];
		addLogProb(MtoC+Mfmatrix[i-1][lenY],tmp);
		//fprintf(stderr,"Xf %d :%f\n",i,tmp);
		XCfmatrix[i] = tmp + emitSingleDefault[(int) (seqX[i - 1])];
		//fprintf(stderr,"Xf %d :%f\n",i,XCfmatrix[i]);
	}
	YCfmatrix[0] = -FLT_MAX;
	for (j = 1; j <= lenY; j++) {
		tmp = selfC + YCfmatrix[j - 1];
		addLogProb(MtoC+Mfmatrix[lenX][j-1],tmp);
		//fprintf(stderr,"Yf %d :%f\n",j,tmp);
		YCfmatrix[j] = tmp + emitSingleDefault[(int) (seqY[j - 1])];
		//fprintf(stderr,"Yf %d :%f\n",j,YCfmatrix[j]);
	}
#endif

#ifdef ENDTRANS
	pxy = Mfmatrix[lenX][lenY] + initDistribDefault[0];
	addLogProb(Yfmatrix[lenX][lenY]+initDistribDefault[1],pxy);
	addLogProb(Xfmatrix[lenX][lenY]+initDistribDefault[2],pxy);
#ifdef LONG
	addLogProb(YLfmatrix[lenX][lenY]+initDistribDefault[3],pxy);
	addLogProb(XLfmatrix[lenX][lenY]+initDistribDefault[4],pxy);
#endif
	//fprintf(stderr,"pxy = %f\n",pxy);
#ifdef NCSTATE
	addLogProb(YCfmatrix[lenY]+initDistribDefault[NUMSTATE-2],pxy);
	addLogProb(XCfmatrix[lenX]+initDistribDefault[NUMSTATE-1],pxy);
#endif
#else
	pxy = Mfmatrix[lenX][lenY];
	addLogProb(Yfmatrix[lenX][lenY],pxy);
	addLogProb(Xfmatrix[lenX][lenY],pxy);
#ifdef LONG
	addLogProb(YLfmatrix[lenX][lenY],pxy);
	addLogProb(XLfmatrix[lenX][lenY],pxy);
#endif
	//fprintf(stderr,"pxy = %f\n",pxy);
#ifdef NCSTATE
	addLogProb(YCfmatrix[lenY],pxy);
	addLogProb(XCfmatrix[lenX],pxy);
#endif
#endif

//	printf("pxy = %f\n", pxy);
//	printf("skip_count = %d / %d = %f", skip_count, total_count, (double) skip_count / (double) total_count);
	//for(i=0;i<=lenY;i++){
	//fprintf(stderr,"%d : %f %f %f\n",i,Mfmatrix[2][i],Xfmatrix[2][i],Yfmatrix[2][i]);
	//}
	// exit(-1);
}

void approx_backward(char *seqX, char *seqY, int lenX, int lenY) {
	int i, j;
	double tmp;

	assert(lenX <= lenY);

#ifdef ENDTRANS
	Mbmatrix[lenX][lenY] = initDistribDefault[0];
	Xbmatrix[lenX][lenY] = initDistribDefault[1];
	Ybmatrix[lenX][lenY] = initDistribDefault[2];
#ifdef LONG
	XLbmatrix[lenX][lenY] = initDistribDefault[3];
	YLbmatrix[lenX][lenY] = initDistribDefault[4];
#endif
#ifdef NCSTATE
	YCbmatrix[lenY] = initDistribDefault[NUMSTATE - 2];
	XCbmatrix[lenX] = initDistribDefault[NUMSTATE - 1];
#endif
#else
	Mbmatrix[lenX][lenY] = 0;
	Xbmatrix[lenX][lenY] = 0;
	Ybmatrix[lenX][lenY] = 0;
#ifdef LONG
	XLbmatrix[lenX][lenY] = 0;
	YLbmatrix[lenX][lenY] = 0;
#endif
#ifdef NCSTATE
	YCbmatrix[lenY] = 0;
	XCbmatrix[lenX] = 0;
#endif
#endif

	for (i = lenX; i >= 0; i--) {
		Mbmatrix[i][lenY + 1] = -FLT_MAX;
		Xbmatrix[i][lenY + 1] = -FLT_MAX;
		Ybmatrix[i][lenY + 1] = -FLT_MAX;
#ifdef LONG
		XLbmatrix[i][lenY + 1] = -FLT_MAX;
		YLbmatrix[i][lenY + 1] = -FLT_MAX;
#endif
	}

	for (j = lenY; j >= 0; j--) {
		Mbmatrix[lenX + 1][j] = -FLT_MAX;
		Xbmatrix[lenX + 1][j] = -FLT_MAX;
		Ybmatrix[lenX + 1][j] = -FLT_MAX;
#ifdef LONG
		XLbmatrix[lenX + 1][j] = -FLT_MAX;
		YLbmatrix[lenX + 1][j] = -FLT_MAX;
#endif

	}

#ifdef NCSTATE
	// Calculate XCbmatrix[i] and YCbmatrix[j]
	for (i = lenX - 1; i >= 0; i--) {
		XCbmatrix[i] = selfC + XCbmatrix[i + 1]
				+ emitSingleDefault[(int) (seqX[i])];

		//fprintf(stderr,"Xb %d :%f\n",i,XCbmatrix[i]);
	}
	for (j = lenY - 1; j >= 0; j--) {
		YCbmatrix[j] = selfC + YCbmatrix[j + 1]
				+ emitSingleDefault[(int) (seqY[j])];
		//fprintf(stderr,"Xb %d :%f\n",j,YCbmatrix[j]);
	}
#endif

	for (i = lenX; i >= 0; i--) {
		for (j = lenY; j >= 0; j--) {
			if (i == lenX && j == lenY) {
				continue;
			}
			// Calculate Mbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				// The approximate bit (Note no offset here)
				if (j < leftBounds[i] || j > rightBounds[i]) {
					Mbmatrix[i][j] = Xbmatrix[i][j] = Ybmatrix[i][j] = -FLT_MAX;
#ifdef LONG
					XLbmatrix[i][j] = YLbmatrix[i][j] = -FLT_MAX;
#endif
					continue;
				}

				tmp = selfM + emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])] + Mbmatrix[i + 1][j + 1];
			}

			if (i != lenX) {
				addLogProb(Xbmatrix[i+1][j]+emitSingleDefault[(int)(seqX[i])]+MtoXY,tmp);
#ifdef LONG
				addLogProb(XLbmatrix[i+1][j]+emitSingleDefault[(int)(seqX[i])]+MtoXYL,tmp);
#endif
			}

			if (j != lenY) {
				addLogProb(Ybmatrix[i][j+1]+emitSingleDefault[(int)(seqY[j])]+MtoXY,tmp);
#ifdef LONG
				addLogProb(YLbmatrix[i][j+1]+emitSingleDefault[(int)(seqY[j])]+MtoXYL,tmp);
#endif
			}

#ifdef NCSTATE
			if (i == lenX && j < lenY) {
				addLogProb(YCbmatrix[j+1]+emitSingleDefault[(int)(seqY[j])]+MtoC,tmp);
			}
			if (j == lenY && i < lenX) {
				addLogProb(XCbmatrix[i+1]+emitSingleDefault[(int)(seqX[i])]+MtoC,tmp);
			}
#endif

			Mbmatrix[i][j] = tmp;

			// Calculate Xbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYtoM
						+ emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (i != lenX) {
				addLogProb(selfXY + Xbmatrix[i+1][j] + emitSingleDefault[(int)seqX[i]],tmp);
			}

			Xbmatrix[i][j] = tmp;
			// Calculate Ybmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYtoM + emitPairsDefault[(int) seqX[i]][(int) seqY[j]]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (j != lenY) {
				addLogProb(selfXY + Ybmatrix[i][j+1] + emitSingleDefault[(int)seqY[j]],tmp);
			}

			Ybmatrix[i][j] = tmp;

#ifdef LONG
			// Calculate XLbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYLtoM
						+ emitPairsDefault[(int) (seqX[i])][(int) (seqY[j])]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (i != lenX) {
				addLogProb(selfXYL + XLbmatrix[i+1][j] + emitSingleDefault[(int)seqX[i]],tmp);
			}

			XLbmatrix[i][j] = tmp;
			// Calculate YLbmatrix[i][j]
			tmp = -FLT_MAX;
			if (i != lenX && j != lenY) {
				tmp = XYLtoM + emitPairsDefault[(int) seqX[i]][(int) seqY[j]]
						+ Mbmatrix[i + 1][j + 1];
			}
			if (j != lenY) {
				addLogProb(selfXYL + YLbmatrix[i][j+1] + emitSingleDefault[(int)seqY[j]],tmp);
			}

			YLbmatrix[i][j] = tmp;
#endif
			//fprintf(stderr,"%d %d: %f %f %f %f %f\n",i,j,Mbmatrix[i][j],Xbmatrix[i][j],Ybmatrix[i][j],XLbmatrix[i][j],YLbmatrix[i][j]);
		}
	}

#ifdef NCSTATE
	// Calculate XNbmatrix[i] and YNbmatrix[j]
	XNbmatrix[lenX] = -FLT_MAX;
	for (i = lenX - 1; i >= 1; i--) {
		tmp = selfN + XNbmatrix[i + 1] + emitSingleDefault[(int) (seqX[i])];
		//fprintf(stderr,"Xb %d :%f\n",i,tmp);
		addLogProb(NtoM+Mbmatrix[i+1][0]+emitSingleDefault[(int)(seqX[i])],tmp);
		XNbmatrix[i] = tmp;
		//fprintf(stderr,"Xb %d :%f\n",i,XNbmatrix[i]);
	}
	YNbmatrix[lenY] = -FLT_MAX;
	for (j = lenY - 1; j >= 1; j--) {
		tmp = selfN + YNbmatrix[j + 1] + emitSingleDefault[(int) (seqY[j])];
		//fprintf(stderr,"Yb %d :%f\n",j,tmp);
		addLogProb(NtoM+Mbmatrix[0][j+1] + emitSingleDefault[(int)(seqY[j])],tmp);
		YNbmatrix[j] = tmp;
		//fprintf(stderr,"Yb %d :%f\n",j,YNbmatrix[j]);
	}
#endif

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
      error("Illegal character <%c>\n",c);
      return -1;
    break;
  }
}

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
