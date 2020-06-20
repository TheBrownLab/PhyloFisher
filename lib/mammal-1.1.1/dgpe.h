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
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
char get_param(int *k, char **arg, char *value);
void winterleavef(FILE *fp, int numsp, int nsite, char name[][11], char *seq);
void dorder(int n, double *x, int *o);
double *dgpe_rate_est(FILE *seqfile, FILE *treefile, FILE *seqfileo,
		      double quant, int *nsites, double *llikg, double *alpha);
void dist_est_paramf(FILE *fp, FILE **treefile, FILE **Qfile, FILE **ratefile, 
		     FILE **seqfile,  char *aaRmat, int *nchar, int *model, 
		     int *nrate, double *kappa, double *ub, int *de);
double lik_zero_rate_site(int nchar, int ntaxa, int *s, int incs, 
			  double *fr);
void numder(double x, double alpha, double *f);
void numderv(int nrates, double *wt, double *wtp, double *wtpp, double *br, 
	    double alpha);
double gamma_mle(int endsite, int nrates, double *rates, int *count,
		 double *likmat, double *alpha_max);
double lik_calc(int endsite, int nrates, int *count, double *likmat, double *wt);
void mean_rate_estf(int endsite, int nrates, double *rate, double *wt, 
		    double *likmat, double *ratee);
void makeprotfreqs(double *fr, double *W, double *Lam);
void chQ_freq(int nbin, double *W, double *Lam, double *fr, double *fro);
double int_elamrt(double t, double lam, double alpha, int deriv, double *de);
void pJCm(double t, int nchar, double *fr, int deriv, 
	    double *M, double *Mp, double *Mpp);
void pJCm_gen(double t, int nchar, double *fr, double *W, double *Lam, 
		int deriv, double *M, double *Mp, double *Mpp);
void pJCma(double t, int nchar, double *fr, int deriv, double alpha,
	   double *M, double *Mp, double *Mpp);
void pJCma_gen(double t, int nchar, double *fr, double *W, double *Lam, 
	       int deriv, double alpha, double *M, double *Mp, double *Mpp);
double muF81f(int nchar, double *fr);
void pF81m(double t, int nchar, double *fr, double mu, int deriv, 
	   double *M, double *Mp, double *Mpp);
void pF81m_gen(double t, int nchar, double *fr, double *mu, double *Lam, 
		 int deriv, double *M, double *Mp, double *Mpp);
void pF81ma(double t, int nchar, double *fr, double mu, int deriv, double alpha,
	   double *M, double *Mp, double *Mpp);
void pF81ma_gen(double t, int nchar, double *fr, double *mu, double *Lam, 
		int deriv, double alpha, double *M, double *Mp, double *Mpp);
int is_transitionf(int i, int j);
double muF84f(double K, double *fr);
void pF84m(double t, double mu, double K, double *fr, int deriv, 
	   double *M, double *Mp, double *Mpp);
void pF84m_gen(double t, int nchar, double *fr, double *mu, double *K, 
	       int deriv, double *M, double *Mp, double *Mpp);
void pF84ma(double t, double mu, double K, double *fr, int deriv, double alpha,
	    double *M, double *Mp, double *Mpp);
void pF84ma_gen(double t, int nchar, double *fr, double *mu, double *K, 
		int deriv, double alpha, double *M, double *Mp, double *Mpp);
double tt2K(double R, double *fr);
double K2tt(double K, double *fr);
double muHKYf(double kappa, double *fr);
void pHKYm(double t, int nchar, double *fr, double mu, double kappa,
	   int deriv, double *M, double *Mp, double *Mpp);
void pHKYm_gen(double t, int nchar, double *fr, double *mu, double *kappa, 
	       int deriv, double *M, double *Mp, double *Mpp);
void pHKYma(double t, int nchar, double *fr, double mu, double kappa,
	    int deriv, double alpha, double *M, double *Mp, double *Mpp);
void pHKYma_gen(double t, int nchar, double *fr, double *mu, double *kappa, 
	       int deriv, double alpha, double *M, double *Mp, double *Mpp);
double tt2kappa(double R, double *fr);
double kappa2tt(double kappa,  double *fr);
void rmnegpv(int m, double *pv, int incp);
void rmnegtrans(int m, double *P);
void Q2frWLam(int nbin, double *Q, double *fr, double *W, double *Lam);
void pGTRm(double t, int nchar, double *fr, double *W, double *Lam, 
	   int deriv, double *M, double *Mp, double *Mpp);
void pGTRma(double t, int nchar, double *fr, double *W, double *Lam, 
	    int deriv, double alpha, double *M, double *Mp, double *Mpp);
int subst_model_setparam(int nchar, int smodel, double *fr, double *Qk,
			 double **W, double **Lam);
int subst_model_setpm(int smodel,
		      void (*(*pm))(double t, int nchar, double *fr, 
				    double *W, double *Lam, int deriv, 
				    double *M, double *Mp, 
				    double *Mpp));
int subst_model_setpma(int smodel,
		       void (*(*pm))(double t, int nchar, double *fr, 
				     double *W, double *Lam, int deriv, 
				     double alpha,
				     double *M, double *Mp, double *Mpp));
void transm_rate(int smodel, int nchar, int ntaxa, double *utreec, 
		 int nrate, double *rate, double alpha,
		 double *fr, double *W, double *Lam, int deriv, 
		 double *M, double *Mp, double *Mpp);
int l2i(char c);
int l2ip(char c);
int rseqf(FILE *fp, char *seq);
char ignore_whitespace_seqfile(FILE *infile);
int rinterleave_block(FILE *infile, int numsp, int nsite, char name[][11], 
		       char *seq, int *initb);
void rinterleavef_nolim(FILE *infile, const int *numsp, const int *nsite, 
			char name[][11], char *seq);
void char_freq_count(int nchar, int numsp, int endsite, int *seq, int *count, 
		     double *freq);
int patt_freq_locs(int numsp, int nsite, int *seq, int *count, int *locs);
int read_tree(FILE *fp, char *tstring);
double get_dlist(char *tstr, int ss, int de, int *ssn, int *ls, int *le);
void underscore2blank(char *s);
void treecsnl(int ntaxa, char *tstr, char names[][11], double *utreec, 
	      int has_names);
char **treecs(int ntaxa, char *tstr, char names[][11], double *utreec, 
	      int has_names);
void read_treecsnl(FILE *treefile, int ntaxa, char names[][11], double *utreec, 
		   int has_names);
void backward_probb(int nchar, int r, double *M, double *f, double *bk,
		   double *b);
void backward_probsb(int nchar, int ntaxa, double *utreec, double *M, double *f, 
		    double *b);
void forward_probsb(int nchar, int ntaxa, int *s, double *utreec, double *M, 
		   double *f);
double lik_deriv_tedge(int nchar, int s, double *fr, double *bk, double *fl, 
		       double *Mp);
double lik_deriv_iedge(int nchar, double *fr, double *bk, double *fl, 
		     double *fc1, double *fc2, double *Mp);
void lik_deriv_calcb(int nchar, int ntaxa, int *s, double *fr, double *utreec, 
		     double *Mp, double *f, double *b, double *lp);
double lik_deriv_site_1rateb(int nchar, int ntaxa, int *s, double *fr, 
			     double *utreec, 
			     double *M, double *Mp, double *lp, int deriv);
double lik_deriv_siteb(int nchar, int ntaxa, int *s, double *fr, double *utreec, 
		      int nrate, double *rate, double *wt, double *M, 
		      double *Mp, double *lp, int deriv);
FILE *fopene(char *filename, char *rw, char *err);
char ctl_paraml(FILE *fp, char *opt, char *value, int *is_option);
int pmodel2smodel(int model, char *aaRmat, int nchar);
int bnewtr(double *a, double *b, double *x, double *f);
int eigenRealSym(double A[], int n, double Root[], double Offdiag[]);
void EigenSort(double d[], double U[], int n);
void HouseholderRealSym(double a[], int n, double d[], double e[]);
int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);
void cdfgam(int *which,double *p,double *q,double *x,double *shape,
	    double *scale,int *status,double *bound);
