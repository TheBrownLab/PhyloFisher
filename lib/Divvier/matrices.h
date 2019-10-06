
EXTERN double *eigmat;
EXTERN double **probmat;
EXTERN int JTT;
EXTERN int PMB;
EXTERN int PAM;
EXTERN int MATRICES;


void init_matrices(double *freqaa);
void makeprotfreqs(double *eigm,double **probm, double *freqaa);
void make_pmatrix(double matrix[20][20], double rat);

