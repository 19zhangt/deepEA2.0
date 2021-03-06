double *expm11 ( int n, double a[] );
double *expm2 ( int n, double a[] );
double *expm3 ( int n, double a[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
double *mexp_a ( int test, int n );
double *mexp_expa ( int test, int n );
int mexp_n ( int test );
void mexp_story ( int test );
int mexp_test_num ( void );
double r8_abs ( double x );
double r8_add ( double x, double y );
double r8_epsilon ( void );
double r8_huge ( void );
double r8_log_2 ( double x );
double r8_max ( double x, double y );
void r8mat_add ( int m, int n, double alpha, double a[], double beta, 
  double b[], double c[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double *r8mat_copy_new ( int m, int n, double a1[] );
double *r8mat_copy_new ( int m, int n, double a1[] );
double *r8mat_fss_new ( int n, double a[], int nb, double b[] );
double *r8mat_identity_new ( int n );
void r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] );
void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] );
double r8mat_norm_l1 ( int m, int n, double a[] );
double r8mat_norm_li ( int m, int n, double a[] );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title );
void r8mat_scale ( int m, int n, double s, double a[] );
int r8mat_significant ( int m, int n, double r[], double s[] );;
double *r8mat_zero_new ( int m, int n );
void timestamp ( void );
