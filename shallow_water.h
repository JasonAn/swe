#define ksyn   0.4      //synchronization parameter k

#define tstart 00    //synchronization time t
#define tsyn   00    //synchronization time t

#define Dm  8           //Time Delay Dimension
#define tau 10          //Time Delay

#define eps 0.001       //delta x

// from f_calc.c

void Fcalc(double *fields_dot, const double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, long int ncycle);

void RK4(double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, long int ncycle);

// from init.c

void initialize_geometry(int xdim, int ydim, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order);

void initialize_fields(double *u, double *v, double *P, double *x_forcing, double *y_forcing, double *z_forcing, int xdim, int ydim, double dx, double dy, double EquilibriumDepth, double A, int **n, unsigned int *lat, unsigned int *lon);

//from io.c

void get_parameters(char *file_name, int *xdim, int *ydim, double *dx, double *dy, long int *n_iter, double *dt, int *write_out_interval, double *f0, double *beta, double *forcingTerm, double *dissipationTerm, double *RayleighFriction, double *EquilibriumDepth, double *A);

void read_field(double *fields_msr, int xdim, int ydim, int t, unsigned int *print_out_order);

void import_field(double *P, int xdim, int ydim, int t, unsigned int *print_out_order);

void print_parameters(int xdim, int ydim, double dx, double dy, long int n_iter, double dt, int write_out_interval, double f0, double beta, double forcingTerm, double dissipationTerm, double RayleighFriction, double EquilibriumDepth, double A);

void print_field(double *field, char *name, int time_step, int xdim, int ydim, unsigned int *print_out_order);

void print_jacobian(double *field, char *name, int time_step, int xdim, int ydim, int m, int n);

void print_svd(double *field, char *name, int time_step,  int xdim, int ydim, int m, int n, int lda );

//from jacobian.c

void jacobian(double *fields_dot, double *fields, const double *parameters, const double *forcing, int xdim, int ydim, double dx, double dy, double dt, int **n, unsigned int *lat, unsigned int *lon, unsigned int *print_out_order, long int ncycle, int write_out_interval);

void dgesvd(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );

void dgelss( int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* s, double* rcond, int* rank, double* work, int* lwork, int* iwork, int* info );



