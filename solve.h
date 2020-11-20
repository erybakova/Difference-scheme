double rho(double t, double x);
double rho_t(double t, double x);
double rho_x(double t, double x);
double u(double t, double x);
double u_t(double t, double x);
double u_x(double t, double x);
double u_xx(double t, double x);
double g(double t, double x);
double g_t(double t, double x);
double g_x(double t, double x);
double f_0(double t, double x);
double f(double t, double x, double C, double mu);
double norm(double *G_curr, int M);
void init_V(double ntau, double *V_curr, double *G_curr,
                     double *a, double *c, double *d, double *b,
                        double C, double mu, double tau, double h, int M);
void init_G(double ntau, double *G_curr, double *V_next, double *V_curr,
                   double *a, double *c, double *d, double *b,
                        double mu, double tau, double h, int M);
void solve_threediag_system(double *a, double *c, double *d, double *b, double *x, int n);
void update(double *curr, double *next, int n);
void scheme(P_diff *p_diff, P_sche *p_sche,
             double *V_curr, double *G_curr,
             double *V_next, double *G_next);
double residual_C(double *V_curr, double *G_curr, double h, double T, int M);
