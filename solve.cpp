#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include <math.h>

#include "init.h"
#include "solve.h"

double rho(double t, double x)
{
    return exp(t) * (cos(M_PI * x / 10.) + 1.5);
}

double rho_t(double t, double x)
{
    return rho(t, x);
}

double rho_x(double t, double x)
{
    return - (M_PI / 10.) * exp(t) * sin(M_PI * x / 10.);
}

double u(double t, double x)
{
    return cos(2. * M_PI * t) * sin(M_PI * x * x / 100.);
}

double u_t(double t, double x)
{
    return - (2. * M_PI) * sin(2. * M_PI * t) * sin(M_PI * x * x / 100.);
}

double u_x(double t, double x)
{
    return (M_PI / 50.) * x * cos(2. * M_PI * t) * cos(M_PI * x * x / 100.);
}

double u_xx(double t, double x)
{
    double coef = (M_PI / 50.) * cos(2. * M_PI * t);
    return coef * cos(M_PI * x * x / 100.) - coef * (M_PI * x * x / 50.) * sin(M_PI * x * x / 100.);
}

double g(double t, double x)
{
    return log(rho(t, x));
}

double g_t(double t, double x)
{
    (void) t;
    (void) x;
    return 1.;
}

double g_x(double t, double x)
{
    return rho_x(t, x) / rho(t, x);
}

double f_0(double t, double x)
{
    return rho_t(t, x) + u(t, x) * rho_x(t, x) + rho(t, x) * u_x(t, x);
}

double f(double t, double x, double C, double mu)
{
    return u_t(t, x) + u(t, x) * u_x(t, x)
            + (C  / rho(t, x)) * rho_x(t,x) - (mu / rho(t, x)) * u_xx(t, x);
}

double norm(double *G_curr, int M)
{
    int m;
    double max = fabs(exp(- G_curr[0]));

    for (m = 1; m < M + 1; m++)
    {
        max = (fabs(exp(- G_curr[m])) > max ? fabs(exp(- G_curr[m])) : max);
    }

    return max;
}

void init_V(double ntau, double *V_curr, double *G_curr,
                     double *a, double *c, double *d, double *b,
                        double C, double mu, double tau, double h, int M)
{
    int m;
    double mh, mu_tilde;

    mu_tilde = mu * norm(G_curr, M);

    a[0] = 1. + 2. * mu_tilde * tau / (h * h);
    c[0] = tau / (6. * h) * (V_curr[1] + V_curr[2]) - mu_tilde * tau / (h * h);
    d[0] = 0.;

    b[0] = V_curr[1]
            - tau / (2. * h) * C * (G_curr[2] - G_curr[0])
                - tau / (h * h) * (mu_tilde - mu * exp(- G_curr[1])) * (V_curr[2] - 2. * V_curr[1] + V_curr[0])
                    + tau * f(ntau, h, C, mu);

    for (m = 2; m < M - 1; m++)
    {
        mh = m * h;

        a[m - 1] = 1. + 2. * mu_tilde * tau / (h * h);
        c[m - 1] = tau / (6. * h) * (V_curr[m] + V_curr[m + 1]) - mu_tilde * tau / (h * h);
        d[m - 1] = - tau / (6. * h) * (V_curr[m] + V_curr[m - 1]) - mu_tilde * tau / (h * h);

        b[m - 1] = V_curr[m]
                    - tau / (2. * h) * C * (G_curr[m + 1] - G_curr[m - 1])
                        - tau / (h * h) * (mu_tilde - mu * exp(- G_curr[m])) * (V_curr[m + 1] - 2. * V_curr[m] + V_curr[m - 1])
                            + tau * f(ntau, mh, C, mu);
    }

    a[M - 2] = 1. + 2. * mu_tilde * tau / (h * h);
    c[M - 2] = 0.;
    d[M - 2] = - tau / (6. * h) * (V_curr[M - 1] + V_curr[M - 2]) - mu_tilde * tau / (h * h);

    b[M - 2] = V_curr[M - 1]
            - tau / (2. * h) * C * (G_curr[M] - G_curr[M - 2])
                - tau / (h * h) * (mu_tilde - mu * exp(- G_curr[M - 1])) * (V_curr[M] - 2. * V_curr[M - 1] + V_curr[M - 2])
                    + tau * f(ntau, (M - 1) * h, C, mu);
}

void init_G(double ntau, double *G_curr, double *V_next, double *V_curr,
                   double *a, double *c, double *d, double *b,
                        double mu, double tau, double h, int M)
{
    int m;
    double mh;

    (void) mu;

    a[0] = 1.;
    c[0] = tau / (2. * h) * V_next[1];
    d[0] = 0.;

    b[0] = G_curr[0]
            + tau / (2. * h) * (G_curr[0] - 2.) * V_next[1]
                + tau / (4. * h) * (4. * G_curr[2] * V_curr[2] - 5. * G_curr[1] * V_curr[1] - G_curr[3] * V_curr[3]
                                   + (2. - G_curr[0]) * (4. * V_curr[2] - 5. * V_curr[1] - V_curr[3]))
                    + tau * f_0(ntau, 0.) * exp(- G_curr[0]);

    for (m = 1; m < M; m++)
    {
        mh = m * h;

        a[m] = 1.;
        c[m] = tau / (4. * h) * (V_next[m] + V_next[m + 1]);
        d[m] =  - tau / (4. * h) * (V_next[m] + V_next[m - 1]);

        b[m] = G_curr[m]
                + tau / (4. * h) * (G_curr[m] - 2.) * (V_next[m + 1] - V_next[m - 1])
                    + tau * f_0(ntau, mh) * exp(- G_curr[m]);
    }

    a[M] = 1.;
    c[M] = 0.;
    d[M] = - tau / (2. * h) * V_next[M - 1];

    b[M] = G_curr[M]
            - tau / (2. * h) * (G_curr[M] - 2.) * V_next[M - 1]
                - tau / (4. * h) * (- 5. * G_curr[M - 1] * V_curr[M - 1] + 4. * G_curr[M - 2] * V_curr[M - 2] - G_curr[M - 3] * V_curr[M - 3]
                                   + (2. - G_curr[M]) * (- 5. * V_curr[M - 1] + 4. * V_curr[M - 2] - V_curr[M - 3]))
                    + tau * f_0(ntau, M * h) * exp(- G_curr[M]);
}

void solve_threediag_system(double *a, double *c, double *d, double *b, double *x, int n)
{
    int i;
    double m;

    for (i = 1; i < n; i++)
    {
        m = d[i] / a[i - 1];
        a[i] = a[i] - m * c[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    x[n - 1] = b[n - 1] / a[n - 1];

    for (i = n - 2; i >= 0; i--)
    {
        x[i] = (b[i] - c[i] * x[i+1]) / a[i];
    }
}

void update(double *curr, double *next, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        curr[i] = next[i];
    }
}

void scheme(P_diff *p_diff, P_sche *p_sche,
             double *V_curr, double *G_curr,
             double *V_next, double *G_next)
{
    double *a, *c, *d, *b;
    double h = p_sche -> h, tau = p_sche -> tau;
    double C = p_diff -> C, mu = p_diff -> mu;
    int M = p_sche -> M, N = p_sche -> N, Dim = p_sche -> Dim;
    double ntau;
    int n;

    /* диагональ */
    a = new double[Dim];
    /* верхняя диагональ */
    c = new double[Dim];
    /* нижняя диагональ */
    d = new double[Dim];

    /* правая часть */
    b = new double[Dim];

    /* цикл по слоям */
    for (n = 0; n < N; n++)
    {
        ntau = n * tau;

        init_V(ntau, V_curr, G_curr, a, c, d, b, C, mu, tau, h, M);
        solve_threediag_system(a, c, d, b, V_next + 1, M - 1);

        V_next[0] = 0.;
        V_next[M] = 0.;

        /*for (int m = 0; m <= M; m++)
        {
            double x = m * h;
            V_next[m] = u(ntau + tau, x);
        }*/

        init_G(ntau, G_curr, V_next, V_curr, a, c, d, b, mu, tau, h, M);
        solve_threediag_system(a, c, d, b, G_next, Dim);

        /*for (int m = 0; m <= M; m++)
        {
            double x = m * h;
            G_next[m] = g(ntau + tau, x);
        }*/

        /*printf("n = %d\n", n);
        for (int m = 0; m < M + 1; m++)
        {
            double curr_v = V_curr[m];
            double curr_g = G_curr[m];
            printf("curr_v %.3f curr_g %.3f\n", curr_v, curr_g);
        }
        printf("\n\n");*/

        update(V_curr, V_next, Dim);
        update(G_curr, G_next, Dim);
    }

    /*for (int m = 0; m < M + 1; m++)
    {
        double curr_v = V_curr[m];
        double curr_g = G_curr[m];
        printf("curr_v %.3f curr_g %.3f\n", curr_v, curr_g);
    }
    printf("\n");*/

    delete [] a;
    delete [] c;
    delete [] d;
    delete [] b;
}

double residual_C(double *V_curr, double *G_curr, double h, double T, int M)
{
    int m;
    double curr_g, curr_v;
    double max_v = 0, max_g = 0;

    for (m = 1; m < M; m++)
    {
        curr_v = fabs(V_curr[m] - u(T, m * h));
        curr_g = fabs(G_curr[m] - g(T, m * h));
        max_v = (curr_v > max_v ? curr_v : max_v);
        max_g = (curr_g > max_g ? curr_g : max_g);
        //printf("curr_v %.3f curr_g %.3f\n", curr_v, curr_g);
    }

    return (max_v > max_g ? max_v : max_g);
}
