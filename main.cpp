#include <stdio.h>
#include <stdlib.h>

#include "init.h"
#include "solve.h"

int main(int argc, char *argv[])
{
    double *V_curr, *G_curr;
    double *V_next, *G_next;
    double mu, C;
    int M, N;
    double r_C, r_L2;
    int i;

    if (argc != 5 || ((M = atoi(argv[1])) <= 0) || ((N = atoi(argv[2])) <= 0)
            || ((C = atof(argv[3])) <= 0) || ((mu = atof(argv[4])) <= 0))
    {
        printf("Usage: %s M N C mu\n"
               "Note:  all parameters must be greater than 0\n", argv[0]);
        return 0;
    }

    P_diff p_diff;
    param_diff(&p_diff, C, mu);

    P_sche p_sche;
    param_sche(&p_sche, &p_diff, M, N);

    V_curr = new double[p_sche.Dim];
    G_curr = new double[p_sche.Dim];
    V_next = new double[p_sche.Dim];
    G_next = new double[p_sche.Dim];

    /* начальные условия */
    for (i = 0; i < p_sche.Dim; i++)
    {
        V_curr[i] = u(0., i * p_sche.h);
        G_curr[i] = g(0., i * p_sche.h);
    }

    /* приближенное решение на N-том слое */
    scheme(&p_diff, &p_sche, V_curr, G_curr, V_next, G_next);

    //! в конечном счете должно быть так
    /*
     * цикл для расчета всех тестов
     * 1) определяются параметры расчета дифференциальной схемы, param_dif
     * 2) определяются параметры схемы, param_sche
     * 3) вызывается подпрограмма расчета разностного решения, scheme
    */

    /* невязка */
    r_C = residual_C(V_curr, G_curr, p_sche.h, p_diff.T, M);
    printf("Residual, C-norm: %.4e\n", r_C);

    delete [] V_curr;
    delete [] V_next;
    delete [] G_curr;
    delete [] G_next;

    return 0;
}
