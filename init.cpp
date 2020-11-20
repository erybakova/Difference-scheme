#include "init.h"

#define T_ 1
#define X_ 10

void param_diff(P_diff *p_diff, double C, double mu)
{
    p_diff -> T = T_;
    p_diff -> X = X_;
    p_diff -> C = C;
    p_diff -> mu = mu;
}

void param_sche(P_sche *p_sche, P_diff *p_diff, int M, int N)
{
    p_sche -> M = M;
    p_sche -> N = N;
    p_sche -> Dim = M + 1;
    p_sche -> h = (p_diff -> X) / M;
    p_sche -> tau = (p_diff -> T) / N;
}
