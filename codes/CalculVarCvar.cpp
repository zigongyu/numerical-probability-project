#include "CalculVarCvar.hpp"
template <class TDistrib, typename TGen>


varcvar CalculVarCvar::algo_naive(unsigned sample_size, double alpha, double f){
    double var=0, cvar=0, xi=0, c=0, xi_next, c_next, gamma;
    double L_IC_var, L_IC_cvar, s_carre=0, s=0;
    for (unsigned k=0; k < sample_size; k++){
        gamma = gamma_n(k+1);
        double x = X(gen);
        s_carre += (x-xi)*(x-xi)*(x>=xi);
        s += (x-xi)*(x>=xi);
        xi_next = xi - gamma * H1(xi,x,alpha);
        c_next = c - gamma * H2(xi,x,c,alpha);
        var -= 1./(k+1)*(var-xi_next);
        cvar -= 1./(k+1)*(cvar-c_next);
        c = c_next;
        xi = xi_next;
    }

    L_IC_var = 2*1.96*sqrt(alpha*(1-alpha)/sample_size)/f;
    L_IC_cvar = 2*1.96*sqrt((s_carre/sample_size-s*s/sample_size/sample_size)/(1-alpha)/(1-alpha)/sample_size)/f;
    varcvar V(var,cvar,L_IC_var,  L_IC_cvar);
    return V;
}
