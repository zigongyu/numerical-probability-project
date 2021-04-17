#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

struct varcvar{
    varcvar(double var, double cvar, double L1, double L2): 
                var(var), cvar(cvar), L_IC_var(L1), L_IC_cvar(L2) {}
    friend std::ostream & operator<<(std::ostream & o, varcvar const & vcv){
        o << "VaR: " << vcv.var << "\tCVaR: " << vcv.cvar << std::endl;
        o << "la longueur de IC de VaR: " << vcv.L_IC_var <<"/(f^2)" << std::endl;
        o << "la longueur de IC de CVaR: " << vcv.L_IC_cvar << std::endl;
        return o;
    }
    protected:
        double var,cvar,L_IC_var, L_IC_cvar;
};

double H1(double xi, double x, double alpha){
    return 1-(x>=xi)/(1-alpha);
}

double H2(double xi, double x, double c, double alpha){
    return c-xi-(x>=xi)*(x-xi)/(1-alpha);
}

double gamma_n(double n){
    return 1./n;
}

template <typename TDistrib, typename TGen>
varcvar algo_naive(TDistrib & X, TGen & gen, unsigned sample_size, double alpha){
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
    L_IC_var = 2*1.96*sqrt(alpha*(1-alpha)/sample_size);
    L_IC_cvar = 2*1.96*sqrt((s_carre/sample_size-s*s/sample_size/sample_size)/(1-alpha)/(1-alpha)/sample_size);
    varcvar V(var,cvar,L_IC_var,  L_IC_cvar);
    return V;
}
