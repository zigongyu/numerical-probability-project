#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

struct varcvar{
    varcvar(double var, double cvar): var(var), cvar(cvar) {}
    friend std::ostream & operator<<(std::ostream & o, varcvar const & vcv){
        return o << "VaR: " << vcv.var << "\tCVaR: " << vcv.cvar;
    }
    protected:
        double var,cvar;
};

double H1(double xi, double x, double alpha){
    return 1-(x>=xi)/(1-alpha);
}

double H2(double xi, double x, double c, double alpha){
    return c-xi-(x>=xi)*(x-xi)/(1-alpha);
}

template <typename TDistrib, typename TGen>
varcvar algo_naive(TDistrib & X, TGen & gen, unsigned sample_size, double alpha){
    double var=0, cvar=0, xi=0, c=0, xi_next, c_next, gamma = 1;
    for (unsigned k=0; k < sample_size; k++){
        double x = X(gen);
        xi_next = xi - gamma * H1(xi,x,alpha);
        c_next = c - gamma * H2(xi,x,c,alpha);
        var -= 1./(k+1)*(var-xi_next);
        cvar -= 1./(k+1)*(cvar-c_next);
        c = c_next;
        xi = xi_next;
        gamma = gamma*(1+k)/(k+2); 
    }
    varcvar s(var,cvar);
    return s;
}