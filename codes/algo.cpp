#include "algo.hpp"
#include <cmath>

std::ostream & operator<<(std::ostream & o, varcvar const & vcv){
    o << "VaR: " << vcv.var << "\tCVaR: " << vcv.cvar << std::endl;
    o << "la longueur de IC de VaR: " << vcv.L_IC_var << std::endl;
    o << "la longueur de IC de CVaR: " << vcv.L_IC_cvar << std::endl;
    // o << "Real position at estimated VaR: " << vcv.estim << std::endl;
    return o;
}

double H1(double xi, double x, double alpha){
    return 1-(x>=xi)/(1-alpha);
}

double H2(double xi, double x, double c, double alpha){
    return c-xi-(x>=xi)*(x-xi)/(1-alpha);
}

double gamma_n(double n){
    return 1./(pow(n,3/4)+100);
}

