#include <iostream>

// class varcvar{
//     public:
//         varcvar(double var, double cvar, double L1, double L2): 
//                     var(var), cvar(cvar), L_IC_var(L1), L_IC_cvar(L2) {}
//         friend std::ostream & operator<<(std::ostream & o, varcvar const & vcv){
//             o << "VaR: " << vcv.var << "\tCVaR: " << vcv.cvar << std::endl;
//             o << "la longueur de IC de VaR: " << vcv.L_IC_var <<"/(f^2)" << std::endl;
//             o << "la longueur de IC de CVaR: " << vcv.L_IC_cvar << std::endl;
//             return o;
//         }
//     protected:
//         double var,cvar,L_IC_var, L_IC_cvar;
// };

template <typename TDistrib, typename TGen>

class CalculVarCvar{
    public :
        CalculVarCvar(TDistrib X, TGen gen): X(X),gen(gen){}
        varcvar algo_naive(unsigned sample_size, double alpha, double f);

    private :
        // double var;
        // double cvar;
        TDistrib X;
        TGen gen;
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