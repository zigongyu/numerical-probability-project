#include <iostream>
#include <stdio.h>
#include <random>
#include <numeric>
#include <typeinfo>
//#include <vector>
//#include <string>
#include "LoiExpon.hpp"
#include "LoiGauss.hpp"
#include "LoiGamma.hpp"

#include "algo.hpp"

using namespace std;

double phi_linear(double x){
    return x+1;
}

double phi_linear_inverse(double x){
    return x-1;
}

double phi(double x){
    return x;
}

double phi_option(double x){
    double x0=100, sigma = 0.2, r = 0.05, T = 1, K = 110, P0 = 10.7;
    double s_T = x0 * exp((r-0.5*sigma*sigma) * T + sigma * sqrt(T) * x);
    double p_f = P0 * exp(r*T);
    return K > s_T ? K - s_T - p_f : -p_f;
}

double function_un(double x){
    return 1;
}


double gNormal(double x){
    return abs(x);
}

double gOption(double x){
    double x0=100, sigma = 0.2, r = 0.05, T = 1, K = 110, P0 = 10.7;
     double s_T = x0 * exp((r-0.5*sigma*sigma) * T + sigma * sqrt(T) * x);
     double p_f = P0 * exp(r*T);
    return K + s_T - p_f;
}


int main(){
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    double lambda = 1;
    double m = 0,sigma = 1;
    double k = 2, theta = 1;

    LoiExpon loiExpon(lambda);

    cout << "alpha = 0.95, phi(X) = X + 1 and X ~ Exp(" << lambda << ")" << endl;
    cout << "VaR of phi(X) is: " << 1-log(1-0.95)/lambda<<endl;
    cout << "CVaR of phi(X) is: " << 1+(1-log(1-0.95))/lambda<<endl;

    LoiGauss loiGauss(m,sigma);

    LoiGamma loiGamma(k,theta);

    std::exponential_distribution<double> Exp(lambda);
    std::normal_distribution<double> Gauss(m,sigma);
    std::gamma_distribution<double> Gamma(k,theta);

    double rou = 0.5, b = 1, c = 1 ;
    double (*G)(double x);
    G = gNormal;
    
    double (*fun)(double x);
    double (*fun_inv)(double x);
    fun = phi_linear;
    fun_inv = phi_linear_inverse;
    // fun = phi;
    // fun_inv = phi;

    double alpha = 0.95;
    unsigned N = 1e6, M = N/100;

    vector<unsigned> N_range{10000, 100000, 1000000};
    vector<double> alpha_range{0.5,0.75, 0.95,0.99,0.999};
    unsigned m_N = 0;
    
    calcul<decltype(Exp), decltype(gen),decltype(loiExpon)> cal_exp(Exp, gen, loiExpon, rou, b, c);
    cout<<"\n" << "Distribution Exp: "<< endl;
    for (unsigned n : N_range){
        cout <<  "N: " << n << endl;
        cout << cal_exp.algo_naive( n, alpha, fun, fun_inv) << "\n\n";

    }
    cout<<"\n\n" << "influence of alpha:"<< endl;
    for (double alpha_n : alpha_range){
        cout <<  "alpha: " << alpha_n << endl;
        cout << cal_exp.algo_naive( N, alpha_n, fun, fun_inv) << "\n";
    }
    cout << "-------------------------------------" <<"\n\n"<< endl;
    
    

    calcul<decltype(Gauss), decltype(gen),decltype(loiGauss)> cal_gauss(Gauss, gen, loiGauss, rou, b, c);
    cout<<"\n" << "Distribution Normale: "<< endl;
    cout<<"\n\n" << "influence of alpha:"<< endl;
    for (double alpha_n : alpha_range){
        cout <<  "alpha: " << alpha_n << endl;
        cout <<  "methode naive: " << endl;
        cout << cal_gauss.algo_naive( N, alpha_n, fun, fun_inv) << "\n";
        cout <<  "methode avance: " << endl;
        cout << cal_gauss.algo_avance(M, N, alpha_n, fun, fun_inv, G) << "\n";
    }

    cout<<"\n\n" << "-------------------------------------" <<"\n\n" <<endl;
  
    calcul<decltype(Gamma), decltype(gen),decltype(loiGamma)> cal_gamma(Gamma, gen, loiGamma, rou, b, c);
    cout<<"\n" << "Distribution Gamma: "<< endl;

    cout<<"\n\n" << "influence of alpha:"<< endl;
    for (double alpha_n : alpha_range){
        cout <<  "alpha: " << alpha_n << endl;
        cout <<  "methode naive: " << endl;
        cout << cal_gamma.algo_naive( N, alpha_n, fun, fun_inv) << "\n";
        cout <<  "methode avance: " << endl;
        cout << cal_gamma.algo_avance(M, N, alpha_n, fun, fun_inv, G) << "\n";
    }

    cout<<"\n\n" << "-------------------------------------" <<"\n\n"<< endl;
  

    cout << "Distribution Normale, option price: "<< endl;
    
    fun = phi_option;
    G = gOption;
    fun_inv = function_un;
    c = 0.5;
    calcul<decltype(Gauss), decltype(gen),decltype(loiGauss)> cal_gauss2(Gauss, gen, loiGauss, rou, b, c);
    for (double alpha_n : alpha_range){
        cout <<  "alpha: " << alpha_n << endl;
        cout <<  "methode naive: " << endl;
        cout << cal_gauss.algo_naive( N, alpha_n, fun, fun_inv) << "\n";
        cout <<  "methode avance: " << endl;
        cout << cal_gauss2.algo_avance(M, N, alpha_n, fun, fun_inv, G) << endl;
        cout << "the IC of VaR here is not accurate" << "\n\n";
    }
    return 0;

}



