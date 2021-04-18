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
// #include "CalculVarCvar.hpp"
#include "algo.hpp"

using namespace std;


double phi(double x){
    return x;
}

double gNormal(double x){
    return abs(x);
}

int main(){
    // LoiExpon Expon(1);
    // cout << "Exp density: " << Expon.density(1) << endl;
    // cout << "Exp cdf: " << Expon.fctRepar(1.69) << endl;
    // cout << Expon.VaR(0.95) << endl;
    // cout << Expon.CVaR(0.95) << endl;
    // cout << endl;



    // LoiGamma Gam(0,1);
    // cout << "Gamma density: " << Gam.density(2) << endl;
    // cout << "Gamma cdf: " << Gam.fctRepar(2.64) << endl;
    // cout << endl;

    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);
    double (*fun)(double x);
    fun = phi;

    double (*G)(double x);
    G = gNormal;

    double alpha = 0.95;
    unsigned N = 1e6, M = 1e4;
    cout << "alpha = " << alpha <<", N = " << N << endl; 

    double lambda = 1;
    double m = 4,sigma = 15;
    // std::exponential_distribution<double> Exp(lambda);
    std::normal_distribution<double> Gauss(m,sigma);
    // for (int i = 0; i < 100; i++){
    //     cout << G(gen) << endl;
    // }

    double rou = 0.5, b = 1, c = 1 ;

    // calcul<decltype(Exp), decltype(gen),decltype(Expon)> cal(Exp, gen, Expon, rou, b, c);
    // cout << "Distribution Exp"<< endl;

    LoiGauss Gau(m,sigma);
    cout <<  "Gauss density: " <<Gau.density(1) << endl;
    cout << "Gauss cdf: " <<Gau.fctRepar(1.64) << endl;
    cout << endl;

    calcul<decltype(Gauss), decltype(gen),decltype(Gau)> cal(Gauss, gen, Gau, rou, b, c);
    cout << "Distribution Normal"<< endl;

    double f = 0.04;
    cout << cal.algo_naive( N, alpha, fun, f) << "\n\n";
    cout << cal.algo_avance(M, N, alpha, fun, f, G) << "\n\n";
    return 0;

}



