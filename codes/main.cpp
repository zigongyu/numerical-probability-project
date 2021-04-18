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


// struct loss_short_put{ // exemple 1 dans BFP09
//     loss_short_put() = default;
//     loss_short_put(double x0, double r, double s, double T, double K, double P0)
//         :x0(x0), r(r), sigma(s), T(T), K(K), P0(P0) {};
//     template<typename TGen>
//     double operator()(TGen & gen){
//         double s = x0 * exp((r-0.5*sigma*sigma) * T + sigma * sqrt(T) * G(gen));
//         double p_f = P0 * exp(r*T);
//         return K > s ? K - s - p_f : -p_f;
//     }
//     private:
//         double x0,sigma,T,r,K,P0;
//         std::normal_distribution<double> G;
// };

double addx(double x){
    return x + 1;
}

int main(){
    LoiExpon Expon(1);
    cout << "density: " << Expon.density(1) << endl;
    cout << "cdf: " << Expon.fctRepar(2.99) << endl;
    cout << Expon.VaR(0.95) << endl;
    cout << Expon.CVaR(0.95) << endl;
    cout << endl;

    LoiGauss Gau(0,1);
    cout <<  "density: " <<Gau.density(1) << endl;
    cout << "cdf: " <<Gau.fctRepar(1) << endl;
    //cout << Gauss.VaR(0.05) << endl;
    cout << endl;

    LoiGamma Gam(1,1);
    cout << "density: " << Gam.density(2) << endl;
    cout <<"cdf: " << Gam.fctRepar(2) << endl;
    //cout << Gauss.VaR(0.05) << endl;
    cout << endl;

    
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    double alpha = 0.95;
    unsigned N = 1e6;
    cout << "alpha = " << alpha <<", N = " << N << endl; 

    double lambda = 1;
    std::exponential_distribution<double> Exp(lambda);
    calcul<decltype(Exp), decltype(gen),decltype(Expon)> cal(Exp, gen, Expon);
    double f = 0.04;
    double (*fun)(double x);
    fun = addx;
    cout << cal.algo_naive( N, alpha, fun, f) << "\n\n";


    return 0;

}



