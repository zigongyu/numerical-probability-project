#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include <random>
#include "algo.hpp"

using namespace std;

struct loss_short_put{ // exemple 1 dans BFP09
    loss_short_put() = default;
    loss_short_put(double x0, double r, double s, double T, double K, double P0)
        :x0(x0), r(r), sigma(s), T(T), K(K), P0(P0) {};
    template<typename TGen>
    double operator()(TGen & gen){
        double s = x0 * exp((r-0.5*sigma*sigma) * T + sigma * sqrt(T) * G(gen));
        double p_f = P0 * exp(r*T);
        return K > s ? K - s - p_f : -p_f;
    }
    private:
        double x0,sigma,T,r,K,P0;
        std::normal_distribution<double> G;
};



int main() {
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    double alpha = 0.5;
    unsigned N = 1e6;
    cout << "alpha = 0.5, N = 1e6" << endl; 
    {
        std::exponential_distribution<double> Exp(2);
        auto vcv = algo_naive(Exp, gen, N, alpha);
        cout << "Distribution Exp(2) :" << endl;
        cout << vcv << "\n\n";

        std::normal_distribution<double> Gauss(3,1);
        auto vcv2 = algo_naive(Gauss, gen, N, alpha);
        cout << "Distribution Normal(3,1) :" << endl;
        cout << vcv2 << "\n\n";

        std::gamma_distribution<double> Gamma(2,4);
        auto vcv3 = algo_naive(Gamma, gen, N, alpha);
        cout << "Distribution Gamma(2,4) :" << endl;
        cout << vcv3 << "\n\n";
    }
    cout << "On refait l'exemple 1 dans BFG90(alpha = 0.95):\n";
    {
        double x0=100, s = 0.2, r = 0.05, T = 1, K = 110, P0 = 10.7;
        loss_short_put loss(x0,r,s,T,K,P0);
        alpha = 0.95;
        auto vcv4 = algo_naive(loss, gen, N, alpha);
        cout << vcv4 << endl;
    }
    
    
    return 0;
}