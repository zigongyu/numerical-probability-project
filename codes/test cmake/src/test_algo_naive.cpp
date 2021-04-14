#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>
#include <random>
#include "algo.cpp"

using namespace std;
int main() {
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    //normal_distribution<double> G;
    exponential_distribution<double> G(1);
    
    unsigned N =1000000;
    double alpha = 0.5;
    auto vcv = algo_naive(G, gen, N, alpha);
    cout << vcv << "\n";
    return 0;
}