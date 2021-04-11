#include "LoiGamma.hpp"
#include "asa147.hpp" 
#include <cmath>
#include <cassert>
//#define _USE_MATH_DEFINES

double LoiGamma::density(double x){
    return pow(x,k-1)*exp(-x/theta)/(tgamma(k)*pow(theta,k));
}


double LoiGamma::fctRepar(double x){
    int ifault = 0;
    return gammds ( x/theta,k, &ifault)/tgamma(k);
}


// double LoiGamma::VaR(double alpha){
//     diff = 10000;
//     xi = 0;
//     while (diff > 0.0001){
//         new_xi = 
//     }
// }