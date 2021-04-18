#include "LoiGauss.hpp"
#include <cmath>
#include <cassert>
//#define _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846264338327950288

double LoiGauss::density(double x){
    return exp(-(x-m)*(x-m)/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}
double LoiGauss::density_grad(double x){
    return - exp(-(x-m)*(x-m)/(2*sigma*sigma)) * (x-m)/(pow(sigma,3)*sqrt(2*M_PI));
}

double LoiGauss::fctRepar(double x){
    return (1+erf((x-m)/(sigma*sqrt(2))))/2;
}


// double LoiGauss::VaR(double alpha){
//     diff = 10000;
//     xi = 0;
//     while (diff > 0.0001){
//         new_xi = 
//     }
// }