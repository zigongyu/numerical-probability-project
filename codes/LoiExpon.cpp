#include "LoiExpon.hpp"
#include <cmath>
#include <cassert>

double LoiExpon::density(double x) {
    if (x<0){
        return 0;
    }
    else{
        return lambda * exp(-lambda*x);
    }
    
}

double LoiExpon::fctRepar(double x) {
    if (x<0){
        return 0;
    }
    else{
        return 1 - exp(-lambda*x);
    }
}

double LoiExpon::VaR(double alpha){
    assert(alpha <=1 && alpha >=0);
    return -log(1-alpha)/lambda;

}