#include <iostream>
#include <stdio.h>
//#include <vector>
//#include <string>
#include "LoiExpon.hpp"
#include "LoiGauss.hpp"
#include "LoiGamma.hpp"
using namespace std;

int main(){
    LoiExpon Expon(1);
    cout << Expon.density(1) << endl;
    cout << Expon.fctRepar(1) << endl;
    cout << Expon.VaR(0.05) << endl;
    cout << endl;


    LoiGauss Gau(0,1);
    cout << Gau.density(1) << endl;
    cout << Gau.fctRepar(1) << endl;
    //cout << Gauss.VaR(0.05) << endl;
    cout << endl;

    LoiGamma Gam(1,1);
    cout << Gam.density(2) << endl;
    cout << Gam.fctRepar(2) << endl;
    //cout << Gauss.VaR(0.05) << endl;
    cout << endl;



}