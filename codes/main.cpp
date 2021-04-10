#include <iostream>
//#include <vector>
//#include <string>
#include "LoiExpon.hpp"
using namespace std;

int main(){
    LoiExpon Expon(1);
    cout << Expon.density(1) << endl;
    cout << Expon.fctRepar(1) << endl;
    cout << Expon.VaR(0.05) << endl;

}