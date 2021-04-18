
class LoiGamma {
    public :
        LoiGamma ( double k, double theta): k(k), theta(theta){}

        double density(double x);
        double fctRepar(double x);
        double density_grad(double x);
  //      double VaR(double alpha);
    private :
        int k;
        int theta;
};