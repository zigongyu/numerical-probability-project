class LoiGauss {
    public :
        LoiGauss ( double m = 0, double sigma = 1): m(m), sigma(sigma){}

        double density(double x);
        double fctRepar(double x);
  //      double VaR(double alpha);
    private :
        int m;
        int sigma;
};