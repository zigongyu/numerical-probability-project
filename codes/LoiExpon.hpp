class LoiExpon {
    public :
        LoiExpon ( double lambda = 0): lambda(lambda){}

        double density(double x);
        double fctRepar(double x);
        double VaR(double alpha);
        double CVaR(double alpha);
    private :
        int lambda;
};