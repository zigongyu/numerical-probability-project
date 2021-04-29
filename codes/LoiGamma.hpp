
class LoiGamma {
    public :
        LoiGamma ( double k, double theta): k(k), theta(theta){}

        double density(double x);
        double fctRepar(double x);
        double density_grad(double x);
    private :
        int k;
        int theta;
};