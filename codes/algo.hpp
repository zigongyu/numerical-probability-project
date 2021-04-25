#include <iostream>
#include <stdio.h>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

using namespace std;

class varcvar{
    public:
        varcvar(double var, double cvar, double L1, double L2, double estim): 
                    var(var), cvar(cvar), L_IC_var(L1), L_IC_cvar(L2), estim(estim){}
        friend std::ostream & operator<<(std::ostream & o, varcvar const & vcv);
        // double get_var();
    protected:
        double var,cvar,L_IC_var, L_IC_cvar, estim;
};



double H1(double xi, double x, double alpha);

double H2(double xi, double x, double c, double alpha);

double gamma_n(double n);

template <typename TDistrib, typename TGen, typename TLoi>
class calcul{
    public:
        calcul(TDistrib X, TGen gen, TLoi loi, double rou, double b, double c):X(X), gen(gen), 
        loi(loi), rou(rou), b(b), c(c){}



        varcvar algo_naive(unsigned sample_size, double alpha, double (*func)(double), double (*func_inv)(double)){
            double var=0, cvar=0, xi=0, cv=0, xi_next, cv_next, gamma;
            double L_IC_var, L_IC_cvar, s_carre=0, s=0;
            for (unsigned k=0; k < sample_size; k++){
                gamma = gamma_n(k+1);
                double x = func(X(gen));
                s_carre += (x-xi)*(x-xi)*(x>=xi);
                s += (x-xi)*(x>=xi);
                xi_next = xi - gamma * H1(xi,x,alpha);
                cv_next = cv - gamma * H2(xi,x,cv,alpha);
                var -= 1./(k+1)*(var-xi_next);
                cvar -= 1./(k+1)*(cvar-cv_next);
                cv = cv_next;
                xi = xi_next;
            }
            L_IC_var = 2*1.96*sqrt(alpha*(1-alpha)/sample_size)/loi.density(func_inv(var));
            L_IC_cvar = 2*1.96*sqrt((s_carre/sample_size-s*s/sample_size/sample_size)/(1-alpha)/(1-alpha)/sample_size);
            varcvar V(var,cvar,L_IC_var,  L_IC_cvar, loi.fctRepar(func_inv(var)) );
            return V;
        }


        
        double L1(double xi, double theta, double x, double alpha, double (*func)(double)){
            double cst = exp(- rou * pow(abs(theta),b) );
            double tmp = loi.density(x+theta)/loi.density(x);
            return cst * (1 - ( func(x+theta) >= xi) *tmp/(1-alpha));
        }

        double L2(double xi, double CV, double mu, double x, double alpha, double (*func)(double)){
            double tmp = loi.density(x+mu)/loi.density(x);
            double new_x = func(x+mu);
            return CV - xi -(new_x - xi) * (new_x >= xi) * tmp / (1-alpha);
        }

        double L3(double xi, double theta, double x, double alpha, double (*func)(double)){
            double new_x = func(x-theta);
            double cst = exp(- 2 * rou * pow(abs(theta),b) );
            double tmp = pow(loi.density(x-theta),2) * loi.density_grad(x-2*theta)/(loi.density(x) * pow(loi.density(x - 2*theta),2) );
            // cout << "loi.density(x), pow(loi.density(x - 2*theta),2)" << loi.density(x) << pow(loi.density(x - 2*theta),2)<< endl;
            // cout << "new_x, cst, temp: " << new_x << cst << tmp << endl;
            return cst * (new_x >= xi) * tmp;
        }

        double L4(double xi, double mu, double x, double alpha, double (*func)(double), double (*G)(double)){
            double new_x = func(x-mu);
            double cst = exp(- 2 * rou * pow(abs(mu),b) ) / (1 + pow(G(-mu),2*c) + pow(xi,2) ) ;
            double tmp = pow(loi.density(x-mu),2) * loi.density_grad(x-2*mu)/(loi.density(x) * pow(loi.density(x - 2*mu),2) );
            return cst * pow(new_x - xi, 2) * (new_x >= xi) * tmp;
        }

        varcvar algo_avance(unsigned M, unsigned N, double alpha, double (*func)(double), double (*func_inv)(double), double (*G)(double)){
            double var=0, cvar=0, gamma;
            double xi, CV=0, xi_old, CV_old, xi_bar = 0, CV_bar = 0;
            double theta, theta_old, mu, mu_old;
            double xi_hat=1, theta_hat=1, mu_hat=1;
            double xi_hat_old, theta_hat_old, mu_hat_old;
            double L_IC_var, L_IC_cvar;
            double alpha_n;
            double s1_carre = 0, s1 = 0, s2_carre = 0, s2 = 0;

            for (int n = 0; n < M; n++){
                if (n < M/3){
                    alpha_n = 0.5;
                }
                else if (n < 2*M/3){
                    alpha_n = 0.8;
                }
                else{
                    alpha_n = alpha;
                }
                gamma = gamma_n(n+1);
                double x = X(gen);

                xi_hat_old = xi_hat;
                theta_hat_old = theta_hat;
                mu_hat_old = mu_hat;

                xi_hat = xi_hat_old - gamma * H1(xi_hat_old, x, alpha_n);
                theta_hat = theta_hat_old - gamma * L3(xi_hat_old, theta_hat_old, x, alpha_n, func);
                mu_hat = mu_hat_old - gamma * L4(xi_hat_old, mu_hat_old, x, alpha_n, func, G);
            }
            
            theta = theta_hat;
            mu = mu_hat;
            xi = xi_hat;
            for (int n = 0; n < N; n++){
                gamma = gamma_n(n+1);
                double x = X(gen);

                xi_old = xi;
                CV_old = CV;
                theta_old = theta;
                mu_old = mu;
                xi = xi_old - gamma * L1(xi_old, theta_old, x, alpha, func);
                CV = CV_old - gamma * L2(xi_old, CV_old, mu_old, x, alpha, func);
                theta = theta_old - gamma * L3(xi_old, theta_old, x, alpha, func);
                mu = mu_old - gamma * L4(xi_old, mu_old, x, alpha, func, G);
                
                xi_bar = xi_bar -(xi_bar - xi) / (n+1);
                CV_bar = CV_bar -(CV_bar - CV) / (n+1);
                s1 = s1 + (func(x+theta) >= xi_bar) * loi.density(x+theta)/loi.density(x);
                s1_carre = s1_carre + pow((func(x+theta) >= xi_bar) * loi.density(x+theta)/loi.density(x), 2);
                s2 = s2 + (func(x+mu) - xi) * (func(x+mu) >= xi) * loi.density(x+mu)/loi.density(x);
                s2_carre = s2_carre + pow( (func(x+mu) - xi) * (func(x+mu) >= xi) * loi.density(x+mu)/loi.density(x), 2);
            }
            
            var = xi_bar;
            cvar = CV_bar;
            L_IC_var = 2*1.96*sqrt(s1_carre/N-pow(s1/N,2))/loi.density(func_inv(var))/sqrt(N) ;
            L_IC_cvar = 2*1.96*sqrt(s2_carre/N-pow(s2/N,2))/(1-alpha)/sqrt(N);
            varcvar V(var,cvar,L_IC_var,  L_IC_cvar, loi.fctRepar(func_inv(var)));
            return V;
        }

    
    private:
        TDistrib X;
        TGen gen;
        TLoi loi;
        double rou,b,c;

};


