#ifndef Heston_temp_h
#define Heston_temp_h

#include<iostream>
#include<cmath>
#include<vector>
using namespace std;
double normcdf(double x){
    return 0.5* erfc(-x/sqrt(2));
}
double HestonDe(double s0, double v0,double r,double kappa, double theta,double sigma,double rho, double T, double k){
    double vsr  = theta *T +(v0-theta)*(1-exp(-kappa*T))/kappa;
    double d1 = (log(s0/k)+r*T+0.5*vsr)/sqrt(vsr);
    double d2 = d1 - sqrt(vsr);
    double dNd1 = exp(-d1*d1/2)/sqrt(2*M_PI);
    //double dNd2 = exp(-d2*d2/2)/sqrt(2*M_PI);
    double simu = exp(-r*T)*k*normcdf(-d2)- s0*normcdf(-d1)-rho*sigma*s0*dNd1*d2*((v0-2*theta)
                                                                                  *(1-exp(-kappa*T))/kappa+T*(theta-(v0-theta)*exp(-kappa*T)))/kappa/(2*vsr)+pow(sigma,2)*s0*dNd1
    *(d1*d2-1)/8.0/pow(kappa,2)/vsr*(3.0/2)*(theta*T+(v0-theta)*(1-exp(-kappa*T))/kappa
                                             -2*theta*(1-exp(-kappa*T))/kappa-2*T*(v0-theta)
                                             *exp(-kappa*T)+theta*(1-exp(-2*kappa*T))/kappa/2
                                             +(v0-theta)*(exp(-kappa*T)-exp(-2*kappa*T))/kappa);
    return simu;
}

#endif /* Heston_temp_h */
