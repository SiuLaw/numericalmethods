#ifndef Heston_Laplace_h
#define Heston_Laplace_h

#include<complex>

const double pi = 3.14159265358979323846;

complex<double> Hestonlaplacetransform(complex<double>& x){
    
    complex<double> res;
    complex<double> kappa (4.0,0);
    complex<double> theta  (0.04,0);
    complex<double> gamma( 0.2,0);
    complex<double> v0 (0.04,0);
    complex<double> T (0.5,0);
    
    complex<double> zeta = sqrt(kappa*kappa +2.0*gamma*gamma*x);
    
    complex<double> B = -2.0*x/(kappa + zeta / tanh(zeta*T/2.0));
    
    complex<double> A = exp(kappa*kappa *theta *T/(gamma*gamma)/cosh(zeta * T/2.0)+kappa/zeta*sinh(zeta*T/2.0))*(2.0*kappa*theta/(gamma*gamma));
    
    res = A*exp(v0*B);
    
    return res;
    
}

double Hestonlaplacetransform(double x){
    
    double res;
    double kappa = 4.0;
    double theta = 0.04;
    double gamma = 0.2;
    double v0 = 0.04;
    double T = 0.5;
    
    double zeta = sqrt(kappa*kappa +2*gamma*gamma*x);
    
    double B = -2*x/(kappa + zeta / tanh(zeta*T/2));
    
    double A = exp(kappa*kappa *theta *T/(gamma*gamma)/cosh(zeta * T/2)+kappa/zeta*sinh(zeta*T/2))*(2*kappa*theta/(gamma*gamma));
    
    
    res = A*exp(v0*B);
    
    return res;
    
    
    
}

double HestonInverseLaplaceFT(double x){
    
    complex<double> res,theta,sigma;
    int M = 15;
    double r = 2*M/(5*x);
    double Sum = 0.0;
    complex<double> i(0.0,1.0);
    complex<double> S;
    
    for (int j = 1;j<M;j++){
        theta = j*pi/M;
        S = r*theta*(1.0/tan(theta)+i);
        sigma = theta+(theta/tan(theta)-1.0)/tan(theta);
        Sum = Sum + real(exp(x*S)* (1.0+i*sigma) * Hestonlaplacetransform(S));
    }
    
    res = r/M* (0.5*exp(r*x)*Hestonlaplacetransform(r) + Sum );
    
    return real(res);
}


double f(double w,double m,double x,double K,double r){
    
    double S0 = 100.0;
    
    return (K-S0*exp(w))*exp(-0.5*w-0.125*x+0.5*r*0.5+r*0.5*w/x - 0.5*r*r*0.25/x)*2*(2*m-w)/sqrt(2*pi)/(x*1.5)*exp(-(2*m-w)*(2*m-w)/2/x);
    
}


double ItgfuncfsrBintegral(double x,double K,double r,double B){
    
    double res = 0;
    double S0 = 100.0;
    
    double L_w = -10.0;
    double U_w = log(K/S0);
    double L_m = log(B/S0);
    double U_m = 10.0;
    
    int n = 1000;
    
    double dw = (U_w-L_w)/n;
    double dm = (U_m-L_m)/n;
    
    for (int i = 0; i < n; i++){
        double recSum = 0;
        for (int j = 0; j < n; j++){
            recSum += f(L_w + j*dw,L_m+i*dm,x,K,r)*dw;
        }
        res += recSum*dm;
    }
    
    return res;
}


double trapeziod(vector<double> &v){
    
    double res = 0.0;
    int s = v.size();
    
    for (int i = 0; i < s; i++){
        res += v[i];
    }
    
    res *= 2;
    res = res - v[0] - v[s-1];
    res /= 2;
    
    return res;
    
}


double HestonbarrierLaplaceFTrB (double K,double r,double B){
    
    double res = 0.0;
    
    vector<double> Y(40,0);
    
    double j;
    
    for (int k=1;k<41;k++){
        j = 1.0*k/400;
        Y[k-1] = HestonInverseLaplaceFT(j) * ItgfuncfsrBintegral(j,K,r,B);
    }
    
    res = 0.0025*trapeziod(Y)*exp(-r*0.5);
    
    return res;
}



#endif /* Heston_Laplace_h */
