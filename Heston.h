#ifndef Heston_h
#define Heston_h

class HT_option : public base_option{
    double normcdf(double x){
        return 0.5* erfc(-x/sqrt(2));
    }
protected:
    double kappa;
    double theta;
    void print_model_para() {
        //  cout << "            = " << << endl;
        cout << "   kappa    = " << kappa << endl;
        cout << "   theta    = " << theta << endl;
    }
    
    void fd_pricing(double rho, double T, vector<double>& res, bool display_flag) {
        // prameter initlisation
        clock_t c;
        double duration;
        double discount = exp(-r * T); // discouting factor
        
        double k = K;
        double sigma = gamma;
        
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
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        
        res.push_back(simu); // price
        res.push_back(duration); //time
    }
    
    
    double NW(double T) {
        cout << "NOT DONE at Heston.h " << endl;
        return 0;
    }
    
    double NN(double T) {
        cout << "NOT DONE at Heston.h " << endl;
        return 0;
    }
public:
    HT_option(const string& type, double s0, double K, double r, double gamma, double v0, int N, double kappa, double theta):base_option(type, s0, K, r, gamma, v0, N),kappa(kappa), theta(theta)  {
        // testing
        model = "Heston";
    }
};

#endif /* Heston_h */
