#ifndef Heston_h
#define Heston_h

class HT_option : public base_option{
protected:
    double kappa;
    double theta;
    void print_model_para() {
        //  cout << "            = " << << endl;
        cout << "   kappa    = " << kappa << endl;
        cout << "   theta    = " << theta << endl;
    }
    
    double v_increment(double v, double W, double dt, double dt_sqrt) {
        return kappa*( theta - v ) * dt + gamma * pow(v,0.5) * W * dt_sqrt;
    }
    
    double s_increment(double s, double Z, double dt, double dt_sqrt, double v) {
        return r * s * dt + pow(v,0.5) * s * Z * dt_sqrt;
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
