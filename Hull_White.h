#ifndef Hull_White_h
#define Hull_White_h

class HW_option : public base_option{
protected:
    double mu;
    void print_model_para() {
        //  cout << "            = " << << endl;
        cout << "   mu       = " << mu  << endl;
    }
    
    double v_increment(double v, double W, double dt, double dt_sqrt) {
        return mu * v * dt + gamma * v * W * dt_sqrt;
    }
    
    double s_increment(double s, double Z, double dt, double dt_sqrt, double v) {
        return r * s * dt + pow(v,0.5) * s * Z * dt_sqrt;
    }
    
    double NW(double T) {
        // PAGE 39
        double term_1 = -8.0 * gamma * pow( v0, 1.5 ) / mu;
        double term_2 = ( exp( 1.5 * mu * T + 9.0/8.0 * pow( gamma, 2.0) * T ) - 1.0 ) / ( 3.0 * ( 3.0 * pow(gamma,2.0) + 4.0 * mu )  );
        double term_3 = exp( mu*T ) * ( 1.0 - exp( 9.0/8.0 * pow(gamma,2.0) * T + 0.5 * mu * T ) ) / ( 9.0 * pow(gamma,2.0) + 4.0 * mu );
        
        return term_1 * ( term_2 + term_3 );
    }
    
    double NN(double T) {
        // PAGE 40
        double term_1 = -pow( v0, 2.0 );
        double term_2 = pow( gamma, 2.0 ) * ( 1.0 - exp(mu*T) ) * ( 1.0 - 3.0 * exp( mu * T ) ) / ( mu * ( 3.0 * pow(gamma,2.0) + mu ) * ( 3.0 * pow(gamma,2.0) + 2.0*mu ) );
        double term_3 = 3.0 * pow( gamma, 4.0 ) * pow( 1.0 - exp( mu * T ), 2.0 ) / ( pow(mu,2.0) * ( 3.0 * pow(gamma,2.0) + mu ) * ( 3.0 * pow(gamma,2.0) + 2.0* mu ) );
        double term_4 = 2.0 * exp( 2.0 * mu * T ) * ( 1.0 - exp( 3.0 * pow(gamma,2.0) * T )) / ( 3.0 * ( 3.0 * pow(gamma,2.0) + mu ) * ( 3.0 * pow(gamma,2.0) + 2.0 * mu ) );
        
        return term_1 * ( term_2 + term_3 + term_4 );
    }
public:
    HW_option(const string& type, double s0, double K, double r, double gamma, double v0, int N, double mu):base_option(type, s0, K, r, gamma, v0, N),mu(mu) {
        // testing
        model = "Hull White";
    }
};

#endif
