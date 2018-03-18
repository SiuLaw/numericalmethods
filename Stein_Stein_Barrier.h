#ifndef Stein_Stein_Barrier_h
#define Stein_Stein_Barrier_h


class SSB_option : public base_option{
protected:
    double kappa;
    double theta;
    double kappa_base;
    double theta_base;
    
    void print_model_para() {
        //  cout << "            = " << << endl;
        cout << "   kappa    = " << kappa_base << endl;
        cout << "   theta    = " << theta_base << endl;
        cout << "   kappa_Q  = " << kappa << endl;
        cout << "   theta_Q  = " << theta << endl;
    }
    
    double v_increment(double v, double W, double dt, double dt_sqrt) {
        return kappa*( theta - v ) * dt + gamma * W * dt_sqrt;
    }
    
    double s_increment(double s, double Z, double dt, double dt_sqrt, double v) {
        return r * s * dt + v * s * Z * dt_sqrt;
    }
    
    double NW(double T) {
        // PAGE 43
        double term_1 = gamma / pow( kappa, 2.0 );
        double term_2 = 0.5 * ( pow(theta,2.0) * ( 4.0*kappa*T - 9.0) + v0*(v0 + 4.0*theta) + pow(gamma,2.0) / kappa * ( kappa * T - 1.0 ) );
        double term_3 = ( 2.0*theta *kappa*T*(theta-v0) + 2.0*theta * ( 3.0*theta - 2.0 * v0 ) ) * exp( -kappa * T );
        double term_4 = 0.5 * ( pow(gamma,2.0)/kappa * ( kappa*T + 1.0 ) - (theta - v0 ) * ( 3.0*theta - v0 ) - 2.0 * kappa * T * pow( theta - v0, 2.0 )  ) * exp( -2.0*kappa*T );
        
        return term_1 * ( term_2 + term_3 + term_4 );
    }
    
    double NN(double T) {
        // PAGE 44
        double term_1 = pow(gamma,2.0) / pow( kappa, 3.0 );
        double term_2 = 0.5 * ( - 5.0 * pow(gamma,2.0)/( 4.0*kappa) + pow(v0,2.0) + 6.0 * theta * v0 - 19.0 * pow(theta, 2.0) + T * ( 8.0*pow(theta,2.0)*kappa + pow(gamma,2.0) ) );
        double term_3 = ( 2.0*theta * ( 7.0*theta - 3*v0 ) + 4.0 * theta * ( theta - v0 ) * kappa * T ) * exp( -kappa*T );
        double term_4 = ( 2.0*theta*( 2.0 *v0 - 3.0 *theta) + T * ( pow(gamma,2.0) - 2.0 * kappa * pow(theta-v0, 2.0) ) + pow(gamma,2.0)/(2.0*kappa) ) * exp( -2.0*kappa*T );
        double term_5 = ( 2.0*theta*(theta - v0) ) * exp( -3.0*kappa*T );
        double term_6 = 0.5 * ( pow(gamma,2.0) / ( 4.0 * kappa ) - pow( theta - v0, 2.0) ) * exp( -4.0*kappa*T );
        
        return term_1 * ( term_2 + term_3 + term_4 + term_5 + term_6 );
    }
    
    void update_Q_measure(double rho) {
        // Update the kappa_Q, theta_Q using the input rho;
        // !!!CHECK: if this depends on RHO before OR after changing sign
        kappa = kappa_base - gamma * rho; // new kappa under measure Q
        theta = kappa_base * theta_base / ( kappa_base - gamma * rho ); // new theta under measure Q
    };
    
public:
    SSB_option(const string& type, double s0, double K, double r, double gamma, double v0, int N, double kappa, double theta, double B):base_option(type, s0, K, r, gamma, v0, N),kappa(kappa),theta(theta),kappa_base(kappa),theta_base(theta) {
        model = "Stein Stein Barrier";
        
        // PAGE 61
        // calculate all the new values for the parameters
        double K_new = B / K; // new stirke
        double r_new = -r; // change the sign of r
        
        // reverse the rho_sign
        rho_sign = -1;
        
        // update the new values
        K = K_new;
        r = r_new;
    }
};

#endif /* Stein_Stein_Barrier_h */
