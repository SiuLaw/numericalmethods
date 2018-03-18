#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>
#include "Normal.h"
#include "misc_functions.h"
#include "base_option.h"
#include "Hull_White.h"
//#define M_PI 3.141592653589793238462643383279502884

using namespace std;

// Stein Stein model
class SS_option : public base_option{
protected:
    double kappa;
    double theta;
    
    void print_model_para() {
        //  cout << "            = " << << endl;
        cout << "   kappa    = " << kappa << endl;
        cout << "   theta    = " << theta << endl;
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
public:
    SS_option(const string& type, double s0, double K, double r, double gamma, double v0, int N, double kappa, double theta):base_option(type, s0, K, r, gamma, v0, N),kappa(kappa),theta(theta) {
        // testing
        model = "Stein Stein";
    }
};


/************************************************************************************************************************************************************************/

// Stein Stein BARRIER model
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


/************************************************************************************************************************************************************************/



int main() {
    // parameters defind;
    double s0,K,r,v0; // these are variables used by both models;
    double rho,T,gamma; // these are variables that will be varied
    double mu,kappa,theta; // theres are model specific variables;
    
    double B; // barrier value for later models;
    
    int N = 10; // number of path steps in the simulation
    
    // PAGE 46
    // create an put option under Hull White model
    s0 = 100;
    K = 97;
    r = 0.01;
    gamma = 0.1;
    v0 = 0.04;
    
    rho = -0.5;
    T = 0.125;
    
    mu = 0.2;
     
    HW_option HW_opt("put",s0,K,r,gamma,v0,N,mu);
    HW_opt.calculate_price("finite_diff", rho, T, true);
    HW_opt.calculate_price("decomp_approx", rho, T, true);
     
    
    /*
    
    // PAGE 49
    // create an put option under Stein Stein model
    s0 = 100;
    K = 97;
    r = 0.01;
    gamma = 0.1;
    v0 = 0.2;
    
    T = 0.5;
    rho = -0.5;
    
    kappa = 4;
    theta = 0.2;
    
    SS_option SS_opt("put",s0,K,r,gamma,v0,N,kappa,theta);
    SS_opt.calculate_price("finite_diff",rho,T,false);
    SS_opt.calculate_price("decomp_approx",rho,T,true);
    
    
    
    
    
    // create an barrier option under Stein Stein model
    // PAGE 62
     B = 95;
     r = 0;
     theta = 0.04;
     T = 0.5;
     
    SSB_option SSB_opt("put",s0,K,r,gamma,v0,N,kappa,theta,B);
    SSB_opt.calculate_price("finite_diff",rho,T,false);
    SSB_opt.calculate_price("decomp_approx",rho,T,true);
    

    
    /*
    // PAGE 52;
    s0 = 100;
    K = 97;
    r = 0.01;
    v0 = 0.2;
    T = 0.5;
    
    gamma = 0.2;
    
    kappa = 4;
    theta = 0.2;
    
    vector<vector<vector<double> > > mg;
    // mesh_grid( rho_start, rho_end, rho_steps,  gamma_start, gamma_end, gamma_steps )
    mg = mesh_grid( 0,-1,2, 0.1,0.4,2 );
    print(mg);
    
    vector<vector<double> > fd_mat, ap_mat;
    vector<double> fd_vec, ap_vec;
    vector<double> fd_res, ap_res ;
    SS_option SS_opt("put",s0,K,r,gamma,v0,N,kappa,theta);
    
    for( unsigned int i = 0; i < mg.size(); i++ ) {
        for( unsigned int j = 0; j < mg[i].size(); j++ ) {
            rho = mg[i][j][0];
            gamma = mg[i][j][1];
            cout << "(" << rho << "," << gamma << ") " << endl;;
            
            SS_opt.set_gamma( gamma );
            fd_res = SS_opt.calculate_price("finite_diff",rho,T,false);
            ap_res = SS_opt.calculate_price("decomp_approx",rho,T,false);
            
            fd_vec.push_back( fd_res[0] );
            ap_vec.push_back( ap_res[0] );
        }
        cout << endl;
        
        fd_mat.push_back( fd_vec );
        ap_mat.push_back( ap_vec );
        
        fd_vec.clear();
        ap_vec.clear();
        
    }
    cout << endl;
    
    cout << "(rho,gamma) = " << endl;
    print(mg);
    cout << endl;
    
    cout << "put value = " << endl;
    print(fd_mat);
    cout << endl;
    
    cout << "approx value = " << endl;
    print(ap_mat);
    cout << endl;
     */
    
	return 0;
}
