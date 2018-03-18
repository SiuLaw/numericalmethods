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
#include "Stein_Stein.h"
#include "Stein_Stein_Barrier.h"
//#define M_PI 3.141592653589793238462643383279502884

using namespace std;

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
