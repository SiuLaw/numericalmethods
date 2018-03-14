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
//#define M_PI 3.141592653589793238462643383279502884

using namespace std;

void hline(int n) {
    for(int i = 0; i < n; i++ ) {
        cout << "-";
    }
    cout << endl;
}

void dotline(int n) {
    for(int i = 0; i < n; i++ ) {
        cout << ".";
    }
    cout << endl;
}


void output_file(string file_name, vector<double> &vec ) {
    ofstream fout;
    fout.open (file_name);
    for(unsigned int i = 0; i < vec.size(); i++ ) {
        fout << vec[i] << endl;
    }
    fout.close();
}

void output_file(string file_name, vector< vector<double> > &mat ) {
    ofstream fout;
    fout.open (file_name);
    for(unsigned int i = 0; i < mat.size(); i++ ) {
        for(unsigned int j = 0; j < mat[i].size() - 1; j++ ) {
            fout << mat[i][j] << ",";
        }
        fout << mat[i][mat[i].size()-1];
        fout << endl;
    }
    fout.close();
}


//model independent string comparision
bool icompare_pred(unsigned char a, unsigned char b){
		return tolower(a) == tolower(b);
}
bool icompare(string const& a, string const& b){

		if (a.length()==b.length()) {
				return equal(b.begin(), b.end(),
													 a.begin(),icompare_pred);
		}
		else {
				return false;
		}
}

void res_print(vector<double> &res) {

    if(res.size() < 3) {
        cout << "   value    = " << res[0] << endl;
        cout << "   time     = " << res[1] << endl;
    } else {
        cout << "   estimate = " << res[0] << endl;
        cout << "   time     = " << res[1] << endl;
        cout << "   err_mean = " << res[2] << endl;
        cout << "   err_var  = " << res[3] << endl;
    }
}

double normalCDF(double value) {
	return 0.5 * erfc(-value / sqrt(2));
}

double normalPDF(double value) {
	return (1 / sqrt(2 * M_PI)) * exp(-0.5 * pow(value, 2));
}

void display_opening(string target_val, string method1, string method2 ) {
    hline(50);
    cout << target_val <<" using method: " << method1 << " with " << method2 << endl;
}

void display_opening(string target_val, string method1) {
    hline(50);
    cout << target_val <<" using method: " << method1 << endl;
}

void display_ending() {
    hline(50);
}

void display_section_split() {
    dotline(50);
}



/******************************************************************************/

//Normal random number generator
Normal normal(Custom);


//class
class vanilla_option{
    
    double s0;
    double K;
    double r;
    double mu;
    double gamma;
    double v0;
    bool is_call;
    
    unsigned int N; // number of time partitions
    

    // this function calculate pay off
    double pay_off(double S){
        if (is_call && S > K) {
            return S - K;
        } else if (!is_call && S < K) {
            return K - S;
        } else {
            return 0;
        }
    }
    
    double cal_Y1( double X1, double X2, double rho ) {
        return X1;
    }
    double cal_Y2( double X1, double X2, double rho ) {
        return rho * X1 + pow(1-rho*rho,0.5) * X2;
    }

    void HW_fd_pricing(double rho, double T, vector<double>& res) {
        // prameter initlisation
        clock_t c;
        double duration;
        double dt = (double)T / N; // calculate fix value
        double dt_sqrt = pow(dt, 0.5); // calculate fix value
        double discount = exp(-r * T); // discouting factor
        
        double this_price, this_price_anti; // itermediate price results
        double price_sum = 0;
        double price_sum_sqr = 0;
        double price;
        
        
        double s = s0; // intermediate price
        double v = v0; // intermediate volatility
        
        // anti variate
        double s_anti = s0;
        double v_anti = v0;
        
        // run the simulation, NO variance reduction
        
        
        double tol = 1/pow(10.0,4); //
        double var = 1; // inital var
        int M = 0; // counter for
        double W, Z; // the values of the correlated Normal samples
        
        c = clock();
        while( var > tol ) {
            
            M += 1;
            
            vector<double> X1 = normal.generate(N); //generate normal vector of size N
            vector<double> X2 = normal.generate(N); //generate normal vector of size N;
            
            for (unsigned int j = 0;j < N;j++) {
                W = cal_Y1( X1[j], X2[j], rho);
                Z = cal_Y2( X1[j], X2[j], rho);
                
                v += mu * v * dt + gamma * v * W * dt_sqrt;
                s += r * s * dt + pow(v,0.5) * s * Z * dt_sqrt;
                
                // v_anti += mu * v_anti * dt + gamma * v_anti * -W * dt_sqrt;
                // s_anti += r * s_anti * dt + pow(v,0.5) * s_anti * -Z * dt_sqrt;
            }
            
            this_price = pay_off( s );
            // this_price_anti = pay_off( s_anti );
            
            price_sum     += this_price;
            // price_sum     += this_price_anti;
            price_sum_sqr += pow( this_price, 2 );
            // price_sum_sqr += pow( this_price_anti, 2 );
            
            s = s0; // rets path
            v = v0; // rets path
            
            // s_anti = s0; // rets path
            // v_anti = v0; // rets path
            
            if( M % 1000 == 0 ) {
                var = ( price_sum_sqr / M  - pow( price_sum / M , 2 ) ) / M;
                tol = discount * price_sum / M / 100 / 1000;
                cout << "M = " << M << "; " << (int)(tol/var * 100)  << "% ; price_est = " << discount * price_sum / M  << endl;
            }
        }
        
        cout << "tol = " << tol << " reached after M = " << M << endl;
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        price = (discount * price_sum / M );

        res.push_back(price); // price
        res.push_back(duration); //time
    }

public:
    // object constructor
    vanilla_option(const string& type, double s0, double K, double r, double mu, double gamma, double v0, int N):s0(s0),K(K),r(r),mu(mu),gamma(gamma),v0(v0), N(N) {
        if(icompare(type,"call")){
            is_call = true;
        } else if(icompare(type,"put")){
            is_call = false;
        } else {
            cout << "call/put string not understood, default to call" << endl;
            is_call = true;
        }
    }

    // this function calculate the option price at time 0 and analytic statistics
    vector<double> calculate_price(const string& method, double rho = -0.5, double T = 1.0, bool display_results = false) {
        
        // return containeer
        vector<double> res;

        //Use selected method to calcuate c_0
        if (icompare(method, "HW_fd")) {
            HW_fd_pricing(rho, T, res);
        } else if (icompare(method,"analytic")){
            // analytic_solution_pricing(S0, r, v, res);
        } else {
            cout << "method " << method << " not recognized, please choose another" << endl;
        }

        string name = "PRICE";
        if(display_results) {
            display_opening(name,method);

            display_section_split();
            
            cout << "   T        = " << T << endl;
            
            cout << "   s0       = " << s0 << endl;
            cout << "   K        = " << K << endl;
            cout << "   r        = " << r  << endl;
            cout << "   mu       = " << mu  << endl;
            cout << "   gamma    = " << gamma  << endl;
            cout << "   v0       = " << v0  << endl;

            if( !icompare(method,"analytic") ) {
                cout << "   N        = " << N << endl;
            }
            display_section_split();
            res_print(res);
            display_ending();
        }

        return res;

    }
};


int main() {
    double s0 = 100;
    double K = 97;
    double r = 0.01;
    double mu = 0.2;
    double gamma = 0.1;
    double v0 = 0.04;
    
    int N = 10000;
    
    double rho = -0.5;
    double T = 0.125;
    
    // create an vanilla option
    vanilla_option opt("put",s0,K,r,mu,gamma,v0,N);
    
    opt.calculate_price("HW_fd", rho, T, true);

	return 0;
}
