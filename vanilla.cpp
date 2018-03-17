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


//base class
class base_option{
protected:
    double s0;
    double K;
    double r;
    double gamma;
    double v0;
    bool is_call;
    
    string model;
    
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
    
    double cal_P_BS(double T, double v0_fut) {
        double d1 = ( log(s0 / K) + ( r + pow(v0_fut,2.0) / 2.0 ) * T ) / ( v0_fut * pow(T,0.5) );
        double d2 = d1 - v0_fut * pow( T,0.5);
        
        // cout << "d1 = " << d1 << endl;
        // cout << "d2 = " << d2 << endl;
        
        return normalCDF(-d2) * K * exp( -r*T) - s0 * normalCDF(-d1);
    }
    
    double H(double T, double v0_fut, double d_plus, double d_mins ) {
        // Since H is model independent, define as a function
        double X_0 = log( s0 );
        return exp( X_0 ) / ( pow( v0_fut, 2.0 ) * T * pow( 2.0 * M_PI, 0.5 )  ) * exp( - pow(d_plus, 2.0)/2.0 ) * ( - d_mins );
    }
    
    double J(double T, double v0_fut, double d_plus, double d_mins ) {
        // Since J is model independent, define as a function
        double X_0 = log( s0 );
        return exp( X_0 ) / ( pow( v0_fut * pow( T, 0.5 ), 3.0 ) * pow( 2.0 * M_PI, 0.5 ) ) * exp( - pow(d_plus, 2.0)/2.0 ) * ( d_plus * d_mins - 1.0 );
    }
    
    virtual void print_model_para() {
        cout << "   ... no additional model para specified" << endl;
    }
    
    virtual double v_increment(double v, double Z, double dt, double dt_sqrt) {};
    virtual double s_increment(double s, double W, double dt, double dt_sqrt, double v) {};
    
    void fd_pricing(double rho, double T, vector<double>& res, bool display_flag) {
        // prameter initlisation
        clock_t c;
        double duration;
        double dt = (double)T / N; // calculate fix value
        double dt_sqrt = pow(dt, 0.5); // calculate fix value
        double discount = exp(-r * T); // discouting factor
        
        double this_price; // itermediate price results
        double price_sum = 0;
        double price_sum_sqr = 0;
        double price;
        
        int percent;
        int freq_display = 20;
        
        
        double s = s0; // intermediate price
        double v = v0; // intermediate volatility
        
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
                
                v += v_increment(v, W, dt, dt_sqrt);
                s += s_increment(s, Z, dt, dt_sqrt, v);
            }
            
            this_price = pay_off( s );
            
            price_sum     += this_price;
            price_sum_sqr += pow( this_price, 2 );
            
            s = s0; // rets path
            v = v0; // rets path
            
            if( M % 1000 == 0 ) {
                var = ( price_sum_sqr / M  - pow( price_sum / M , 2 ) ) / M;
                tol = discount * price_sum / M / 100 / 1000;
                percent = (int)(tol/var * 100);
                if(display_flag && percent % freq_display == 0 ) {
                    cout << "M = " << M << "; " << percent << "% ; price_est = " << discount * price_sum / M  << endl;
                }
                
            }
        }
        
        if(display_flag) {
            cout << "tol = " << tol << " reached after M = " << M << endl;
        }
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        price = (discount * price_sum / M );
        
        res.push_back(price); // price
        res.push_back(duration); //time
    }
    
    double vt(double T) {
        // prameter initlisation
        clock_t c;
        double duration;
        double dt = (double)T / N; // calculate fix value
        double dt_sqrt = pow(dt, 0.5); // calculate fix value
        
        double this_price; // itermediate price results
        double vol_sum = 0;
        double vol_sum_sqr = 0;
        double vol;
        double this_vol_sum = 0;
        
        double v = v0;
        
        double tol = 1/pow(10.0,10); //
        double var = 1; // inital var
        int M = 0; // counter for
        double W, Z; // the values of the correlated Normal samples
        
        c = clock();
        while( var > tol ) {
            M += 1;
            
            vector<double> X1 = normal.generate(N); //generate normal vector of size N
            
            for (unsigned int j = 0;j < N;j++) {
                W = X1[j];
                v += v_increment(v, W, dt, dt_sqrt);
                this_vol_sum += pow(v,2.0);
            }
            vol_sum     += this_vol_sum;
            vol_sum_sqr += pow( this_vol_sum, 2 );
            
            // reset parameters for next loop
            v = v0; // resets initial v0
            this_vol_sum = 0; // resets sum of volatility
            
            if( M % 1000 == 0 ) {
                var = ( vol_sum_sqr / M  - pow( vol_sum / M , 2 ) ) / M;
                tol = vol_sum / M / 100 / 1000;
                // cout << "M = " << M << "; " << (int)(tol/var * 100)  << "% ; future_average_vol = " << pow( vol_sum / M / N , 0.5 ) << endl;
            }
        }
        
        // cout << "tol = " << tol << " reached after M = " << M << endl;
        
        // record time
        return pow( vol_sum / M / N , 0.5 );
        
    }
    
    virtual double NW(double T) {};
    virtual double NN(double T) {};
    
    void ap_pricing(double rho, double T, vector<double> &res ) {
        // prameter initlisation
        clock_t c;
        double duration;
        double price;
        
        c = clock();
        
        double v0_est = vt(T);
        double v0_fut = v0_est;
        
        // cout << endl;
        // cout << "v0 = " << v0 << endl;
        // cout << "v0_est = " << v0_est << endl;
        
        double P_BS = cal_P_BS(T, v0_fut );
        
        double d_plus = ( log(s0 / K ) + ( r + pow(v0_fut,2.0) / 2.0 ) * T ) / ( v0_fut * pow(T, 0.5) );
        double d_mins = d_plus - v0_fut * pow( T, 0.5 );
        
        double first_bracket  = rho/2.0 * H(T,v0_fut,d_plus,d_mins) * NW(T);
        double second_bracket = 1.0/8.0 * J(T,v0_fut,d_plus,d_mins) * NN(T);
        
        // record time
        duration = (clock() - c) / (double)CLOCKS_PER_SEC;
        price = P_BS + first_bracket + second_bracket;
        
        res.push_back(price); // price
        res.push_back(duration); //time
        
        // cout << "P_BS = " << P_BS << endl;
        
    }

public:
    base_option(const string& type, double s0, double K, double r, double gamma, double v0, int N):s0(s0),K(K),r(r),gamma(gamma),v0(v0), N(N) {
        if(icompare(type,"call")){
            is_call = true;
        } else if(icompare(type,"put")){
            is_call = false;
        } else {
            cout << "call/put string not understood, default to call" << endl;
            is_call = true;
        }
    }
    
    vector<double> calculate_price(const string& method, double rho = -0.5, double T = 1.0, bool display_results = false) {
        
        // return containeer
        vector<double> res;
        
        //Use selected method to calcuate c_0
        if (icompare(method, "finite_diff")) {
            fd_pricing(rho, T, res, display_results);
        } else if (icompare(method,"decomp_approx")){
            ap_pricing(rho, T, res);
        } else {
            cout << "ERROR: method " << method << " not recognized, please choose another" << endl;
            display_results = false;
        }
        
        string name = "PRICE";
        if(display_results) {
            display_opening(name,method,model);
            
            display_section_split();
            
            cout << "   T        = " << T << endl;
            
            cout << "   s0       = " << s0 << endl;
            cout << "   K        = " << K << endl;
            cout << "   r        = " << r  << endl;
            cout << "   gamma    = " << gamma  << endl;
            cout << "   v0       = " << v0  << endl;
            cout << "   rho      = " << rho << endl;
            if( !icompare(method,"ap") ) {
                cout << "   N        = " << N << endl;
            }
            print_model_para(); // additional parameters specific to the model
            
            display_section_split();
            res_print(res);
            display_ending();
        }
        
        return res;
        
    }
};


/************************************************************************************************************************************************************************/

// Hull White model
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


/************************************************************************************************************************************************************************/

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

int main() {
    double s0 = 100;
    double K = 97;
    double r = 0.01;
    double mu = 0.2;
    double gamma = 0.1;
    double v0 = 0.04;
    
    int N = 10;
    
    double rho = -0.5;
    double T = 0.125;
    
    // create an option under Hull White model
    // PAGE 46
    HW_option HW_opt("put",s0,K,r,gamma,v0,N,mu);
    HW_opt.calculate_price("finite_diff", rho, T, true);
    HW_opt.calculate_price("decomp_approx", rho, T, true);
    
    double kappa = 4;
    double theta = 0.2;
    v0 = 0.2;
    
    // create an option under Stein Stein model
    // PAGE 49
    SS_option SS_opt("put",s0,K,r,gamma,v0,N,kappa,theta);
    SS_opt.calculate_price("finite_diff",rho,T,true);
    SS_opt.calculate_price("decomp_approx",rho,T,true);

	return 0;
}
