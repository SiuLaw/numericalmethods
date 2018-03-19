#ifndef base_option_h
#define base_option_h

//base class
class base_option{
protected:
    double s0;
    double K;
    double r;
    double gamma;
    double v0;
    bool is_call;
    
    double rho_sign;
    
    Normal normal;
    
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
        return rho * rho_sign * X1 + pow(1-rho*rho,0.5) * X2;
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
    
    virtual double v_increment(double v, double Z, double dt, double dt_sqrt) { return 0; }
    virtual double s_increment(double s, double W, double dt, double dt_sqrt, double v) { return 0; }
    
    void fd_pricing(double rho, double T, vector<double>& res, bool display_flag) {
        // prameter initlisation
        clock_t c;
        double duration;
        double dt = (double)T / N; // calculate fix value
        double dt_sqrt = pow(dt, 0.5); // calculate fix value
        double discount = exp(-r * T); // discouting factor
        
        update_Q_measure(rho); // update rho if needed
        
        double this_price; // itermediate price results
        double price_sum = 0;
        double price_sum_sqr = 0;
        double price;
        
        int percent;
        int freq_display = 20;
        
        
        double s = s0; // intermediate price
        double v = v0; // intermediate volatility
        
        
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
    
    virtual double NW(double T) {return 0;}
    virtual double NN(double T) {return 0;}
    
    virtual void update_Q_measure(double rho) {}; // Update the kappa_Q, theta_Q using the input rho;
    
    void ap_pricing(double rho, double T, vector<double> &res ) {
        // prameter initlisation
        clock_t c;
        double duration;
        double price;
        
        c = clock();
        
        update_Q_measure(rho);
        
        double v0_est = vt(T);
        double v0_fut = v0_est;
        
        // cout << endl;
        // cout << "v0 = " << v0 << endl;
        // cout << "v0_est = " << v0_est << endl;
        
        double P_BS = cal_P_BS(T, v0_fut );
        
        double d_plus = ( log(s0 / K ) + ( r + pow(v0_fut,2.0) / 2.0 ) * T ) / ( v0_fut * pow(T, 0.5) );
        double d_mins = d_plus - v0_fut * pow( T, 0.5 );
        
        double first_bracket  = rho*rho_sign/2.0 * H(T,v0_fut,d_plus,d_mins) * NW(T);
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
        
        rho_sign = 1; // by default, the rho is not reversed;
        Normal new_normal(Custom);
        normal = new_normal;
    }
    
    void set_gamma(double new_gamma ) {
        gamma = new_gamma;
    }
    
    void set_K( double new_K ) {
        K = new_K;
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
            cout << "   rho      = " << rho*rho_sign << endl;
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


#endif
