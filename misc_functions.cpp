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
using namespace std;

vector<double> linspace( double start, double end, int num_of_steps ) {
    vector<double> res;
    double step_size = (end - start) / ( num_of_steps );
    double cur_val = start;
    
    res.push_back(cur_val);
    for( int i = 0; i < (num_of_steps); i++ ) {
        cur_val += step_size;
        res.push_back( cur_val );
    }
    return res;
}

void print( vector<double> vec ) {
    for( unsigned int i = 0; i < vec.size(); i++ ) {
        cout << vec[i] << ", ";
    }
}

void print( vector<vector<double> > mat ) {
    for( unsigned int i = 0; i < mat.size(); i++ ) {
        for( unsigned int j = 0; j < mat[i].size(); j++ ) {
            cout << mat[i][j] << ", ";
        }
        cout << endl;
    }
}

void print( vector<vector<vector<double> > > mat ) {
    unsigned int k;
    for( unsigned int i = 0; i < mat.size(); i++ ) {
        for( unsigned int j = 0; j < mat[i].size(); j++ ) {
            cout << "(";
            for( k = 0; k < mat[i][j].size()-1; k++ ) {
                cout << mat[i][j][k] << ",";
            }
            cout << mat[i][j][k] << ") ";
        }
        cout << endl;
    }
}

vector<vector<vector<double> > > mesh_grid( double x_start, double x_end, int x_steps, double y_start, double y_end, int y_steps ) {
    vector<double> x_vec, y_vec;
    vector<vector<vector<double> > > res_mat;
    vector<vector<double> > row;
    vector<double> cell;
    
    x_vec = linspace( x_start, x_end, x_steps );
    y_vec = linspace( y_start, y_end, y_steps );
    
    for( int i = 0; i < y_vec.size(); i++ ) {
        for( int j = 0; j < x_vec.size(); j++ ) {
            cell.push_back( x_vec[j] );
            cell.push_back( y_vec[i] );
            row.push_back( cell );
            cell.clear();
        }
        res_mat.push_back(row);
        row.clear();
    }
    
    return res_mat;
}

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


