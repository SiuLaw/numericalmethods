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
