
#ifndef SIMULATIONMETHODS_NORMAL_H
#define SIMULATIONMETHODS_NORMAL_H

#include <random>
#include <chrono>
#include "Random.h"

using namespace std;

class Normal: public Random {
private:
    default_random_engine generator;
protected:
    double m_Mean;
    double m_Variance;
    vector<double> custom_generate(unsigned int n){
        double v1;
        double v2;
        double w;
        vector<double> v;
        for(unsigned int i=0;i<(n/2)+1;i++){
            v1 = 2.0*rand()/RAND_MAX -1;
            v2 = 2.0*rand()/RAND_MAX -1;
            w = pow(v1,2)+pow(v2,2);
            if(w<=1){
                v.push_back(sqrt(-2*log(w)/w)*v1);
                v.push_back(sqrt(-2*log(w)/w)*v2);
            } else {
                --i;
            }
        }
        //size adjustment
        if(v.size()>n+1) {
            v.pop_back();
            v.pop_back();
        } else {
            v.pop_back();
        }
        
        return move(v);
    }
    
    vector<double> standard_generate(unsigned int n) {
        int j;
        normal_distribution<double> distribution(m_Mean, m_Variance);
        vector<double> rand_numbers;
        for (j = 0; j < n; j = j + 1) {
            double random_number = distribution(generator);
            rand_numbers.push_back(random_number);
        }
        return move(rand_numbers);
    }
public:
    Normal(const GeneratorType generatorType = Custom, int seed = std::chrono::system_clock::now().time_since_epoch().count(),  double mean=0.0, double variance=1.0): Random(generatorType), generator(seed) {
        this->m_Mean = mean;
        this->m_Variance = variance;
    }
    vector<double> generate(unsigned int n) {
        //This function defines which implementation to use - a custom one or a c++ standard based one
        if(m_GeneratorType==Custom) {
            return custom_generate(n);
        } else {
            return standard_generate(n);
        }
    }
};

#endif //SIMULATIONMETHODS_NORMAL_H
