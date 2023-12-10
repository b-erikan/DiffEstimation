#pragma once
#include <vector>
#include <cmath>
#include "/home/berke/Downloads/boost_1_83_0/boost/math/special_functions/jacobi.hpp" //Careful with this!!
/*
#define MID_POINT 0
#define EULER 1
#define TRAPEZOIDAL 2
#define SIMPSON_RULE 3
#define SIMPSON_38_RULE 4
#define BOOLE_RULE 5
*/
enum method {mid_point,euler,trapezoidal,simpson_rule,simpson_38_rule,boole_rule};


class AlgDiff{
public:
    AlgDiff(float ts, float alpha,float beta,int N,float T,float wc, bool corr);
    float get_cutoffFreq(void);
    void computeTfromWc(float wc);
    void discretize(int der, bool reduceFilLength, float redTol, bool discreteSpectrum,method mtd);
    float get_delay(void);
    std::vector<float> evalKernelDer(std::vector<float> &t, int k);
    std::vector<float> estimateDer(const int der,const std::vector<float> &x);
    
private:
    
    std::vector<float> convolution(const std::vector<float>& x,const std::vector<float>& weights,const std::string& conv_type);
    //std::vector<float> newton_cotes_rules(const std::vector<float>& p, int order, int L);
    std::vector<float> weightFcn(float a, float b, std::vector<float> &t);
    std::vector<float> timeShift(std::vector<float> &t);
    std::vector<float> __w;
    float __ts,__alpha,__beta,__T,__theta,__wc,delay;
    int __N;
    bool correction,__thetaBool;
    char __setUp;
};