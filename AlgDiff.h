#pragma once
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
    AlgDiff::AlgDiff(float ts, float alpha,float beta,int N,float T,float wc, bool corr);
    float get_cutoffFreq(void);
    void computeTfromWc(float wc);
    void discretize(int der, bool reduceFilLength, float redTol, bool discreteSpectrum,method mtd);
private:
    float __ts,__alpha,__beta,__T,__theta,__wc;
    int __N;
    bool correction,__thetaBool;
    char __setUp;
};