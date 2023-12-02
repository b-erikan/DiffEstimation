#include <iostream>
#include "AlgDiff.h"
#include "AlgDiff.cpp"
int main(void){
    /*This is a testbench for the AlgDiff class. When it's done, I will give the same inputs to both algebraicDiferentiator.py and AlgDiff.cpp and compare the outputs.*/


    //double val = tgamma(0.32551);
    float ts=0.001, T=0.068, alpha=4.0, beta=4.0, wc=100;
    int N=0;
    bool corr=false;
    AlgDiff diff(ts,alpha,beta,N,T,wc,corr);


    std::vector<float> t={0.2,0.3,0.4,0.5};
    /*

    timeShift works.

    std::vector<float> shifted_t=diff.timeShift(t);

    for(float x:shifted_t){
        std::cout<<x<<std::endl;
    } 
    
    weightFcn works.

    std::vector<float> w=diff.weightFcn(alpha,beta,t);
    for(float x:w){
        std::cout<<x<<std::endl;
    } 

    */
    
    



    return 0;
}