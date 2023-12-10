#include <iostream>
#include "AlgDiff.h"
#include "AlgDiff.cpp"
int main(void){
    /*This is a testbench for the AlgDiff class.*/


    //double val = tgamma(0.32551); Gamma works.
    float ts=0.001, T=0.009, alpha=4.0, beta=4.0, wc=1000;
    int N=0;
    bool corr=false;
    AlgDiff diff(ts,alpha,beta,N,T,wc,corr);


    std::vector<float> x={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    
    std::vector<float> dx=diff.estimateDer(1,x);

    for(float x:dx){
        std::cout<<x<<std::endl;
    } 
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