#include "AlgDiff.h"

AlgDiff::AlgDiff(float ts, float alpha,float beta,int N,float T,float wc, bool corr)
        :__ts {ts}, __alpha{alpha}, __beta{beta}, __N {N}, correction {corr}, __wc {wc}
    {
        if(0==__N){
            __theta=0;
            __thetaBool=false;
        }
        else{
            //Calculate roots of jacobi here later and get max value as theta.
            __thetaBool=true;
        }

        //Get corner frequency if window length is given. Get window length if corner frequency is given.

        if (T){
            __T=T;
            __wc=get_cutoffFreq();
            __setUp='T';
        }

    }
    
    float get_cutoffFreq(void){
        //Calculate wc based on T here.
        return 0;
    }