class AlgDiff{
    AlgDiff::AlgDiff(float ts, float alpha,float beta,int N,float T,float wc, bool corr);
private:
    float __ts,__alpha,__beta,__T,__theta,__wc;
    int __N;
    bool correction,__thetaBool;
    char __setUp;
};