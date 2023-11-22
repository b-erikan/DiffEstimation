#include "AlgDiff.h"
#include <iostream>
AlgDiff::AlgDiff(float ts, float alpha, float beta, int N, float T, float wc, bool corr)
    : __ts{ts}, __alpha{alpha}, __beta{beta}, __N{N}, correction{corr}, __wc{wc}
{
    if (0 == __N)
    {
        __theta = 0;
        __thetaBool = false;
    }
    else
    {
        // Calculate roots of jacobi here later and get max value as theta.
        __thetaBool = true;
    }

    // Get corner frequency if window length is given. Get window length if corner frequency is given.

    if (T)
    {
        __T = T;
        __wc = get_cutoffFreq();
        __setUp = 'T';
    }
    else if (wc)
    {
        __wc = wc;
        computeTfromWc(__wc);
        __setUp = 'w';
    }
    else
        throw std::invalid_argument("T or wc must be given."); // Maybe unnecessary to check?
}

float AlgDiff::get_delay()
{
    if (0 == __N)
    {
        delay = (__alpha + 1) / (__alpha + __beta + 2) * __T;
    }
    else
    {
        delay = (1 - __theta) / 2 * __T;
    }
    return delay;
}

float AlgDiff::get_cutoffFreq(void)
{
    // Calculate wc based on T here.
    return 0;
}
void AlgDiff::computeTfromWc(float wc)
{
    // Calculate T based on wc here.
    return;
}

void AlgDiff::discretize(int der, bool reduceFilLength, float redTol, bool discreteSpectrum, method mtd)
{
    int L0 = static_cast<int>(__T / __ts);
}

std::vector<double> AlgDiff::weightFcn(double a, double b, std::vector<double> &t)
{
    std::vector<double> w(t.size());

    if (a < 0 || b < 0)
    {
        for (size_t i = 0; i < t.size(); ++i)
        {
            if (t[i] >= 1.0)
                t[i] = 1.0;
            if (t[i] <= -1.0)
                t[i] = -1.0;

            w[i] = std::pow(1.0 - t[i], a) * std::pow(t[i] + 1.0, b);

            if (t[i] >= 1.0 || t[i] <= -1.0)
                w[i] = 0.0;
        }
    }
    else
    {
        for (size_t i = 0; i < t.size(); ++i)
        {
            if (t[i] > 1.0)
                t[i] = 1.1;
            if (t[i] < -1.0)
                t[i] = -1.1;

            w[i] = std::pow(1.0 - t[i], a) * std::pow(t[i] + 1.0, b);

            if (t[i] > 1.0 || t[i] < -1.0)
                w[i] = 0.0;
        }
    }

    return w;
}

std::vector<float> AlgDiff::newton_cotes_rules(const std::vector<float> &p, int order, int L)
{
    // We don't have numpy arrays so I used vectors instead. May not need it.
    order -= 1;
    std::vector<float> out(L, 0.0);
    out[0] = p[0];

    for (int i = 1; i < L - 1; ++i)
    {
        if (i % order == 0)
        {
            out[i] = p[0] + p.back();
        }
        else
        {
            out[i] = p[i % order];
        }
    }

    out.back() = p.back();
    return out;
}
