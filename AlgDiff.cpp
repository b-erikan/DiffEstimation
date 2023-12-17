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

std::vector<float> AlgDiff::timeShift(std::vector<float> &t)
{
    std::vector<float> result;
    result.reserve(t.size());

    for (float val : t)
    {
        result.push_back(1.0 - 2.0 / __T * val);
    }

    return result;
}
std::vector<float> AlgDiff::specialConvolve(const std::vector<float> &x)
{
    size_t wSize = __w.size();
    size_t xSize = x.size();
    size_t outputSize = wSize + xSize - 1;
    std::vector<float> result;
    float sum = 0;

    for (int i = 0; i < outputSize; i++)
    {
        int iter = i;
        for (int j = 0; j <= i; j++)
        {
            if (xSize <= j || wSize <= iter)
                sum;
            else
            {
                sum += x.at(j) * __w.at(iter);
            }
            iter--;
        }
        result.push_back(sum);
        sum = 0;
    }
    std::vector<float> output(result.cbegin() + (wSize / 2), result.cend() - (wSize / 2));
    for (int i = 0; i < wSize / 2; i++)
    {
        output.insert(output.begin(), 0);
        output.pop_back();
    }
    return output;
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

    switch (mtd)
    {
    case mid_point:

    case euler:
        break;
    case trapezoidal:
        break;
    case simpson_rule:
        break;
    case simpson_38_rule:
        break;
    case boole_rule:
        break;
    }
}

std::vector<float> AlgDiff::weightFcn(float a, float b, std::vector<float> &t)
{
    std::vector<float> w(t.size());

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

std::vector<float> AlgDiff::evalKernelDer(std::vector<float> &t, int k)
{
    // This function needs work!!! Totally wrong...
    float a = __alpha;
    float b = __beta;
    // float dg = 0;
    std::vector<float> shifted_t = timeShift(t);
    std::vector<float> w = weightFcn(a - k, b - k, shifted_t);
    std::vector<float> out;
    for (int i = 0; i <= __N; ++i)
    {
        double h = pow(2, a + b + 1) * tgamma(i + a + 1) * tgamma(i + b + 1) / (tgamma(i) * (2 * i + a + b + 1) * tgamma(i + a + b + 1));
        std::vector<float> P;
        for (float var : shifted_t)
        {
            P.push_back(boost::math::jacobi(i + k, a - k, b - k, var));
        }
        out = P;
        // dg += 1 / h * boost::math::jacobi(i, a, b, __theta) * P * tgamma(i + 1 + k)/tgamma(i + 1);
    }

    return out;
}
std::vector<float> AlgDiff::estimateDer(const unsigned der, const std::vector<float> &x)
{
    
    if (0 == der)
    {
        __w = {0.00209468, 0.08364409, 0.28054109, 0.38788526, 0.27344734,
               0.05872276, -0.05375467, -0.03154697, -0.00103358};
    }
    else if (1 == der)
    {
        __w = {15.30847194, 158.49429167, 191.74630732, -3.36629456,
               -200.59170407, -188.09790359, -30.50232806, 47.87577424,
               7.35648253};
    }

    return specialConvolve(x);
    // return convWithWeights(x);
}
/*
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
}*/
