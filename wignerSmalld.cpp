
#include "eulerAngles.h"
#include <algorithm>
#include <cmath>
#include "main.h"


double factorial(double number) 
{
    return  tfactorial[int(number)];
}

double w(double &j, double &m, double &k, double &n, double &numerator)
{
    
    double denominator = factorial(j-m-n) * factorial(j+k-n) * factorial(n+m-k) * factorial(n);
    if(numerator > 1e300 || denominator > 1e300)
    {
        writeToErrorLogFile("Factorial in weigner function is too high");
    }
    double W = numerator/denominator;
    return W;
}


double wignerSmalld(double &j, double &m, double &k, eulerAngles &ea)
{
    double nmin = std::max(0.0,k-m);
    double nmax = std::min(j-m,j+k);
    double d=0.0;
    double numerator = sqrt(factorial(j+m)) * sqrt(factorial(j-m))
                                * sqrt(factorial(j+k)) * sqrt(factorial(j-k));
    
    for(double n=nmin;n<=nmax;++n)
    {
        double cosPow = double(2)*j + k - m - double(2)*double(n);
        double sinPow = m - k + double(2)*double(n);
        
        d += pow(-1,n) * w(j, m, k, n, numerator) * pow(cos(ea.beta/double(2)), cosPow ) * 
                pow(-sin(ea.beta/double(2)), sinPow );
    }
    return d;
}

