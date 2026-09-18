#include "abstract_integrator.hpp"
#include <cmath>

class gauss_integrator : public abstract_integrator
{
    public:

    gauss_integrator(double a, double b, pfn f):abstract_integrator(a,b,f){};
    gauss_integrator(gauss_integrator const& integ):abstract_integrator(integ){};
    ~gauss_integrator(){};

    double compute_integral(int Nint) const;

};

double gauss_integrator::compute_integral(int Nint) const
{
    double sum{0.};
    double x[5];
    double w[5];

    x[0] = 0;
    x[1] = 1/3. *std::sqrt(5 - 2*std::sqrt(10/7.));
    x[2] = -1/3. *std::sqrt(5 - 2*std::sqrt(10/7.));
    x[3] = 1/3. *std::sqrt(5 + 2*std::sqrt(10/7.));
    x[4] = -1/3. *std::sqrt(5 + 2*std::sqrt(10/7.));

    w[0] = 128/225.0;
    w[1] = (322 + 13*std::sqrt(70))/900.0;
    w[2] = (322 + 13*std::sqrt(70))/900.0;
    w[3] = (322 - 13*std::sqrt(70))/900.0;
    w[4] = (322 - 13*std::sqrt(70))/900.0;
    
    for(int ii=0; ii<5; ii++)
    {
        sum += w[ii]*integrand_((upper_ - lower_)/2.0 * x[ii] + (upper_ + lower_)/2.0 );
    }
    
    return 0.5*(upper_ - lower_)*sum;
}