#include "abstract_integrator.hpp"

class trapezoidal_integrator : public abstract_integrator
{
    public:

    trapezoidal_integrator(double a, double b, pfn f):abstract_integrator(a,b,f){};
    trapezoidal_integrator(trapezoidal_integrator const& integ):abstract_integrator(integ){};
    ~trapezoidal_integrator(){};

    double compute_integral(int Nint) const;

};


double trapezoidal_integrator::compute_integral(int Nint) const
{
    double sum{0.};
    double h = (upper_ - lower_)/Nint;

    for(int ii=0; ii<Nint; ii++)
    {
        sum += 0.5*(integrand_(lower_ + ii*h) + integrand_(lower_ + (ii+1)*h));
    }

    return h*sum;
}