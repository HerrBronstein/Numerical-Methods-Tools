#include "abstract_integrator.hpp"

class rectangle_integrator : public abstract_integrator
{
    public:

    rectangle_integrator(double a, double b, pfn f):abstract_integrator(a,b,f){};
    rectangle_integrator(rectangle_integrator const& integ):abstract_integrator(integ){};
    ~rectangle_integrator(){};


    double compute_integral(int Nint) const;
};


double rectangle_integrator::compute_integral(int Nint) const
{
    double sum{0.};
    double h = (upper_ - lower_)/Nint;

    for(int ii=0; ii<Nint; ii++)
    {
        sum += integrand_(lower_ + ii*h);
    }

    return h*sum;
}