#pragma once

typedef double (*pfn) (double);


class abstract_integrator
{
    public:

    abstract_integrator(double a, double b, pfn f);
    abstract_integrator(abstract_integrator const&);

    double lwbd() const {return lower_;};
    double upbd() const {return upper_;};
    void change_bounds(double a, double b);
    void change_integrand(pfn f){integrand_ = f;};
    virtual double compute_integral(int Nint) const = 0;

    virtual ~abstract_integrator(){};


    protected:

    double lower_;
    double upper_;
    pfn integrand_;
};

abstract_integrator::abstract_integrator(double a, double b, pfn f)
{
    lower_ = a;
    upper_ = b;
    integrand_ = f;
}

abstract_integrator::abstract_integrator(abstract_integrator const& integ)
{
    integrand_ = integ.integrand_;
    lower_ = integ.lower_;
    upper_ = integ.upper_;
}

void abstract_integrator::change_bounds(double a, double b)
{
    lower_ = a;
    upper_ = b;
}
