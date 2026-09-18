#include "abstract_integrator.hpp"
#include "rectangle_integrator.hpp"
#include "trapezoidal_integrator.hpp"
#include "gauss_integrator.hpp"
#include <iostream>
#include <iomanip>

double square(double x)
{
    return x*x;
}

int main()
{

    double a = 0;
    double b = 1;
    int Nint = 1000;
    int n_precision = 12;

    rectangle_integrator integ_rectangle(a,b,nullptr);
    trapezoidal_integrator integ_trapezoidal(a,b,nullptr);
    gauss_integrator integ_gauss(a,b,nullptr);

    integ_rectangle.change_integrand(square);
    integ_trapezoidal.change_integrand(square);
    integ_gauss.change_integrand(square);

    double rectangle = integ_rectangle.compute_integral(Nint);
    double trapezoidal = integ_trapezoidal.compute_integral(Nint);
    double gauss = integ_gauss.compute_integral(Nint);

    std::cout<< "Nint = " << Nint << std::endl;
    std::cout<< "Rectangle method : " << std::setprecision(n_precision) << rectangle << std::endl;
    std::cout<< "Trapezoidal method : "<< std::setprecision(n_precision) << trapezoidal << std::endl;
    std::cout<< "Gauss-Legendre quadrature : "<< std::setprecision(n_precision) << gauss << std::endl;
    return 0;
}