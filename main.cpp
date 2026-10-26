#include <iostream>
#include <cmath>
#include <iomanip>
#include "time_scheme.hpp"
#include "explicit_euler.hpp"
#include "RK2_scheme.hpp"
#include "RK4_scheme.hpp"

double f(double t, double u)
{
    return u;
}

int main()
{

    double T_final = 1.;
    size_t Nint = 50;
    size_t n = 16;
    double dt = (T_final - 0)/(double)Nint;
    double initial_data = 1.;
    ExplicitEuler<double,double> Euler;
    RK2<double, double> RK2_scheme;
    RK4<double,double> RK4_scheme;

    TimeIntegrator<double,double> solver(T_final,dt,initial_data,f);


    std::cout << "For " << Nint << " cells, " << std::endl; 
    std::cout << std::endl;
    std::cout << "Exact solution, exp("<<T_final<<") = " << std::setprecision(n) << std::exp(T_final) << std::endl;    

    double Euler_solution = solver.run(Euler);
    std::cout << "Explicit Euler, u(" << T_final <<") = "<< std::setprecision(n) << Euler_solution << std::endl;
    
    double RK2_solution = solver.run(RK2_scheme);
    std::cout << "RK2, u(" << T_final <<") = " << std::setprecision(n) << RK2_solution << std::endl;

    double RK4_solution = solver.run(RK4_scheme);
    std::cout << "RK4, u(" << T_final <<") = " << std::setprecision(n) << RK4_solution << std::endl;


    return 0;
}