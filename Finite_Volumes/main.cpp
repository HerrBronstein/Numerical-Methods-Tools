#include <iostream>
#include <cmath>
#include "time_loop.hpp"
#include "output_writer.hpp"
#include "burgers_flux.hpp"
#include "LF_flux.hpp"
#include "local_LF_flux.hpp"
#include "BBM_flux.hpp"

//---------------------------------------------------------------------------------------
//--------------------------------- Initial Data Examples--------------------------------

double gauss_pulse(double x)
{
    return std::exp(- std::pow(x-5.,2));
}

double sinusoid(double x)
{
    return std::sin(x);
}

double soliton(double x)
{
    double c = 1.2; //Celerity
    double a = 3 * (c-1.); //Amplitude
    double delta = 7; //Phase shift

    double X = 0.5 * std::sqrt((c - 1.)/c) * (x - delta) ; //Argument

    double sech = 1/std::cosh(X);
    return  a * sech * sech ;
}

int main()
{
    // ----------------------------------------------------------------------------------
    //----------------------------------- Parameters ------------------------------------

    double x_l = - 3.14;
    double x_r = 3.14;
    size_t n_cells = 1000;
    double T_final = 4;


    //-----------------------------------------------------------------------------------
    //-------------------------- Definitions & initializations --------------------------

    mesh<double> msh(x_l,x_r,n_cells); //Definition of the mesh
    output_writer<double> output(msh,"FV_gauss_pulse"); //Initialize output
    field<double> uh(sinusoid,msh); //Initialize data

    output.write_solution(uh,"initial"); //First output

    burgers_flux<double> flux; //Initializing physical flux
    BBM_flux<double> flux_BBM;

    LF_flux<double> num_flux; //Initializing numerical flux
    // local_LF_flux<double> num_local_flux;
 

    //-----------------------------------------------------------------------------------
    //-------------------------------- Time solving ------------------------------------- 
    
    time_loop<double> time(msh,flux_BBM,num_flux,T_final);
    time.run(uh); //Time loop


    //----------------------------------- Output ----------------------------------------
    output.write_solution(uh,"final"); //Final output

    return 0;
}   
