#pragma once
#include "physical_flux.hpp"
#include "numerical_flux.hpp"
#include "field.hpp"

template<typename T>
class LF_flux : public numerical_flux<T>
{
public:
    T operator()(field<T> const& uh, size_t cell, physical_flux<T> const& flux) const;

};

template<typename T>
T LF_flux<T>::operator()(field<T> const& uh, size_t cell, physical_flux<T> const& flux) const
{
    T max_vel=0;
    for(auto ii=0; ii<=uh.n_cells(); ii++)
    {
        if(max_vel <= std::abs(flux.derivative(uh(ii))))
        {
            max_vel = std::abs(flux.derivative(uh(ii)));
        }
    }

    T ul, ur;

    // We use periodic boundary conditions for example
    if(cell == 0)
    {
        ul = uh(uh.n_cells()-1);
        ur = uh(0);
    }
    else if(cell == uh.n_cells())
    {
        ul = uh(uh.n_cells()-1);
        ur = uh(0);
    }
    else
    {
        ul = uh(cell-1);
        ur = uh(cell);
    }

    return 0.5 * ( flux(ul) + flux(ur) - max_vel * (ur - ul) );
}