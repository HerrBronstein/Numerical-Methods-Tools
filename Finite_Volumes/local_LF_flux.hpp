#pragma once
#include "physical_flux.hpp"
#include "numerical_flux.hpp"
#include "field.hpp"

template<typename T>
class local_LF_flux : public numerical_flux<T>
{
public:
    T operator()(field<T> const& uh, size_t cell, physical_flux<T> const& flux) const;

};

template<typename T>
T local_LF_flux<T>::operator()(field<T> const& uh, size_t cell, physical_flux<T> const& flux) const
{
    T ul, ur, gamma;

    // We use ghost cells with the neighboring data as boundary condition
    if(cell == 0)
    {
        ul = uh(0);
        ur = uh(0);
    }
    else if(cell == uh.n_cells()-1)
    {
        ul = uh(uh.n_cells()-1);
        ur = uh(uh.n_cells()-1);
    }
    else
    {
        ul = uh(cell-1);
        ur = uh(cell);
    }
    gamma = std::max(std::abs(ul),std::abs(ur));

    return 0.5 * ( flux(ul) + flux(ur) - gamma * (ur - ul) );
}