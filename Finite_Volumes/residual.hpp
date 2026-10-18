#pragma once 
#include "field.hpp"
#include "physical_flux.hpp"
#include "numerical_flux.hpp"

// template<typename T>
// T flux_burgers(T const& u)
// {
//     return 0.5 * u * u;
// }

// template<typename T>
// T flux_LF(T const& ul, T const& ur, T max_vel, physical_flux<T> const& flux)
// {
//     return 0.5 * (flux(ul) + flux(ur) - max_vel * (ur - ul) );
// }

template<typename T>
class residual
{
    public:

    residual(mesh<T> const& msh):mesh_(msh),values_(msh.n_cells()){};

    inline void clear() {values_.clear();};
    void assemble_from_field(field<T>& uh, physical_flux<T> const& flux, numerical_flux<T> const& num_flux);


    const mesh<T>& getMesh() const {return mesh_;};
    inline const size_t n_cells() const {return mesh_.n_cells();}; 

    T operator()(size_t ii) const {return values_[ii];};
    T& operator()(size_t ii) {return values_[ii];};

    template<typename S>
    friend residual<S> operator*(S lambda, residual<S> const& res);

    protected:

    std::vector<T> values_;
    mesh<T> const& mesh_;

};

template<typename S>
residual<S> operator*(S lambda, residual<S> const& res)
{
    residual<S> prod(res.getMesh());
    for(auto ii=0; ii<res.n_cells(); ii++)
    {
        prod(ii) = lambda * res(ii);
    }

    return prod;
} 

template<typename T>
void residual<T>::assemble_from_field(field<T>& uh, physical_flux<T> const& flux, numerical_flux<T> const& num_flux)
{
    for(auto ii=0; ii<mesh_.n_cells(); ii++)
    {
        values_[ii] = - ( num_flux(uh,ii+1,flux) - num_flux(uh,ii,flux) )/ mesh_.dx() ;
    }

}