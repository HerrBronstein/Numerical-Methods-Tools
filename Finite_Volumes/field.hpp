#pragma once
#include "functional"
#include "mesh_1d.hpp"

template<typename T>
class residual;

template<typename T>
class field
{
    public:

    field(std::function<T(T)> initial_value, mesh<T> const& msh);

    inline const size_t n_cells() const {return n_cells_;};
    inline const std::vector<T>& operator()() {return values_;};
    inline double operator()(size_t ii) const {return values_[ii];};

    field<T>& operator+=(residual<T> res);


    protected:
    std::vector<T> values_;
    size_t n_cells_{0};
};

template<typename T>
field<T>::field(std::function<T(T)> initial_value, mesh<T> const& msh):n_cells_(msh.n_cells())
{
    values_.resize(n_cells());
    for(auto ii=0; ii<n_cells(); ii++)
    {
        values_[ii] = initial_value(msh.xc(ii));
    }
}

template<typename T>
field<T>& field<T>::operator+=(residual<T> res)
{
    for(auto ii=0; ii<n_cells(); ii++)
    {
        values_[ii] += res(ii);
    }
    return *this;
}