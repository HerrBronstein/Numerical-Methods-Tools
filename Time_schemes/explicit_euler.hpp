#pragma once
#include "abstract_scheme.hpp"

template<typename T, typename V>
class ExplicitEuler : public AbstractScheme<T,V>
{
    public:

    ~ExplicitEuler(){};


    V step(T t,T dt,V const& u,std::function<V(T, V const&)> const& func) const;

};

template<typename T, typename V>
V ExplicitEuler<T,V>::step(T t, T dt, V const& u, std::function<V(T,V const&)> const& func) const
{
    return u + dt * func(t,u);
}