#pragma once
#include "abstract_scheme.hpp"

template<typename T, typename V>
class RK2 : public AbstractScheme<T,V>
{
    public:

    ~RK2(){};


    V step(T t,T dt,V const& u,std::function<V(T, V const&)> const& func) const;

};

template<typename T, typename V>
V RK2<T,V>::step(T t, T dt, V const& u, std::function<V(T,V const&)> const& func) const
{
    V u_predicted = u + 0.5 * dt * func(t,u);
    return u + dt * func(t + dt * 0.5 ,u_predicted);
}