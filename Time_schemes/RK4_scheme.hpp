#pragma once
#include "abstract_scheme.hpp"

template<typename T, typename V>
class RK4 : public AbstractScheme<T,V>
{
    public:

    ~RK4(){};


    V step(T t,T dt,V const& u,std::function<V(T, V const&)> const& func) const;

};

template<typename T, typename V>
V RK4<T,V>::step(T t, T dt, V const& u, std::function<V(T,V const&)> const& func) const
{
    V k1 = func(t,u);
    V k2 = func(t + 0.5*dt ,u + 0.5*dt*k1);
    V k3 = func(t + 0.5*dt, u + 0.5*dt*k2);
    V k4 = func(t + dt, u + dt*k3);

    return u + dt/6. * (k1 + 2*k2 + 2*k3 + k4);
}