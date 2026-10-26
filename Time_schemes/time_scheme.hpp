#pragma once
#include <functional>
#include "abstract_scheme.hpp"
#include <algorithm>

template<typename T,typename V>
class TimeIntegrator
{
    public:
    TimeIntegrator(T final, T dt, V initial_data, std::function<V(T,V const&)> const& f):
    final_time_(final),dt_(dt),initial_data_(initial_data),f_(f){};

    ~TimeIntegrator(){};

    T final_time() const {return final_time_;};
    T const& dt() const {return dt_;};

    V const& initial_data() const {return initial_data_;};
    V& initial_data() {return initial_data_;};

    std::function<V(T,V const&)> get_f() const {return f_;};

    void change_dt(T dt){dt_ = dt;};
    V run(AbstractScheme<T,V> const& scheme);

    protected:

    T final_time_;
    T dt_;
    T physical_time_{0};
    V initial_data_;

    std::function<V(T,V const&)> f_;
};


template<typename T, typename V>
V TimeIntegrator<T,V>::run(AbstractScheme<T,V> const& scheme)
{
    V solution = initial_data_;
    physical_time_ = T(0);
    while(physical_time_ < final_time_)
    {
        T dt_current = std::min(dt_, final_time_ - physical_time_);

        solution = scheme.step(physical_time_,dt_current,solution,f_);

        physical_time_ += dt_current;
    }
    return solution;
}