#pragma once
#include <functional>


template<typename T, typename V>
class AbstractScheme
{
    public:
    
    virtual ~AbstractScheme() = default;

    virtual V step(T t,T dt,V const& u,std::function<V(T, V const&)> const& func) const = 0;
};