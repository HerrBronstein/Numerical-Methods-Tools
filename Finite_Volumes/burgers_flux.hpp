#pragma once
#include "physical_flux.hpp"

template<typename T>
class burgers_flux : public physical_flux<T>
{
public:
    T operator()(T const& u) const override {return 0.5 * u * u;};

    T derivative(T const& u) const override {return u;};
};