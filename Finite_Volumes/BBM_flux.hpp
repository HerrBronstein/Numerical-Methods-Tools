#pragma once
#include "physical_flux.hpp"

template<typename T>
class BBM_flux : public physical_flux<T>
{
public:
    T operator()(T const& u) const override {return u + 0.5 * u * u;};

    T derivative(T const& u) const override {return 1 + u;};
};