#pragma once

template<typename T>
class physical_flux
{
public:
    virtual T operator()(T const& u) const = 0;

    virtual T derivative(T const& u) const = 0;

    virtual ~physical_flux() = default;
};