#pragma once
#include "field.hpp"

template<typename T>
class numerical_flux
{
public:
    virtual T operator()(field<T> const& uh, size_t cell, physical_flux<T> const& flux) const = 0;

    virtual ~numerical_flux() = default;
};