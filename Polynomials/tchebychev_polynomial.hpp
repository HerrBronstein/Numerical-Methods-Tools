#pragma once
#include "abstract_polynomial.hpp"

class tchebychevPolynomial : public polynomial
{
    public:

    tchebychevPolynomial(size_t degree);
    ~tchebychevPolynomial(){};

    tchebychevPolynomial& operator=(polynomial const& p);

};
