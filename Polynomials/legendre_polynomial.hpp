#pragma once
#include "abstract_polynomial.hpp"


class LegendrePolynomial : public polynomial
{
    public:

    LegendrePolynomial(size_t degree);
    ~LegendrePolynomial(){};

    LegendrePolynomial& operator=(polynomial const& p);

};