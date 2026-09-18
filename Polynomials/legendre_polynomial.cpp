#include "legendre_polynomial.hpp"

LegendrePolynomial::LegendrePolynomial(size_t degree):polynomial(degree,0.)
{
    if(degree == 0)
    {
        (*this)[0] = 1;
    }
    else if(degree == 1)
    {
        (*this)[1] = 1;
    }
    else
    {
        double v[2] = {0.,1.};
        polynomial p_prev_prev(0,1);
        polynomial p_prev(1,v);
        polynomial p_current(2);

        for(auto ii=1; ii<degree; ii++)
        {
            p_current = (2*ii + 1.)/(ii + 1.) * p_prev.multiply_by_x() - ii/(ii + 1.)*p_prev_prev;
            p_prev_prev = p_prev;
            p_prev = p_current;
        }
        *this = p_current;
    }
}


LegendrePolynomial& LegendrePolynomial::operator=(polynomial const& p)
{
    this->polynomial::operator=(p);
    return *this;
}