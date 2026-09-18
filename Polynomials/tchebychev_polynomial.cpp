#include "tchebychev_polynomial.hpp"

tchebychevPolynomial::tchebychevPolynomial(size_t degree):polynomial(degree,0.)
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
        polynomial T_prev_prev(0,1.);
        double v[2] = {0.,1.};
        polynomial T_prev(1,v);
        polynomial T_current(2);

        for(auto ii=1; ii<degree; ii++)
        {
            T_current = 2*T_prev.multiply_by_x() - T_prev_prev;
            T_prev_prev = T_prev;
            T_prev = T_current;
        }
        *this = T_current;
    }
}

tchebychevPolynomial& tchebychevPolynomial::operator=(polynomial const& p)
{
    this->polynomial::operator=(p);
    return *this;
}