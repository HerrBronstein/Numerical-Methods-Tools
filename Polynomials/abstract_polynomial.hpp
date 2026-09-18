#pragma once
#include <cmath>
#include <stdexcept>
#include <iostream>


class polynomial
{
    public:

    polynomial(size_t degree, double a = 0);
    polynomial(size_t degree, double tab[]);
    polynomial(polynomial const& p);
    ~polynomial(){delete[] coefficients_;};

    polynomial& operator=(polynomial const& p);

    double& operator[](size_t ii);
    double operator[](size_t ii) const;
    size_t degree() const {return degree_;};

    friend double calculate_polynomial(polynomial const& p, double x);
    double operator()(double x) const;

    friend polynomial operator+(polynomial const& p, polynomial const& q);
    friend polynomial operator-(polynomial const& p, polynomial const& q);
    friend polynomial operator-(polynomial const& p);
    friend polynomial operator*(double lambda, polynomial const& p);
    friend polynomial operator*(polynomial const& p, polynomial const& q);
    friend polynomial& operator+=(polynomial& p, polynomial const& q);
    friend std::ostream& operator<<(std::ostream& out, polynomial const& p);
    
    polynomial multiply_by_x() const;
    polynomial derivative() const;

    protected:

    double* coefficients_;
    size_t degree_;
};
