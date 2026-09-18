#pragma once
#include "abstract_polynomial.hpp"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <Eigen/Dense>

class LagrangePolynomial : public polynomial
{
    public:

    LagrangePolynomial(const std::vector<double>& x_values, const std::vector<double>& y_values);
    ~LagrangePolynomial(){};

    std::vector<double> compute_Lagrange_coef(const std::vector<double>& x_values, const std::vector<double>& y_values);

    void addPoint(double x, double y);
    double x_val(int ii){return x_val_[ii];};
    double y_val(int ii){return y_val_[ii];};

    protected:

    std::vector<double> x_val_;
    std::vector<double> y_val_;
};