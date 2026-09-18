#include "lagrange_polynomial.hpp"

 
std::vector<double> LagrangePolynomial::compute_Lagrange_coef(const std::vector<double>& x_values, const std::vector<double>& y_values)
{

    if(x_values.size() != y_values.size())
    {
        throw std::invalid_argument("x_values and y_values should be of the same size");
    }
    else
    {
        int n = x_values.size();
        Eigen::MatrixXd V(n,n);
        for(auto ii=0; ii<n; ii++)
        {
            for(auto jj=0; jj<n; jj++)
            {
                V(ii,jj) = std::pow(x_values[ii],jj);
            }
        }
        Eigen::VectorXd Y(n);
        for(auto ii=0; ii<n; ii++)
        {
            Y(ii) = y_values[ii];
        }

        Eigen::VectorXd A = V.colPivHouseholderQr().solve(Y);

        std::vector<double> A_vector(n);
        for(auto ii=0; ii<n; ii++)
        {
            A_vector[ii] = A(ii);
        }
        return A_vector;
    }
}

LagrangePolynomial::LagrangePolynomial(const std::vector<double>& x_values, const std::vector<double>& y_values):polynomial(x_values.size()-1,compute_Lagrange_coef(x_values,y_values).data()),x_val_(x_values),y_val_(y_values){}


void LagrangePolynomial::addPoint(double x, double y)
{
    x_val_.push_back(x);
    y_val_.push_back(y);

    LagrangePolynomial p(x_val_,y_val_);

    *this = p;
}

// double LagrangePolynomial::evaluate(double x) const
// {
//     double res{0.};
//     int n = x_values_.size();
//     for(int ii = 0; ii<n; ii++)
//     {
//         double prod{1.};

//         for(int jj=0; jj<ii; jj++)
//         {
//             prod *= (x - x_values_[jj])/(x_values_[ii] - x_values_[jj]);
//         }
//         for(int jj=ii+1; jj<n; jj++)
//         {
//             prod *= (x - x_values_[jj])/(x_values_[ii] - x_values_[jj]);
//         }

//         res += y_values_[ii]*prod;
//     }
//     return res;
// }

// double LagrangePolynomial::operator()(double x) const
// {
//     return evaluate(x);
// }
