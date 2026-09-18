#include "basis.hpp"

basis::basis(int order)
{
    if (order < 0) 
    {
    throw std::invalid_argument("order must be a positive integer");
    }
    size_ = order + 1;
    std::vector<double> x_val(size_);
    std::vector<double> y_val(size_);

    if(order == 0)
    {
        x_val[0] = 0.;
        y_val[0] = 1.;
        LagrangePolynomial p(x_val,y_val);
        reference_nodes_.push_back(x_val[0]);
        basis_.push_back(p);
    }
    else
    {
        double h = 2.0/order; //Reference element is [-1,1]
        for(auto ii = 0; ii<size_; ii++)
        {
            x_val[ii] = -1 + h*ii;
            reference_nodes_.push_back(x_val[ii]);
        }
        for(auto ii=0; ii<size_; ii++)
        {
            y_val[ii] = 1.;
            basis_.emplace_back(x_val,y_val);
            y_val[ii] = 0.; 
        }
    }
}
