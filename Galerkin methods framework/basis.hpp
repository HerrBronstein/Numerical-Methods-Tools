#pragma once
#include "lagrange_polynomial.hpp"

class basis
{
    public:

    basis(int order);
    ~basis(){};

    LagrangePolynomial& operator[](int ii){return basis_[ii];};
    const LagrangePolynomial& operator[](int ii) const {return basis_[ii];};
    const std::vector<LagrangePolynomial>& getPoly() const {return basis_;};
    std::vector<LagrangePolynomial>& getPoly() {return basis_;};
    int getSize(){return size_;};
    double getNode(size_t ii){return reference_nodes_[ii];};

    private:

    std::vector<LagrangePolynomial> basis_;
    std::vector<double> reference_nodes_;
    int size_;
};
