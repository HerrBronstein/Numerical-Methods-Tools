#pragma once
#include "basis.hpp"
#include "mesh_1d.hpp"
#include "gauss_integrator.hpp"

template<typename T>
class FiniteElement
{
    public: 

    FiniteElement(const mesh<T>& msh, size_t cell, int order);
    ~FiniteElement(){};

    double mapToReference(double x) const;
    double mapToReal(double xi) const;
    double evaluate(int ii, double x);

    const std::vector<LagrangePolynomial>& getPoly() const {return basis_.getPoly();};
    std::vector<LagrangePolynomial>& getPoly(){return basis_.getPoly();};
    const basis& getBasis(){return basis_;};
    T getLeft() const {return x_left_;};
    T getRight() const {return x_right_;};
    T dx() const {return dx_;};
    double getNode(size_t ii){return nodes_[ii];};
    size_t getDof(){return nodes_.size();};

    polynomial derivative(int ii) const;

    Eigen::MatrixXd local_stiffness(gauss_integrator& integ) const;
    Eigen::MatrixXd local_mass(gauss_integrator& integ) const;
    Eigen::VectorXd local_load(std::function<double(double)> const& source_term, gauss_integrator& integ) const;

    private:

    basis basis_;
    size_t cell_;
    std::vector<double> nodes_;
    T x_left_;
    T x_right_;
    T dx_;
};


template<typename T>
FiniteElement<T>::FiniteElement(const mesh<T>& msh, size_t cell, int order):basis_(order),cell_(cell),dx_(msh.dx())
{
    x_left_ = msh.xc(cell_) - 0.5*msh.dx();
    x_right_ = msh.xc(cell_) + 0.5*msh.dx();
    for(auto ii=0; ii<=order; ii++)
    {
        nodes_.push_back(mapToReal((basis_.getNode(ii))));
    }
}

template<typename T>
double FiniteElement<T>::mapToReference(double x) const
{
    return 2 * (x - x_left_)/dx_ - 1;
}

template<typename T>
double FiniteElement<T>::mapToReal(double xi) const
{
    return 0.5*(x_left_ + x_right_ + xi*dx_);
}

template<typename T>
double FiniteElement<T>::evaluate(int ii, double x)
{
    return getPoly()[ii](mapToReference(x));
}

template<typename T>
polynomial FiniteElement<T>::derivative(int ii) const
{
    double a = 2./dx_;
    return a*basis_[ii].derivative();
}

template<typename T>
Eigen::MatrixXd FiniteElement<T>::local_stiffness(gauss_integrator& integ) const
{
    size_t dof = nodes_.size();
    double J = 0.5*dx_;
    Eigen::MatrixXd A(dof,dof);

    for(auto ii=0; ii<dof; ii++)
    {
        for(auto jj=0; jj<dof; jj++)
        {
            integ.change_integrand(J * derivative(ii) * derivative(jj));
            A(ii,jj) = integ.compute_integral(1);
        }
    }
    return A;
}

template<typename T>
Eigen::MatrixXd FiniteElement<T>::local_mass(gauss_integrator& integ) const
{
    size_t dof = nodes_.size();
    double J = 0.5*dx_;
    Eigen::MatrixXd M(dof,dof);

    for(auto ii=0; ii<dof; ii++)
    {
        for(auto jj=0; jj<dof; jj++)
        {
            integ.change_integrand(J * getPoly()[ii] * getPoly()[jj]);
            M(ii,jj) = integ.compute_integral(1);
        }
    }
    return M;   
}

template<typename T>
Eigen::VectorXd FiniteElement<T>::local_load(std::function<double(double)> const& source_term, gauss_integrator& integ) const
{
    size_t dof = nodes_.size();
    double J = 0.5*dx_;
    Eigen::VectorXd B(dof);

    for(auto ii=0; ii<dof; ii++)
    {
        const LagrangePolynomial& phi = getPoly()[ii];
        std::function<double(double)> load_entry = [source_term, &phi, J, this](double x)
        {
            return J * source_term(mapToReal(x)) * phi(x);
        };
        integ.change_integrand(load_entry);
        B(ii) = integ.compute_integral(1);
    }
    return B;
}