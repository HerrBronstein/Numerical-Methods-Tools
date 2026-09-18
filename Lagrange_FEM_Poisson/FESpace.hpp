#pragma once
#include "finite_element.hpp"

template<typename T>
class FESpace
{
    public:

    FESpace(mesh<T> const& msh, int order);
    ~FESpace(){};

    const FiniteElement<T>& getFE(size_t ii){return FiniteElements_[ii];};
    const mesh<T>& getMesh() const {return mesh_;};

    LagrangePolynomial& getPolynomial(size_t cell, size_t num_poly){return FiniteElements_[cell].getPoly()[num_poly];};

    size_t numElement() const {return FiniteElements_.size(); };
    size_t numDof() const;
    size_t numLocalDof() const {return order_ + 1;};
    size_t getGlobalDof(size_t ElementIndex, size_t localDof) const;
    const double getNode(size_t ii) const {return nodes_[ii];};
    
    Eigen::MatrixXd assembleStiffnessMatrix(gauss_integrator& integ);
    Eigen::VectorXd assembleLoadVector(const std::function<T(T)>& source_term, gauss_integrator& integ);

    void applyDirichletCondition(Eigen::MatrixXd& A, Eigen::VectorXd& B, double left_val, double right_val) const;

    double evaluate_projection(std::vector<T> const& uh, double x);

    private:

    const mesh<T>& mesh_;
    std::vector<FiniteElement<T>> FiniteElements_;
    std::vector<double> nodes_;
    int order_;
};


template<typename T>
FESpace<T>::FESpace(mesh<T> const& msh, int order):mesh_(msh),order_(order)
{
    FiniteElements_.reserve(msh.n_cells());
    for(auto ii=0; ii< msh.n_cells(); ii++)
    {
        FiniteElements_.emplace_back(msh,ii,order);
        
        for(auto jj=0; jj< numLocalDof()-1; jj++)
        {
            nodes_.emplace_back(FiniteElements_[ii].getNode(jj));
        }
    }
    nodes_.emplace_back(FiniteElements_[msh.n_cells()-1].getNode(numLocalDof()-1));
}

template<typename T>
double FESpace<T>::evaluate_projection(std::vector<T> const& uh, double x)
{
    double res{0.};
    size_t cell = mesh_.getMesh(x);

    for(auto ii=0; ii<numLocalDof(); ii++)
    {
        res += uh[getGlobalDof(cell,ii)]*FiniteElements_[cell].evaluate(ii,x);
    }

    return res;
}

template<typename T>
size_t FESpace<T>::getGlobalDof(size_t ElementIndex, size_t localDof) const
{
    return ElementIndex * order_ + localDof;
}

template<typename T>
size_t FESpace<T>::numDof() const
{
    if(order_ == 0)
    {
        return mesh_.n_cells();
    }
    else
    {
        return order_ * mesh_.n_cells() + 1;
    }
}

template<typename T>
Eigen::MatrixXd FESpace<T>::assembleStiffnessMatrix(gauss_integrator& integ)
{
    size_t Dof = numDof();
    Eigen::MatrixXd A(Dof,Dof);

    for(auto cell=0; cell < numElement(); cell++)
    {
        Eigen::MatrixXd A_loc = FiniteElements_[cell].local_stiffness(integ);

        for(auto ii=0; ii< numLocalDof() ; ii++)
        {
            size_t ii_global = getGlobalDof(cell,ii);

            for(auto jj=0; jj< numLocalDof(); jj++)
            {
                size_t jj_global = getGlobalDof(cell,jj);

                A(ii_global,jj_global) += A_loc(ii,jj);
            }
        }
    }
    return A;
}

template<typename T>
Eigen::VectorXd FESpace<T>::assembleLoadVector(const std::function<T(T)>& source_term, gauss_integrator& integ)
{
    size_t Dof = numDof();
    Eigen::VectorXd B = Eigen::VectorXd::Zero(Dof);

    for(auto cell=0; cell< numElement(); cell++)
    {
        Eigen::VectorXd B_loc = FiniteElements_[cell].local_load(source_term,integ);
        for(auto ii=0; ii< numLocalDof(); ii++)
        {
            size_t ii_global = getGlobalDof(cell,ii);
            B(ii_global) += B_loc(ii);
        }
    }

    return B;
}

template<typename T>
void FESpace<T>::applyDirichletCondition(Eigen::MatrixXd& A, Eigen::VectorXd& B, double left_val, double right_val) const
{
    for(auto ii=0; ii<numDof(); ii++)
    {
        A(0,ii) = 0.;
        A(ii,0) = 0.;
        A(numDof() - 1, ii) = 0.;
        A(ii,numDof() -1) = 0.;
    }
    A(0,0) = 1.;
    A(numDof()-1,numDof()-1) = 1.;

    B(0) = left_val;
    B(numDof()-1) = right_val;
}
