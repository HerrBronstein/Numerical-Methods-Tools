#pragma once
#include "finite_element.hpp"

template<typename P,typename T>
class dGSpace
{
    public:

    dGSpace(size_t order, const mesh<T>& msh);
    ~dGSpace(){};

    const size_t getOrder() const {return order_;};
    const FiniteElement<P,T>& getElement(size_t ii) const {return Elements_[ii];};
    FiniteElement<P,T>& getElement(size_t ii){return Elements_[ii];};

    size_t getGlobalDof(size_t cell, size_t localDof) const {return numLocalDof() * cell + localDof;};
    size_t getNumDof() const {return numLocalDof() * mesh_.n_cells();};
    size_t numLocalDof() const {return order_ + 1;};

    Eigen::SparseMatrix<double> assembleSIPMatrix(gauss_integrator& integ, double SIP_penalty) const;
    Eigen::SparseMatrix<double> assembleMassMatrix(gauss_integrator& integ) const;
    Eigen::VectorXd assembleLoadVector(std::function<double(double)> const& source_term, gauss_integrator& integ) const;

    void applyDirichletCondition(Eigen::VectorXd& B, double penalty, double x_left,double x_right) const;
    double evaluate_projection(std::vector<T> const& uh, size_t cell, double x) const;

    private:

    size_t order_;
    const mesh<T>& mesh_;
    std::vector<FiniteElement<P,T>> Elements_;
};

template<typename P, typename T>
dGSpace<P,T>::dGSpace(size_t order, const mesh<T>& msh):order_(order),mesh_(msh)
{
    FiniteElement<P,T> FE(order,msh,0);
    Elements_.push_back(FE);
    for(auto ii=1; ii<mesh_.n_cells(); ii++)
    {
        FE.changeCell(ii);
        Elements_.push_back(FE);
    }
}

template<typename P,typename T>
double dGSpace<P,T>::evaluate_projection(std::vector<T> const& uh,size_t cell, double x) const
{
    double res{0.};

    for(auto ii=0; ii<numLocalDof(); ii++)
    {
        res += uh[getGlobalDof(cell,ii)]*Elements_[cell].evaluate(ii,x);
    }

    return res;
}


template<typename P, typename T>
Eigen::SparseMatrix<double> dGSpace<P,T>::assembleSIPMatrix(gauss_integrator& integ,double SIP_penalty) const
{
    size_t dof = getNumDof();
    Eigen::SparseMatrix<double> A(dof,dof);

    Eigen::MatrixXd left_boundary = Elements_[0].local_SIP_boundary_left(SIP_penalty);
    Eigen::MatrixXd right_boundary = Elements_[mesh_.n_cells()-1].local_SIP_boundary_right(SIP_penalty);

    for(auto ii=0; ii<numLocalDof(); ii++)
    {
        for(auto jj=0; jj<numLocalDof(); jj++)
        {
            A.coeffRef(getGlobalDof(0,ii),getGlobalDof(0,jj)) += left_boundary(ii,jj);
            A.coeffRef(getGlobalDof(mesh_.n_cells()-1,ii), getGlobalDof(mesh_.n_cells()-1,jj)) += right_boundary(ii,jj);
        }
    }

    for(auto cell=0; cell<mesh_.n_cells(); cell++)
    {
        Eigen::MatrixXd local_volumic = Elements_[cell].local_SIP_volumic(integ);

        for(auto ii=0; ii< numLocalDof(); ii++)
        {
            size_t ii_global = getGlobalDof(cell,ii);
            for(auto jj=0; jj< numLocalDof() ; jj++)
            {
                size_t jj_global = getGlobalDof(cell,jj);
                A.coeffRef(ii_global,jj_global) += local_volumic(ii,jj);
            }
        }
    }

    for(auto cell=0; cell< mesh_.n_cells()-1; cell++)
    {

        Eigen::MatrixXd local_interface = Elements_[cell].local_SIP_faces_right();
        Eigen::MatrixXd local_penalty= Elements_[cell].local_SIP_penalty_right(SIP_penalty);

        for(auto ii=0; ii< 2 * numLocalDof(); ii++)
        {
            size_t I;

            if(ii < numLocalDof())
            {
                I = getGlobalDof(cell,ii);
            }
            else
            {
                I = getGlobalDof(cell+1, ii-numLocalDof());
            }

            for(auto jj=0; jj< 2*numLocalDof(); jj++)
            {
                size_t J;

                if(jj < numLocalDof())
                    J = getGlobalDof(cell, jj);
                else
                    J = getGlobalDof(cell+1, jj-numLocalDof());

                A.coeffRef(I,J) += local_interface(ii,jj) + local_penalty(ii,jj);            
            }
        }
    }

    // std::cout << "local volume\n"
    //       << Elements_[0].local_SIP_volumic(integ)
    //       << "\n";

    // std::cout << "right interface\n"
    //         << Elements_[1].local_SIP_faces_right()
    //         << "\n";

    // std::cout << "riht penalty\n"
    //         << Elements_[1].local_SIP_penalty_right(SIP_penalty)
    //         << "\n";

    return A;
}

template<typename P, typename T>
Eigen::SparseMatrix<double> dGSpace<P,T>::assembleMassMatrix(gauss_integrator& integ) const
{
    size_t dof = getNumDof();
    Eigen::SparseMatrix<double> A(dof,dof);

    for(auto cell=0; cell<mesh_.n_cells(); cell++)
    {
        Eigen::SparseMatrix<double> A_loc = Elements_[cell].local_mass(integ);
        for(auto ii=0; ii< numLocalDof(); ii++)
        {
            size_t ii_global = getGlobalDof(cell,ii);
            for(auto jj=0; jj< numLocalDof(); jj++)
            {
                size_t jj_global = getGlobalDof(cell,jj);
                A.coeffRef(ii_global,jj_global) += A_loc.coeffRef(ii,jj);
            }
        }
    }
    return A;
}

template<typename P, typename T>
Eigen::VectorXd dGSpace<P,T>::assembleLoadVector(std::function<double(double)> const& source_term, gauss_integrator& integ) const
{
    size_t dof = getNumDof();
    Eigen::VectorXd B(dof);

    for(auto cell=0; cell<mesh_.n_cells(); cell++)
    {
        Eigen::VectorXd B_loc = Elements_[cell].local_load(source_term, integ);
        for(auto ii=0; ii< numLocalDof(); ii++)
        {
            size_t ii_global = getGlobalDof(cell,ii);
            B(ii_global) = B_loc(ii);
        }
    }
    return B;
}

template<typename P, typename T>
void dGSpace<P,T>::applyDirichletCondition(Eigen::VectorXd& B, double penalty, double g_l,double g_r) const
{
    FiniteElement<P,T> Element_l = Elements_[0];
    FiniteElement<P,T> Element_r = Elements_[mesh_.n_cells()-1];

    double x_l = Element_l.getLeft();
    double x_r = Element_r.getRight();

    for(auto ii=0; ii< numLocalDof(); ii++)
    {
        B(ii) +=  g_l * Element_l.derivative(ii)(Element_l.mapToReference(x_l)) 
                    + penalty/mesh_.dx() * g_l * Element_l.evaluate(ii,x_l);

        B(getGlobalDof(mesh_.n_cells() - 1, ii)) += - g_r * Element_r.derivative(ii)(Element_l.mapToReference(x_r)) 
                                                    + penalty/mesh_.dx() * g_r * Element_r.evaluate(ii,x_r) ;
    }    
}