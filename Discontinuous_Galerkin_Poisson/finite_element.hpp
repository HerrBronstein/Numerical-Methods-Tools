#pragma once
#include "abstract_polynomial.hpp"
#include "reference_element.hpp"
#include "mesh_1d.hpp"
#include "gauss_integrator.hpp"
#include <Eigen/Sparse>

template<typename P,typename T>
class FiniteElement
{
    public:

    FiniteElement(size_t order,mesh<T> const& msh,size_t cell);
    ~FiniteElement(){};

    double mapToReference(double x) const;
    double mapToReal(double xi) const;
    double evaluate(size_t ii, double x) const;

    polynomial derivative(size_t ii) const;

    T getLeft() const {return x_left_;};
    T getRight() const {return x_right_;};
    const P& getBasis(size_t ii) const {return reference_.getBasis(ii);};
    const size_t& getOrder() const {return order_;};

    const size_t Dof() const {return order_ + 1;};

    void changeCell(size_t newCell);

    Eigen::SparseMatrix<double> local_mass(gauss_integrator& integ) const;
    Eigen::VectorXd local_load(std::function<double(double)> const& source_term, gauss_integrator& integ) const;
    
    Eigen::MatrixXd local_SIP_volumic(gauss_integrator& integ) const;

    Eigen::MatrixXd local_SIP_penalty_right(double eta) const;
    Eigen::MatrixXd local_SIP_penalty_left(double eta) const;

    Eigen::MatrixXd local_SIP_faces_left() const;
    Eigen::MatrixXd local_SIP_faces_right() const;

    Eigen::MatrixXd local_SIP_boundary_left(double eta) const;
    Eigen::MatrixXd local_SIP_boundary_right(double eta) const;

    private:

    const referenceElement<P> reference_; 
    size_t order_;
    const mesh<T>& mesh_;

    size_t cell_;
    T x_left_;
    T x_right_;
};

template<typename P,typename T>
FiniteElement<P,T>::FiniteElement(size_t order, mesh<T> const& msh, size_t cell):reference_(order),order_(order),mesh_(msh),cell_(cell)
{
    x_left_ = msh.xc(cell) - 0.5*msh.dx();
    x_right_ = msh.xc(cell) + 0.5*msh.dx();
}

template<typename P,typename T>
void FiniteElement<P,T>::changeCell(size_t newCell)
{
    cell_ = newCell;
    x_left_ = mesh_.xc(newCell) - 0.5*mesh_.dx();
    x_right_ = mesh_.xc(newCell) + 0.5*mesh_.dx();
}

template<typename P, typename T>
double FiniteElement<P,T>::mapToReference(double x) const
{
    return 2*(x - x_left_)/mesh_.dx() - 1;
}

template<typename P, typename T>
double FiniteElement<P,T>::mapToReal(double xi) const
{
    return mesh_.xc(cell_) + 0.5 * mesh_.dx() * xi ;
}

template<typename P, typename T>
double FiniteElement<P,T>::evaluate(size_t ii, double x) const
{
    return reference_.getBasis(ii)(mapToReference(x));
}

template<typename P, typename T>
polynomial FiniteElement<P,T>::derivative(size_t ii) const
{
    return 2/mesh_.dx() * reference_.getBasis(ii).derivative();
}


template<typename P,typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_volumic(gauss_integrator& integ) const
{
    Eigen::MatrixXd A(Dof(),Dof());
    double J = 0.5 * mesh_.dx();

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            integ.change_integrand(J * derivative(ii) * derivative(jj));
            A(ii,jj) = integ.compute_integral(1);
        }
    }

    return A;
}

template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_penalty_right(double eta) const
{
    size_t n = 2*Dof();

    Eigen::MatrixXd penalty = Eigen::MatrixXd::Zero(n,n);
    double K = eta / mesh_.dx();

    if(cell_ == mesh_.n_cells() - 1)
    {
        return penalty;
    }

    FiniteElement<P,T> right_element(getOrder(), mesh_, cell_ + 1);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            penalty(ii,jj) = K* evaluate(ii,getRight()) * evaluate(jj,getRight()); // - -

            penalty(Dof() + ii, jj) = - K * right_element.evaluate(ii,right_element.getLeft()) * evaluate(jj,getRight()); // - +

            penalty(ii,Dof() + jj) = - K * evaluate(ii,getRight()) * right_element.evaluate(jj,right_element.getLeft()); // + -

            penalty(Dof() + ii, Dof() + jj) = K * right_element.evaluate(ii,right_element.getLeft()) * right_element.evaluate(jj,right_element.getLeft()); // + +
        }
    }


    return penalty;
}

template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_penalty_left(double eta) const
{
    size_t n = 2*Dof();

    Eigen::MatrixXd penalty = Eigen::MatrixXd::Zero(n,n);
    double K = eta / mesh_.dx();

    if(cell_ == 0)
    {
        return penalty;
    }

    FiniteElement<P,T> left_element(getOrder(), mesh_, cell_ -1);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            penalty(ii,jj) = K* left_element.evaluate(ii,left_element.getRight()) * left_element.evaluate(jj,left_element.getRight()); // - -

            penalty(Dof() + ii, jj) = - K * evaluate(ii,getLeft()) * left_element.evaluate(jj,left_element.getRight()); // + -

            penalty(ii,Dof() + jj) = - K * left_element.evaluate(ii,left_element.getRight()) * evaluate(jj,getLeft()); // + -

            penalty(Dof() + ii, Dof() + jj) = K * evaluate(ii,getLeft()) * evaluate(jj,getLeft()); // + +
        }
    }

    return penalty;
}


template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_faces_right() const
{
    size_t n = 2*Dof();
    Eigen::MatrixXd faces(n,n);

    FiniteElement<P,T> right_element(getOrder(),mesh_, cell_+1);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            faces(ii,jj) = - 0.5 * ( derivative(ii)(mapToReference(getRight())) * evaluate(jj,getRight()) 
                                    + evaluate(ii,getRight()) * derivative(jj)(mapToReference(getRight())) );

            faces(Dof() + ii, jj) = - 0.5 * ( right_element.derivative(ii)(right_element.mapToReference(right_element.getLeft())) * evaluate(jj,getRight()) 
                                            - right_element.evaluate(ii,right_element.getLeft()) * derivative(jj)(mapToReference(getRight())) );

            faces(ii, Dof() + jj) = - 0.5 * ( - derivative(ii)(mapToReference(getRight())) * right_element.evaluate(jj,right_element.getLeft())
                                            + evaluate(ii,getRight()) * right_element.derivative(jj)(right_element.mapToReference(right_element.getLeft())) );

            faces(Dof() + ii, Dof() + jj) = - 0.5 * ( - right_element.derivative(ii)(right_element.mapToReference(right_element.getLeft())) * right_element.evaluate(jj,right_element.getLeft())
                                                    - right_element.evaluate(ii,right_element.getLeft()) * right_element.derivative(jj)(right_element.mapToReference(right_element.getLeft())));
        }
    }
    return faces;
}

template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_faces_left() const
{
    size_t n = 2*Dof();
    Eigen::MatrixXd faces(n,n);

    FiniteElement<P,T> left_element(getOrder(),mesh_, cell_-1);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            faces(ii,jj) = - 0.5 * ( left_element.derivative(ii)(left_element.mapToReference(left_element.getRight())) * left_element.evaluate(jj,left_element.getRight()) 
                                    + left_element.evaluate(ii,left_element.getRight()) * left_element.derivative(jj)(left_element.mapToReference(left_element.getRight())) );

            faces(Dof() + ii, jj) = - 0.5 * ( derivative(ii)(mapToReference(getLeft())) * left_element.evaluate(jj,left_element.getRight()) 
                                            - evaluate(ii,getLeft()) * left_element.derivative(jj)(left_element.mapToReference(left_element.getRight())) );

            faces(ii, Dof() + jj) = - 0.5 * ( - left_element.derivative(ii)(left_element.mapToReference(left_element.getRight())) * evaluate(jj,getLeft()) 
                                            + left_element.evaluate(ii,left_element.getRight()) * derivative(jj)(mapToReference(getLeft())) );

            faces(Dof() + ii, Dof() + jj) = - 0.5 * ( - derivative(ii)(mapToReference(getLeft())) * evaluate(jj,getLeft())
                                                    - evaluate(ii,getLeft()) * derivative(jj)(mapToReference(getLeft())));
        }
    }
    return faces;
}

template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_boundary_right(double eta) const
{
    Eigen::MatrixXd boundary(Dof(),Dof());
    T x_r = mesh_.vertice_coord(mesh_.vertices() - 1);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            boundary(ii,jj) = - ( derivative(ii)(mapToReference(x_r)) * evaluate(jj,x_r) 
                                + evaluate(ii,x_r) * derivative(jj)(mapToReference(x_r)) )
                                + eta/mesh_.dx() * evaluate(ii,x_r)*evaluate(jj,x_r);
        }
    }
    return boundary;
}

template<typename P, typename T>
Eigen::MatrixXd FiniteElement<P,T>::local_SIP_boundary_left(double eta) const
{
    Eigen::MatrixXd boundary(Dof(),Dof());
    T x_l = mesh_.vertice_coord(0);

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj<Dof(); jj++)
        {
            boundary(ii,jj) = - ( - derivative(ii)(mapToReference(x_l)) * evaluate(jj,x_l) 
                                - evaluate(ii,x_l) * derivative(jj)(mapToReference(x_l)) )
                                + eta/mesh_.dx() * evaluate(ii,x_l) * evaluate(jj,x_l);
        }
    }
    return boundary;
}

template<typename P, typename T>
Eigen::SparseMatrix<double> FiniteElement<P,T>::local_mass(gauss_integrator& integ) const
{
    Eigen::SparseMatrix<double> A(Dof(),Dof());
    double J = 0.5 * mesh_.dx();

    for(auto ii=0; ii<Dof(); ii++)
    {
        for(auto jj=0; jj< Dof(); jj++)
        {
            integ.change_integrand(J * reference_.getBasis(ii) * reference_.getBasis(jj));
            A.coeffRef(ii,jj) = integ.compute_integral(1);
        }
    }
    return A;
}

template<typename P, typename T>
Eigen::VectorXd FiniteElement<P,T>::local_load(std::function<double(double)> const& source_term, gauss_integrator& integ) const
{
    Eigen::VectorXd B(Dof());
    double J = 0.5 * mesh_.dx();

    for(auto ii=0; ii< Dof(); ii++)
    {
        const P& phi = reference_.getBasis(ii);
        std::function<double(double)> load_entry = [source_term, &phi, J, this](double x)
        {
            return J * source_term(mapToReal(x)) * phi(x);
        };
        integ.change_integrand(load_entry);
        B(ii) = integ.compute_integral(1);
    }
    return B;
}