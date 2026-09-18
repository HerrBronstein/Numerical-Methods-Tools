#include <iostream>
#include <cmath>
#include <cstdlib>
#include "dGSpace.hpp"
// #include "finite_element.hpp"
#include "legendre_polynomial.hpp"
#include "tchebychev_polynomial.hpp"
#include "mesh_1d.hpp"
#include "gauss_integrator.hpp"
#include "solver.hpp"
#include "output_writer.hpp"

double f(double x)
{
    return 1;
}

double gauss_bell(double x)
{
    return -std::exp(- x * x);
}

double log(double x)
{
    return std::log(std::abs(x));
}

double absolute(double x)
{
    if(x>= 0)
    {
        return x;
    }
    else
    {
        return -x;
    }
}

double polynome(double x)
{
    return x*x;
}


int main()
{
    //------------------------------Test template base polynômiale-------------------------------
    // size_t order = 1;

    // referenceElement<LegendrePolynomial> Ref_legendre(order);

    // std::cout << "Legendre basis of order " << Ref_legendre.getOrder() <<" on reference element [-1,1] : " << std::endl;
    // std::cout << std::endl;
    // for(auto ii=0; ii<=order; ii++)
    // {
    //     std::cout << "hat_phi_" << ii <<"(x) = " << Ref_legendre.getBasis(ii) << std::endl;
    // }

    // std::cout << std::endl;

    // referenceElement<tchebychevPolynomial> Ref_tchebychev(order);

    // std::cout << "Tchebychev basis of order " << Ref_tchebychev.getOrder() <<" on reference element [-1,1] : " << std::endl;
    // std::cout << std::endl;
    // for(auto ii=0; ii<=order; ii++)
    // {
    //     std::cout << "hat_phi_" << ii <<"(x) = " << Ref_tchebychev.getBasis(ii) << std::endl;
    // }

    // std::cout << std::endl;

    // ----------------------------------Test FiniteElement-------------------------------------

    // size_t order = 1;
    // referenceElement<LegendrePolynomial> Ref_legendre(order);
    // double x_left = 0.;
    // double x_right = 1.;
    // size_t n_cells = 10;
    
    // mesh<double> msh(x_left,x_right,n_cells);
    // msh.print();
    // std::cout << std::endl;

    // FiniteElement<LegendrePolynomial,double> FE(order, msh, 0);

    // gauss_integrator gauss(-1,1,nullptr);
    // double SIP_penalty = 10;


    // std::cout << "Left bound of FE : " << FE.getLeft() << std::endl;
    // std::cout << "Right bound of FE : " << FE.getRight() << std::endl;
    // std::cout << std::endl;

    // for(auto ii=0; ii<=order; ii++)
    // {
    //     std::cout << "phi_" << ii << "(0) = " << FE.evaluate(ii,0) << std::endl;
    //     std::cout << "hat_phi_" << ii << "(-1) = " << Ref_legendre.getBasis(ii)(-1) << std::endl;
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // FE.changeCell(1);

    // std::cout << "After cell change :" << std::endl;

    // std::cout << "Left bound of FE : " << FE.getLeft() << std::endl;
    // std::cout << "Right bound of FE : " << FE.getRight() << std::endl;
    // std::cout << std::endl;

    // for(auto ii=0; ii<=order; ii++)
    // {
    //     std::cout << "phi_" << ii << "(0.5) = " << FE.evaluate(ii,0.5) << std::endl;
    //     std::cout << "hat_phi_" << ii << "(-1) = " << Ref_legendre.getBasis(ii)(-1) << std::endl;
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Eigen::MatrixXd stiffness = FE.local_SIP_volumic(gauss);
    // Eigen::MatrixXd boundary_right = FE.local_SIP_boundary_right(SIP_penalty);
    // Eigen::MatrixXd face_left = FE.local_SIP_faces_left();    
    // Eigen::MatrixXd penalty_left = FE.local_SIP_penalty_left(SIP_penalty);
    // Eigen::SparseMatrix<double> mass = FE.local_mass(gauss);

    // std::cout << "Local SIP volumic Stiffness Matrix : " << std::endl;
    // for(auto ii=0; ii<= FE.getOrder(); ii++)
    // {
    //     for(auto jj=0; jj<= FE.getOrder(); jj++)
    //     {
    //         std::cout << stiffness(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl; 

    // std::cout << "Local SIP right boundary matrix : " << std::endl;
    // for(auto ii=0; ii<= FE.getOrder(); ii++)
    // {
    //     for(auto jj=0; jj<= FE.getOrder(); jj++)
    //     {
    //         std::cout << boundary_right(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl; 

    // std::cout << "Local SIP left interface Matrix : " << std::endl;
    // for(auto ii=0; ii< 2 * FE.Dof(); ii++)
    // {
    //     for(auto jj=0; jj< 2 * FE.Dof(); jj++)
    //     {
    //         std::cout << face_left(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
     

    // std::cout << "Local SIP left penalty Matrix : " << std::endl;
    // for(auto ii=0; ii< 2 * FE.Dof(); ii++)
    // {
    //     for(auto jj=0; jj< 2 * FE.Dof(); jj++)
    //     {
    //         std::cout << penalty_left(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl; 

    // FE.changeCell(0);

    // Eigen::MatrixXd penalty_right = FE.local_SIP_penalty_right(SIP_penalty);
    // Eigen::MatrixXd face_right = FE.local_SIP_faces_right();

    // std::cout << "Local SIP right interface Matrix : " << std::endl;
    // for(auto ii=0; ii< 2 * FE.Dof(); ii++)
    // {
    //     for(auto jj=0; jj< 2 * FE.Dof(); jj++)
    //     {
    //         std::cout << face_right(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "Local SIP right penalty Matrix : " << std::endl;
    // for(auto ii=0; ii< 2 * FE.Dof(); ii++)
    // {
    //     for(auto jj=0; jj< 2 * FE.Dof(); jj++)
    //     {
    //         std::cout << penalty_right(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl; 



    // std::cout << "Local Mass Matrix : " << std::endl;
    // for(auto ii=0; ii<= FE.getOrder(); ii++)
    // {
    //     for(auto jj=0; jj<= FE.getOrder(); jj++)
    //     {
    //         std::cout << mass.coeffRef(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // ---------------------------------------Test dGSpace----------------------------------

    // size_t order = 1;
    // referenceElement<LegendrePolynomial> Ref_legendre(order);
    // double x_left = 0.;
    // double x_right = 1.;
    // size_t n_cells = 2;
    // mesh<double> msh(x_left,x_right,n_cells);
    // msh.print();
    // FiniteElement<LegendrePolynomial,double> FE(order, msh, 0);
    // gauss_integrator gauss(-1,1,nullptr);
    // double SIP_penalty = 10;
    
    // dGSpace<LegendrePolynomial,double> FEM(order,msh);
    // Eigen::SparseMatrix<double> global_stiffness = FEM.assembleSIPMatrix(gauss,SIP_penalty);
    // Eigen::SparseMatrix<double> global_mass = FEM.assembleMassMatrix(gauss);
    // Eigen::VectorXd load = FEM.assembleLoadVector(f,gauss);

    // std::cout << "Global SIP stiffness matrix : " << std::endl;
    // for(auto ii=0; ii<FEM.getNumDof(); ii++)
    // {
    //     for(auto jj=0; jj<FEM.getNumDof(); jj++)
    //     {
    //         std::cout << global_stiffness.coeffRef(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // std::cout << "Global mass matrix : " << std::endl;
    // for(auto ii=0; ii<FEM.getNumDof(); ii++)
    // {
    //     for(auto jj=0; jj<FEM.getNumDof(); jj++)
    //     {
    //         std::cout << global_mass.coeffRef(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // std::cout << "Global load vector : " << std::endl;
    // for(auto ii=0; ii < FEM.getNumDof(); ii++)
    // {
    //     std::cout << load(ii) << ", " << std::endl;
    // }



    //---------------------------------Résolution problème de Poisson----------------------------

    // Paramètres :
    double x_left = -10;
    double x_right = 10;
    size_t n_cells = 1000;
    size_t order = 6;
    double penalty = 10.;

    // Conditions de Dirichlet :
    double g_l = 0;
    double g_r = 0;

    // Maillage, espace et quadrature :
    mesh<double> msh(x_left,x_right,n_cells);
    dGSpace<LegendrePolynomial,double> FEM(order,msh);
    gauss_integrator gauss(-1,1,nullptr);

    // Définition du problème :
    Eigen::SparseMatrix<double> StiffnessMatrix = FEM.assembleSIPMatrix(gauss,penalty);
    Eigen::VectorXd LoadVector = FEM.assembleLoadVector(gauss_bell,gauss);
    FEM.applyDirichletCondition(LoadVector,penalty,g_l,g_r);
    // std::cout << "||load|| = " << LoadVector.norm() << std::endl;
    // std::cout << "sum load = " << LoadVector.sum() << std::endl;

    // Résolution :
    std::string method = "LU";
    LinearSolver Solver(StiffnessMatrix,method);
    std::vector<double> solution = Solver.solve(LoadVector);

    // Export :
    output_writer<LegendrePolynomial,double> output(msh,"Poisson");
    output.write_solution(solution,FEM,"solution");


    return 0;
}