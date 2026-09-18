#include <iostream>
#include "FESpace.hpp"
#include "solver.hpp"
#include "output_writer.hpp"

double f(double x){return std::exp(-x*x);}

int main()
{

    //---------------------------TEST POLYNOMES DE LAGRANGE--------------------------------
    // int Nint = 4;
    // std::vector<double> x_val(Nint);
    // std::vector<double> y_val(Nint);
    // double x = 1.4;

    // for(auto ii=0; ii<Nint; ii++)
    // {
    //     x_val[ii] = 0.5 + ii/4.;
    //     y_val[ii] = std::exp(x_val[ii]);
    // }

    // LagrangePolynomial p(x_val,y_val);

    // std::cout << "Interpolation of exp(" << x <<") : " << p(x) << std::endl;
    // std::cout << "Real value of exp("<<x<<") : " << std::exp(x) << std::endl;
    // std::cout << "p(x) = " << p << std::endl;

    // p.addPoint(2.,std::exp(2.));

    // std::cout << "Interpolation of exp(" << x <<") : " << p(x) << std::endl;
    // std::cout << "Real value of exp("<<x<<") : " << std::exp(x) << std::endl;
    // std::cout << "p(x) = " << p << std::endl;

    // int order = 2;

    //----------------------------------TEST BASE LAGRANGE----------------------------------
    
    // basis phi(order);

    // std::cout << "phi_0(x) = " << phi[0] << std::endl;
    // std::cout << "phi_1(x) = " << phi[1] << std::endl;
    // std::cout << "phi_2(x) = " << phi[2] << std::endl;
    // std::cout << "phi_0(-1) = " << phi[0](-1) << std::endl;
    // std::cout << "phi_0(0) = " << phi[0](0) << std::endl;
    // std::cout << "phi_1(0) = " << phi[1](0) << std::endl;
    // std::cout << "phi_2(1) = " << phi[2](1) << std::endl;

    // double x_left = -1;
    // double x_right = 1;
    // size_t n_cells = 2;

    // mesh<double> msh(x_left, x_right, n_cells);

    // msh.print();

    // double tab[10] = {1.,2.,5.,7.,3.,9.,0.,88.,28.,61.};

    // polynomial p1(15,0.23);
    // polynomial p2(9,tab);

    // std::cout << "p1(x) = " << p1 << std::endl;
    // std::cout << "p2(x) = " << p2 << std::endl;

    // std::cout << "p1 * p1 (x) = " << p1*p2 << std::endl;

    //-------------------------------------TEST ELEMENT FINI-----------------------------------

    // double x_left = 0.;
    // double x_right = 1.;
    // size_t n_cells = 1;
    // int order = 0;

    // mesh<double> msh(x_left,x_right,n_cells);

    // msh.print();

    // FiniteElement<double> FE(msh,0,order);

    // std::cout << "FE_left = " << FE.getLeft() << std::endl;
    // std::cout << "FE_right = " << FE.getRight() << std::endl;
    // std::cout << "FE.dx() = " << FE.dx() << std::endl;
    // std::cout << std::endl;

    // std::vector<LagrangePolynomial> poly = FE.getPoly();

    // for(auto ii=0; ii<poly.size(); ii++)
    // {
    //     for(auto jj=0; jj<poly.size(); jj++)
    //     {
    //         std::cout << "phi_"<< ii <<"("<< FE.getNode(jj) <<") = " << FE.evaluate(ii,FE.getNode(jj)) << std::endl;
    //     }
    // }
    // std::cout << std::endl;

    // basis base(order);

    // for(auto ii=0; ii<base.getSize(); ii++)
    // {
    //     for(auto jj=0; jj<base.getSize(); jj++)
    //     {
    //         std::cout << "hat_phi_"<< ii <<"("<< base.getNode(jj) <<") = " << base[ii](base.getNode(jj)) << std::endl;
    //     }
    // }
    // std::cout << std::endl;

    // double x = 0.8;

    // size_t cell = msh.getMesh(x);
    // std::cout << x << " is located in cell number "<< cell << " : [" << msh.vertice_coord(cell) << ", " << msh.vertice_coord(cell+1) << "]" << std::endl;
    // std::cout << std::endl;

    // gauss_integrator gauss(-1,1,nullptr);
    // Eigen::MatrixXd stiffness = FE.local_stiffness(gauss);
    // Eigen::MatrixXd mass = FE.local_mass(gauss);
    // size_t n = FE.getDof();

    // for(auto ii=0; ii<n; ii++)
    // {
    //     for(auto jj=0; jj<n; jj++)
    //     {
    //         std::cout << stiffness(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;
    // for(auto ii=0; ii<n; ii++)
    // {
    //     for(auto jj=0; jj<n; jj++)
    //     {
    //         std::cout << mass(ii,jj) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Eigen::VectorXd load = FE.local_load(f,gauss);

    // for(auto ii=0; ii<n; ii++)
    // {
    //     std::cout << load(ii) << ", " << std::endl;
    // }

    //-------------------------------------TEST FESpace--------------------------------------

    // double x_l = 0;
    // double x_r = 1;
    // int Nint = 2;
    // mesh<double> msh(x_l,x_r,Nint);

    // int order = 1;

    // FESpace<double> FEM(msh,order);

    // size_t dof = FEM.numDof();
    // gauss_integrator gauss(-1,1,nullptr);

    // Eigen::MatrixXd stiffness = FEM.assembleStiffnessMatrix(gauss);
    // Eigen::VectorXd load = FEM.assembleLoadVector(f,gauss);


    // FEM.applyDirichletCondition(stiffness,load,0,0);

    // for(auto ii=0; ii<dof; ii++)
    // {
    //     for(auto jj=0; jj<dof; jj++)
    //     {
    //         std::cout << stiffness(ii,jj) << ", ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;
    // for(auto ii=0; ii<dof; ii++)
    // {
    //     std::cout << load(ii) << ", " << std::endl;
    // }
    // std::cout << std::endl;
    

    //--------------------------------------------------------------------------------------
    //-------------------------------TEST PROBLEME POISSON----------------------------------
    //--------------------------------------------------------------------------------------
    //Paramètres :
    double x_left = -20;
    double x_right = 20;
    size_t n_cells = 2000;
    int order = 1;

    double left_val = 0;
    double right_val = 0;
    std::string method = "QR";
    //-------------------------------------------------------------------------------------- 
    // Maillage, espace, intégrateur et affichage :
    mesh<double> msh(x_left,x_right,n_cells);
    FESpace<double> FEM(msh,order);
    gauss_integrator gauss0(-1,1,nullptr);
    output_writer<double> writer(msh,"Poisson");
    //--------------------------------------------------------------------------------------
    //Assemblage du système et conditions limites :
    Eigen::MatrixXd A = FEM.assembleStiffnessMatrix(gauss0);
    Eigen::VectorXd B = FEM.assembleLoadVector(f,gauss0);
    FEM.applyDirichletCondition(A,B,left_val,right_val);
    //--------------------------------------------------------------------------------------
    // Résolution :
    LinearSolver solver(A,B,method);
    std::vector<double> solution = solver.solve();
    //--------------------------------------------------------------------------------------
    //Export :
    writer.write_solution(solution,FEM,"0");
    
    

    return 0;
}