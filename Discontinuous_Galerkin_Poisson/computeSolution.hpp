#pragma once
#include "dGSpace.hpp"
#include "output_writer.hpp"
#include "solver.hpp"


template<typename P, typename T>
std::vector<T> computeSolution(dGSpace<P,T>& FEM, std::function<double(double)> load, std::string method, double dirichlet_L, double dirichlet_R)
{
    double penalty = 20.;

    // dGSpace<P,T> FEM(order,msh);
    gauss_integrator gauss(-1,1,nullptr);

    // Définition du problème :
    Eigen::SparseMatrix<double> StiffnessMatrix = FEM.assembleSIPMatrix(gauss,penalty);
    Eigen::VectorXd LoadVector = FEM.assembleLoadVector(load,gauss);
    FEM.applyDirichletCondition(LoadVector,penalty,dirichlet_L,dirichlet_R);
    // std::cout << "||load|| = " << LoadVector.norm() << std::endl;
    // std::cout << "sum load = " << LoadVector.sum() << std::endl;

    // Résolution :
    LinearSolver Solver(StiffnessMatrix,method);
    std::vector<double> solution = Solver.solve(LoadVector);

    return solution;
}