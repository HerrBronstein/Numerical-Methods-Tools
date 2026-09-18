#pragma once
#include <eigen/Dense>
#include <eigen/Sparse>
#include <vector>
#include <string>
#include <stdexcept>

class LinearSolver
{
    public:

    LinearSolver(Eigen::SparseMatrix<double> const& A, std::string method):A_(A),method_(method){};
    ~LinearSolver(){};

    void changeMethod(std::string const& newMethod){method_ = newMethod;};
    void changeMatrix(Eigen::SparseMatrix<double> const& A){A_ = A;};

    Eigen::SparseMatrix<double>& matrix(){return A_;};
    std::string& method(){return method_;};

    std::vector<double> solve(Eigen::VectorXd& B) const;

    private:

    Eigen::SparseMatrix<double> A_;
    std::string method_;
};


std::vector<double> LinearSolver::solve(Eigen::VectorXd& B) const
{
    Eigen::VectorXd vec;
    if(method_ == "LU")
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A_);
        vec = solver.solve(B);
    }
    else if(method_ == "Choleski")
    {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A_);
        vec = solver.solve(B);
    }
    else
    {
        throw std::invalid_argument("Unknown solver method");
    }

    std::vector<double> sol(vec.data(), vec.data() + vec.size());

    return sol;
}