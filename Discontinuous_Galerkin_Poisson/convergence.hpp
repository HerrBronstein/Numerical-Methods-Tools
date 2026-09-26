#include "output_writer.hpp"
#include "dGSpace.hpp"
#include <vector>
#include <functional>
#include "solver.hpp"
#include "computeSolution.hpp"


template<typename P, typename T>
class convergence
{
    public:

    convergence(std::function<double(double)> const& u_exact, mesh<T> msh, size_t order, double g_l, double g_r, std::function<double(double)> const& load, std::string method, size_t nb_raf) 
        :u_exact_(u_exact),mesh_(msh),order_(order),g_l_(g_l),g_r_(g_r),load_(load),method_(method),nb_raf_(nb_raf){};

    ~convergence(){};

    void compute_L2_error();
    // void compute_H1_error();

    const std::vector<double> getL2error() const {return L2err_;};
    const std::vector<double> getH1error() const {return H1err_;};

    // void export_error();

    private:

    const std::function<double(double)>& u_exact_;
    const std::function<double(double)>& load_;
    std::string method_;
    mesh<T> const& mesh_;
    size_t order_;
    double g_l_;
    double g_r_;
    const size_t nb_raf_;

    std::vector<double> L2err_;
    std::vector<double> H1err_;
};


template<typename P,typename T>
void convergence<P,T>::compute_L2_error()
{
    T xl = mesh_.vertice_coord(0);
    T xr = mesh_.vertice_coord(mesh_.n_cells());
    std::vector<T> h_raf;
    gauss_integrator integ_(-1,1,nullptr);

    for(auto ii=0; ii<=nb_raf_; ii++)
    {
        double error = 0.;
        size_t n_cell = mesh_.n_cells() * std::pow(2,ii);
        mesh<T> raf_mesh = mesh<T>(xl,xr,n_cell);
        T dx = raf_mesh.dx();

        dGSpace<P,T> space(order_,raf_mesh);

        double J = 0.5 * dx;
    
        std::vector<T> uh = computeSolution(space,load_,method_,g_l_,g_r_);

        for(auto cell=0; cell < raf_mesh.n_cells(); cell++)
        {
            T x_l = raf_mesh.vertice_coord(cell);
            T x_r = raf_mesh.vertice_coord(cell+1);

            std::function<double(double)> diff = [this,uh, J, x_l, x_r, dx, &space, cell](double xi)
            {
                double x = 0.5*(x_r + x_l) + 0.5 * dx * xi;
                double e = u_exact_(x) - space.evaluate_projection(uh,cell,x);
                return J * e * e ;
            };

            integ_.change_integrand(diff);
            error += integ_.compute_integral(1);
            
        }
        double err = std::sqrt(error);
        std::cout << "For h = "<< dx <<", || u - uh ||_L^2 = " << err << std::endl;
        L2err_.push_back(err);
        h_raf.push_back(dx);
    }
    output_writer<P,T> writer(mesh_,"Poisson");

    writer.write(h_raf,L2err_,"L2_error");
}

