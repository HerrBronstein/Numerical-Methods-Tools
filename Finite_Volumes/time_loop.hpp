#pragma once
#include "field.hpp"
#include "residual.hpp"

template<typename T>
class time_loop
{
    public:

    time_loop(mesh<T> const& msh, physical_flux<T> const& flux, numerical_flux<T> const& num_flux, T final_time)
    :mesh_(msh),final_time_(final_time),residual_(residual<T>(msh)),flux_(flux),num_flux_(num_flux){};

    ~time_loop(){};

    inline const T final_time() const {return final_time_;};
    void run(field<T>& uh);
    void compute_dt(field<T>& uh);

    protected:

    size_t iteration_counter_{0};
    T dt_{0.};
    T physical_time_{0.};
    T final_time_{0.};
    bool last_iteration_{false};
    double cfl_number_ = 0.5 ;

    residual<T> residual_;

    mesh<T> const& mesh_;

    physical_flux<T> const& flux_;
    numerical_flux<T> const& num_flux_;


};


template<typename T>
void time_loop<T>::compute_dt(field<T>& uh)
{
    T max = 0;
    for(auto ii=0; ii<mesh_.n_cells(); ii++)
    {
        if( max <= std::abs( flux_.derivative( uh(ii) ) ) ) 
        {
            max = std::abs( flux_.derivative( uh(ii) ) );
        }
    }

    dt_ = cfl_number_ * mesh_.dx() / max;
}


template<typename T>
void time_loop<T>::run(field<T>& uh)
{

    while( physical_time_ < final_time() )
    {
        compute_dt(uh);
        if( physical_time_ + dt_ >= final_time() )
        {
            dt_ = final_time() - physical_time_;
            last_iteration_ = true;
        }

        residual_.assemble_from_field(uh,flux_,num_flux_);
        uh += dt_ * residual_;

        physical_time_ += dt_;
        iteration_counter_ += 1;

        if(last_iteration_) break;

    }

}