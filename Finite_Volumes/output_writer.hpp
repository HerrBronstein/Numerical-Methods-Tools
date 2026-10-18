#pragma once
#include <fstream>
#include "mesh_1d.hpp"
#include "field.hpp"

template<typename T>
class output_writer
{
    public:

    output_writer(mesh<T> const& msh, std::string filename):mesh_(msh),radical_name_(filename){};

    void write_solution(field<T> const& solution, std::string loop_counter);
    void write(std::vector<T> const& f, std::vector<T> const& g, std::string name);
    inline const mesh<T>& getMesh() const {return mesh_;};

    protected:

    mesh<T> const& mesh_;
    std::string radical_name_;
};


template<typename T>
void output_writer<T>::write_solution(field<T> const& uh, std::string loop_counter)
{
    std::string filename = radical_name_ + "_" + loop_counter + ".txt";

    std::ofstream file(filename);

    if (file.is_open())
    {
        for(auto cell=0; cell < mesh_.n_cells(); cell++)
        {
            T x_i = mesh_.xc(cell);
            file << x_i << " " << uh(cell) << "\n";
        }
        file.close();
        std::cout << "Fichier " << filename << " enregistre avec succes." << std::endl;
    }
    else
    {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << std::endl;
    }
}

template<typename T>
void output_writer<T>::write(std::vector<T> const& f, std::vector<T> const& g, std::string name)
{
    std::string filename = radical_name_ + "_" + name + ".txt";

    std::ofstream file(filename);

    if (file.is_open())
    {
        for(auto ii=0; ii < f.size(); ii++)
        {
            file << f[ii] << " " << g[ii] << "\n";
        }
        file.close();
        std::cout << "Fichier " << filename << " enregistre avec succes." << std::endl;
    }
    else
    {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << std::endl;
    }
}