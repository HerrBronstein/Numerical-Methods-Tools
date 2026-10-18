#pragma once
#include <vector>

template<typename T>
class mesh
{
    public:

    mesh(T x_left, T x_right, size_t n_cells);
    ~mesh(){};

    inline const size_t vertices() const {return n_vertices_;};
    inline const size_t n_cells() const {return n_cells_;};
    inline const T dx() const {return dx_;};
    size_t getMesh(T x) const;
    T vertice_coord(size_t ii) const {return vertices_coordinates_[ii];};

    inline T xc(size_t ii) const {return cells_centers_[ii];};
    void print();


    protected:
    
    std::vector<T> vertices_coordinates_;
    std::vector<T> cells_centers_;
    size_t n_vertices_{0};
    size_t n_cells_{0};
    T dx_;
};

template<typename T>
mesh<T>::mesh(T x_left, T x_right, size_t n_cells)
{
    dx_ = (x_right - x_left)/n_cells;
    n_cells_ = n_cells;
    n_vertices_ = n_cells + 1;

    cells_centers_.resize(n_cells_);
    vertices_coordinates_.resize(n_vertices_);

    for(auto ii=0; ii<n_cells_; ii++)
    {
        cells_centers_[ii] = x_left + (ii + 0.5)*dx_;
    }
    for(auto ii=0; ii<n_vertices_; ii++)
    {
        vertices_coordinates_[ii] = x_left + ii*dx_;
    }
}

template<typename T>
void mesh<T>::print()
{
    std::cout << "|";
    for(auto ii=0; ii<n_cells_; ii++)
    {
        std::cout << vertices_coordinates_[ii] << "[---(" << cells_centers_[ii] << ")---]";
    }
    std::cout << vertices_coordinates_[n_cells_] << "|" << std::endl;
}

template<typename T>
size_t mesh<T>::getMesh(T x) const
{
    for(auto ii=0; ii < n_vertices_; ii++)
    {
        if(x >= vertices_coordinates_[ii] && x < vertices_coordinates_[ii+1])
        {
            return ii;
        }
    }

    throw std::out_of_range("Point x is outside the mesh");
}