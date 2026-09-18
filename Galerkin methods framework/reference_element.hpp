#pragma once
#include <iostream>
#include <vector>


template<typename P>
class referenceElement
{
    public:

    referenceElement(size_t order);
    ~referenceElement(){};

    const P& getBasis(size_t ii) const {return basis_[ii];};
    const size_t& getOrder() const {return order_;};
    size_t& getOrder(){return order_;};

    private:

    size_t order_;
    std::vector<P> basis_;
};


template<typename P>
referenceElement<P>::referenceElement(size_t order):order_(order)
{
    for(auto ii=0; ii<=order; ii++)
    {
        basis_.push_back(P(ii));
    }
}
