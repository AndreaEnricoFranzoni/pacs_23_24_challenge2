#include "matrix.hpp"

/*!
* Definition of the method norms, in all of its three specializations: One, Infinity, Frobenius.
* Specialization is done using tag dispatching due to std::integral_constant.
*/

/*!
* Norm One: max within the cols of the sum of the absolute values of every element of the col.
*/
template<class T, StorageOrder O>
auto
algebra::Matrix<T,O>::norm(NT<NormType::One>) const{

    //if the matrix has no elements at all: return 0 casted to a double
    if (this->m()==0 and this->n()==0)
    {
        return static_cast<double>(0);
    }
    
    //storing the sum of each col
    std::vector<double> tmp;
    tmp.resize(m_n);

    //checking all the cols
    std::iota(tmp.begin(),tmp.end(),0);

    //for each col, in parallel, the sum of the absolute values of the elements of the col is done
    //and then stored in the temporary container
    std::transform(std::execution::par,
                   tmp.begin(),
                   tmp.end(),
                   tmp.begin(),
                   [this](auto el)
                   {    auto row_i = this->get_col(el).second; 
                        return std::transform_reduce(std::execution::par,row_i.cbegin(),row_i.cend(),0.0,std::plus{},[](auto i){return std::abs(i);}); 
                    } 
                    );
    
    //returning the max
    return std::ranges::max(tmp);
}


/*!
* Norm Infinity: max within the rows of the sum of the absolute values of every element of the row.
*/
template<class T, StorageOrder O>
auto
algebra::Matrix<T,O>::norm(NT<NormType::Infinity>) const{

    //if the matrix has no elements at all: return 0 casted to a double
    if (this->m()==0 and this->n()==0)
    {
        return static_cast<double>(0);
    }
    
    //storing the sum of each row
    std::vector<double> tmp;
    tmp.resize(m_m);

    //checking all the rows
    std::iota(tmp.begin(),tmp.end(),0);

    //for each row, in parallel, the sum of the absolute values of the elements of the row is done
    //and then stored in the temporary container
    std::transform(std::execution::par,
                   tmp.begin(),
                   tmp.end(),
                   tmp.begin(),
                   [this](auto el)
                   {    auto col_i = this->get_row(el).second; 
                        return std::transform_reduce(std::execution::par,col_i.cbegin(),col_i.cend(),0.0,std::plus{},[](auto i){return std::abs(i);}); 
                    } 
                    );
    
    //returning the max
    return std::ranges::max(tmp);
}


/*!
* Norm Frobenius: squared root of the sum of the squared absolute values of each element of the matrix.
*/
template<class T, StorageOrder O>
auto
algebra::Matrix<T,O>::norm(NT<NormType::Frobenius>) const{
    
    //if the matrix has no elements at all: return 0 casted to a double
    if (this->m()==0 and this->n()==0)
    {
        return static_cast<double>(0);
    }
    
    //compressed: relaying on parallel since values are in a sequential container (std::vector)
    if (this->is_compressed())
    {
        auto res = std::transform_reduce(std::execution::par,
                                         m_val_comp.cbegin(),
                                         m_val_comp.cend(),
                                         0.0,
                                         std::plus{},
                                         [](T el){return std::pow(std::abs(el),2);}
                                         );
        return std::sqrt(res);
        
    }
    //uncompressed: not relaying on parallel since values are in an associative container (std::map)
    else
    {
        
        auto res = std::transform_reduce(m_mat_uncomp.cbegin(),
                                        m_mat_uncomp.cend(),
                                        0.0,
                                        std::plus{},
                                        [](auto el){return std::pow(std::abs(el.second),2);}
                                        );
        return std::sqrt(res);
    }
}