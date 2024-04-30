#include "matrix.hpp"

/*!
* Definition of compress and uncompress methods and const and non-const call operators.
*/

/*!
* Compress method: looping all along the map, with the ordering that will depend on the storaging, to retrieve the nnz elements.
*/
template<class T, StorageOrder O>
void
algebra::Matrix<T,O>::compress(){

    //matrix has to be uncompressed
    if (!(this->is_compressed()))
    {   
        //storing the nonzero values:
        //cannot use std::ranges::copy since elements of map are pair: transform needed
        m_val_comp.reserve(m_nnz);
        std::ranges::transform(m_mat_uncomp,std::back_inserter(m_val_comp),[](std::pair<key_type,T> p){return p.second;});
        
        if constexpr(O==StorageOrder::RowWise)  //RowWise
        {   
            //storing the columns in which there are non zero elements
            m_outer.reserve(m_nnz);
            std::ranges::transform(m_mat_uncomp,std::back_inserter(m_outer),[](std::pair<key_type,T> p){return p.first[1];});

            //storing for each row how many cumulative nnz elements there are
            m_inner.reserve(m_m+1);
            m_inner.push_back(0);    //first element is always 0
            for (size_t i = 0; i < m_m; ++i)
            {   
                m_inner.push_back(std::ranges::count_if(m_mat_uncomp,[&i](std::pair<key_type,T> p){return p.first[0]<=i;}));
            }
        }   
        else                                    //ColumnWise
        {
            //storing the rows in which there are non zero elements
            m_outer.reserve(m_nnz);
            std::ranges::transform(m_mat_uncomp,std::back_inserter(m_outer),[](std::pair<key_type,T> p){return p.first[0];});

            //now I have to store for each column how many cumulative nnz elements I have
            m_inner.reserve(m_n+1);
            m_inner.push_back(0);    //first element is always 0
            for (size_t i = 0; i < m_n; ++i)
            {   
                m_inner.push_back(std::ranges::count_if(m_mat_uncomp,[&i](std::pair<key_type,T> p){return p.first[1]<=i;}));
            }
        }   
        this->is_compressed() = true;   //matrix has been compressed
        this->clear_buffer();           //free memory       
    }
}


/*!
* Uncompress method: unsqueezing the vector to fill the map. Storage order will take account of the ordering in the map.
*/ 
template<class T, StorageOrder O>
void
algebra::Matrix<T,O>::uncompress() {

    //matrix has to be compressed
    if (this->is_compressed())
    {
        if constexpr(O==StorageOrder::RowWise){             //RowWise
            
            //checking every row
            for (size_t i = 0; i < m_m; ++i)
            {   
                //checking if there is any element in that row
                if (m_inner[i]<m_inner[i+1])    
                {   
                    //indeces of row i-th in values and outer
                    auto begin_row_i = m_inner[i];
                    std::vector<size_t> col_idx;
                    col_idx.resize(m_inner[i+1]-m_inner[i]);
                    //indeces of the elements in m_val_comp and m_outer that are respective for row i-th:
                    //their values are delimited by the values in m_inner
                    std::iota(col_idx.begin(),col_idx.end(),begin_row_i);

                    //inserting in the map(no check needed for presence)
                    //no std::ranges version in order to relay on parallel operations (since the operation is on a std::vector)
                    std::for_each(std::execution::par,
                                  col_idx.begin(),col_idx.end(),
                                  [&i,this](std::size_t j){ key_type key_j{i,this->m_outer[j]}; this->m_mat_uncomp.insert({key_j,this->m_val_comp[j]});}
                                  );
                    //clearing
                    col_idx.clear();
                }
            } 
        }  
        else                                                //ColumnWise
        {
            //same code(different names for interpretability), except for the order of inserting in the pair for the map key
            //checking every col
            for (size_t i = 0; i < m_n; ++i)
            {   
                //checking if there is any element in that col
                if (m_inner[i]<m_inner[i+1])    
                {   
                    //indeces of col i-th in values and outer
                    auto begin_col_i = m_inner[i];
                    std::vector<size_t> row_idx;
                    row_idx.resize(m_inner[i+1]-m_inner[i]);
                    std::iota(row_idx.begin(),row_idx.end(),begin_col_i);

                    //inserting in the map(no check needed for presence)
                    //no std::ranges version in order to relay on parallel operations
                    std::for_each(std::execution::par,
                                  row_idx.begin(),row_idx.end(),
                                  [&i,this](std::size_t j){ key_type key_j{this->m_outer[j],i}; this->m_mat_uncomp.insert({key_j,this->m_val_comp[j]});}
                                  );
                    //clearing
                    row_idx.clear();
                }
            } 

        }
        this->is_compressed()=false;    //matrix has been compressed
        this->clear_buffer();           //free memory
    } 
}


/*!
Const call operator: if there is no element, 0 casted to type T is returned. Else the value of elem (i,j).
*/
template<class T, StorageOrder O>
T
algebra::Matrix<T,O>::operator()(std::size_t i, std::size_t j) const
{   
    //checking if the element is present: if not, 0 of the right type is returned
    if (!this->check_presence(i,j)) {return (static_cast<T>(0));}

    //uncompressed
    if(!(this->is_compressed())){                   //way of storaging is not important since the ordering of the map will take care of everything
        key_type elem{i,j};
        return m_mat_uncomp.find(elem)->second;  //.find() method exploited 
    }

    //compressed
    if constexpr(O == StorageOrder::RowWise)        //RowWise storage
    {   
        //using std::find instead of std::ranges::find to relay on parallelization
        //range considered in the find: only one value of j in the range defined by
        //iterator at the beginning of the cols in a row and
        //sentinel at the element past the last one at the cols of row i-th
        auto elem = std::find(std::execution::par,
                              m_outer.cbegin() + m_inner[i],m_outer.cbegin() + m_inner[i+1],
                              j);
        //taking the distance from the begin to retrain the position in the values' vector
        return m_val_comp[std::distance(m_outer.cbegin(),elem)];        
    }
    else                                            //ColWise storage
    {   
        //same way of proceeding as above, but using the col index
        auto elem = std::find(std::execution::par,
                              m_outer.cbegin() + m_inner[j],m_outer.cbegin() + m_inner[j+1],
                              i);
        //taking the distance from the begin to retrain the position in the value's vector
        return m_val_comp[std::distance(m_outer.cbegin(),elem)];         
    }
}

/*!
Non-const call operator: if compressed, cannot add/remove elements, but only modify already present nnz elements. If uncompressed, can modify, add or remove elements.
*/
template<class T, StorageOrder O>
T &
algebra::Matrix<T,O>::operator()(std::size_t i, std::size_t j, const T newvalue)
{   
    //checking if the indexes passed are ok
    assert(i<m_m and j<m_n);
    //checking if the element is already present
    bool check_ij = check_presence(i,j);
    //checking if newvalue is null
    T zerovalue = 0;
    bool new_zero = (newvalue==zerovalue);
    //No non-null elements and newvalue null means nothing to be done
    if (new_zero and m_nnz==0){return zerovalue;}

    //uncompressed
    if (!(this->is_compressed()))
    {
        key_type key_elem{i,j};
        if(new_zero){   //adding a 0 means removing an element
            m_nnz--;
            m_mat_uncomp.erase(key_elem);
            return zerovalue;
        } else{         //add if not already present, if not substitute: not using .insert_or_assign() method since the check has already been done
        if (check_ij){
            m_mat_uncomp[key_elem]=newvalue;
        } else{
            m_nnz++;
            m_mat_uncomp.insert({key_elem,newvalue});
        }              
        return m_mat_uncomp[key_elem];
        }
    }
    
    //compressed
    //Cannot add/remove an element in compress state: 0 is returned
    if (this->is_compressed() and (!check_ij or new_zero))
    {
        std::cerr<<"Cannot add/remove an element in compressed state"<<std::endl;
        return zerovalue;
    }
    //modify only   
    if constexpr(O == StorageOrder::RowWise)            //RowWise storage
    {       
            //finding the position in a constrained range, amd then modifying the correspective value
            auto elem = std::find(std::execution::par,
                                  m_outer.cbegin() + m_inner[i],m_outer.cbegin() + m_inner[i+1],
                                  j);
            auto index = std::distance(m_outer.cbegin(),elem);
            m_val_comp[index] = newvalue;
            return m_val_comp[index];      
    }   else                                            //ColWise storage
    {
        auto elem = std::find(std::execution::par,
                              m_outer.cbegin() + m_inner[j],m_outer.cbegin() + m_inner[j+1],
                              i);
        auto index = std::distance(m_outer.cbegin(),elem);
        m_val_comp[index] = newvalue;
        return m_val_comp[index]; 
    }
}