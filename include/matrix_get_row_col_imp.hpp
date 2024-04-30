#include "matrix.hpp"

/*!
* Definition of check_presence_row, check_presence_col, check_presence, get_row and get_col method.
*/

/*!
* Checking presence of row idx-th: at least a nnz element in that row.
*/
template<class T, StorageOrder O>
bool
algebra::Matrix<T,O>::check_presence_row(const std::size_t &idx) const
{
    assert(idx<m_m);            //checking if passed index is ok wrt the dimension of matrix

    if constexpr(O==StorageOrder::RowWise)
    {
        key_type lb{idx,0};             //lower bound of the row
        key_type up{idx+1,0};           //upper bound of the row
        //if compressed: is sufficient that between row i-th and i-th there is an increment in m_inner
        //if uncompressed: relaying on upper_bound and lower_bound methods. To handle the element at the beginning of the row is necessary to check that specific row
        return (this->is_compressed()) ? (m_inner[idx]<m_inner[idx+1]) : (!(std::distance(m_mat_uncomp.upper_bound(lb),m_mat_uncomp.lower_bound(up))==0) or m_mat_uncomp.contains(lb));
    }
    else
    {   
        return (this->is_compressed()) ? (std::find(std::execution::par,m_outer.cbegin(),m_outer.cend(),idx)!=m_outer.cend()) : (std::ranges::find_if(m_mat_uncomp,[&idx](auto x){return x.first[0]==idx;})!=m_mat_uncomp.cend());
        //if compressed: is sufficient that in the vector containing the rows there is once idx: search can be parallelized
        //if uncompressed: is sufficient finding a key such that the first element of it is idx: cannot be parallelized
    }
}


/*!
* Checking presence of col idx-th: at least a nnz element in that col.
*/
template<class T, StorageOrder O>
bool
algebra::Matrix<T,O>::check_presence_col(const std::size_t &idx) const
{
    assert(idx<m_n);            //checking if passed index is ok wrt the dimension of matrix

    if constexpr(O==StorageOrder::RowWise)
    {
        return (this->is_compressed()) ? (std::find(std::execution::par,m_outer.cbegin(),m_outer.cend(),idx)!=m_outer.cend()) : (std::ranges::find_if(m_mat_uncomp,[&idx](auto x){return x.first[1]==idx;})!=m_mat_uncomp.cend());
        //if compressed: is sufficient that in the vector containing the rows there is once idx: search can be parallelized 
        //if uncompressed: is sufficient finding a key such that the first element of it is idx: cannot be parallelized
    }   else
    {   
        //same considerations of check_row in RowWise storaging
        key_type lb{0,idx};
        key_type up{0,idx+1};
        return (this->is_compressed()) ? (m_inner[idx]<m_inner[idx+1]) : (!(std::distance(m_mat_uncomp.upper_bound(lb),m_mat_uncomp.lower_bound(up))==0) or m_mat_uncomp.contains(lb));
    }
}


/*!
* Checking for the presence of the element in position (i,j).
*/
template<class T, StorageOrder O>
bool
algebra::Matrix<T,O>::check_presence(const std::size_t &i, const std::size_t &j) const
{
    assert(i<m_m and j < m_n);
    //checking if there is any element in that row and in that col
    if (!this->check_presence_row(i) or !this->check_presence_col(j)){
        return false;
    }
    //checking in uncompressed state, the operator less will take account for it
    if (!this->is_compressed())
    {
        key_type elem{i,j};
        return m_mat_uncomp.contains(elem);
    }
    //checking in compressed state, and this depends on the storaging
    if constexpr(O==StorageOrder::RowWise)
    {   //searching if within the column indeces of the row there is the one requested
        auto start_search = m_outer.cbegin() + m_inner[i];
        auto sentinel_search = m_outer.cbegin() + m_inner[i+1];
        return std::binary_search(start_search,sentinel_search,j);
    }   
    else
    {   //searching if within the row indeces of the col there is the one requested
        auto start_search = m_outer.begin() + m_inner[j];
        auto sentinel_search = m_outer.begin() + m_inner[j+1];
        return std::binary_search(start_search,sentinel_search,i);
    }
}


/*!
* Getting row idx-th.
*/
template<class T, StorageOrder O>
std::pair<std::vector<std::size_t>,std::vector<T>>
algebra::Matrix<T,O>::get_row(const std::size_t &idx) const{

    //not checking the presence: it will not store the index: ok
    //functions that uses it will check for row presence

    //storing col indeces of the row
    std::vector<std::size_t> col_r_idx;         
    //storing values of the row
    std::vector<T> val_r_idx;
    //number of elements in the row
    std::size_t nr;

    if constexpr (O == StorageOrder::RowWise)                           //RowWise
    {
        if (this->is_compressed())                                      //compressed
        {   
            //iterators to take the range of the row in the 
            //compressed representation, for col indeces and values
            auto start_col = m_outer.cbegin() + m_inner[idx];           
            auto sentinel_col = m_outer.cbegin() + m_inner[idx+1];
            auto start_val = m_val_comp.cbegin() + m_inner[idx];
            auto sentinel_val = m_val_comp.cbegin() + m_inner[idx+1];

            nr = std::distance(start_col,sentinel_col);

            //memory reserved
            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            //copy of the value: done in parallel since operating on std::vector
            std::copy(std::execution::par,
                      start_col,sentinel_col,std::back_inserter(col_r_idx));
            std::copy(std::execution::par,
                      start_val,sentinel_val,std::back_inserter(val_r_idx));
        }   
        else                                                            //uncompressed
        {   
            //problem if looking for the first row, since upper bound retrieve the first greater
            //and not greater_equal: cannot give {-1,0} as key, since these elements are size_t, it is necessary to
            //take account for the element (0,0)
            if (idx==0)
            { 
                key_type ub{idx+1,0};                                   //upper bound of the row
                auto start_r = m_mat_uncomp.cbegin();                   //beginning of the map
                auto sentinel_r = m_mat_uncomp.lower_bound(ub);         //sentinel for the row

                nr = std::distance(start_r,sentinel_r);                 //number of the elements in the row
                col_r_idx.reserve(nr);                                  //reserving memory
                val_r_idx.reserve(nr);

                //putting elements in the right vector
                std::for_each(start_r,sentinel_r,[&col_r_idx,&val_r_idx](auto p){col_r_idx.push_back(p.first[1]); val_r_idx.emplace_back(p.second);});
            }
            else
            {
                key_type lb{idx-1,this->n()+1};                         //lower bound of the row
                key_type ub{idx+1,0};                                   //upper bound of the row
                auto start_r = m_mat_uncomp.upper_bound(lb);
                auto sentinel_r = m_mat_uncomp.lower_bound(ub);         //sentinel for the row

                nr = std::distance(start_r,sentinel_r);                 //number of the elements in the row
                col_r_idx.reserve(nr);                                  //reserving memory for the storaging
                val_r_idx.reserve(nr);

                //putting elements in the right vector
                std::for_each(start_r,sentinel_r,[&col_r_idx,&val_r_idx](auto p){col_r_idx.push_back(p.first[1]); val_r_idx.emplace_back(p.second);});
            }
        }
    }   
    else                                                                //ColWise
    {   

        if (this->is_compressed())                                      //compressed
        {   
            //number of element in the row: can be done in parallel since it is std::vector
            nr = std::count(std::execution::par,                        
                            m_outer.cbegin(),m_outer.cend(),
                            idx);

            //memory reserving
            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            //variable for looping
            std::size_t k{0};
            //idea: -the value can be retrieved triavially since it is only matter of 
            //       accessing the position of m_value_comp in which m_outer has idx as element
            //      -the indeces of the cols can be retrieved as:
            //       the j-th elem of m_inner contains the number of elem up to the col j-1-th
            //       so, looping over all the cols, starting from the last one in which an elem is
            //       found, until the corrispettive elem in m_inner is bigger of the position
            //       of the elem in m_outer, that indicate which elems are in the idx-th row
            //Looping until all the elemnts have been found
            std::size_t check_point_row_cols{0};

            while (val_r_idx.size()<nr)
            {
                if (m_outer[k]==idx)
                {   
                    //pushing the value
                    val_r_idx.emplace_back(m_val_comp[k]);

                    //searching the first col available as
                    auto col_ptr = std::find_if(m_inner.cbegin()+check_point_row_cols,
                                                m_inner.cend(),
                                                [&k](auto elem){return elem>k;});
                    col_r_idx.push_back(std::distance(m_inner.cbegin(),col_ptr-1));
                    check_point_row_cols+=col_r_idx.back();
                }
                k++;
            }
        }
        else                                                            //uncompressed
        {
            //map to temporary store row idx-th
            MatrixUncompressed<T,O> temp_map;
            //copy only the element of the row idx-th
            std::ranges::copy_if(m_mat_uncomp,std::inserter(temp_map,temp_map.begin()),[&idx](auto p){return p.first[0]==idx;});
            
            //number of elements in the row idx-th
            nr = temp_map.size();

            //reserving memory
            col_r_idx.reserve(nr);
            val_r_idx.reserve(nr);

            //pushing back the cols and the values in the row idx-th
            std::ranges::for_each(temp_map,[&col_r_idx,&val_r_idx](auto p){col_r_idx.push_back(p.first[1]); val_r_idx.emplace_back(p.second);});

            //clearing memory
            temp_map.clear();      
        }
    }
    return std::make_pair(col_r_idx,val_r_idx); 
}


/*!
* Getting col idx-th.
*/
template<class T, StorageOrder O>
std::pair<std::vector<std::size_t>,std::vector<T>>
algebra::Matrix<T,O>::get_col(const std::size_t &idx) const{

    //not checking the presence: it will not store the index: ok
    //functions that uses it will check for col presence

    //storing row indeces of the col
    std::vector<std::size_t> row_c_idx;         
    //storing values of the col
    std::vector<T> val_c_idx;
    //number of elements in the row
    std::size_t nc;

    if constexpr (O == StorageOrder::ColumnWise)                           //ColWise
    {
        if (this->is_compressed())                                          //compressed
        {   
            //iterators to take the range of the col in the 
            //compressed representation, for row indeces and values
            auto start_col = m_outer.cbegin() + m_inner[idx];           
            auto sentinel_col = m_outer.cbegin() + m_inner[idx+1];
            auto start_val = m_val_comp.cbegin() + m_inner[idx];
            auto sentinel_val = m_val_comp.cbegin() + m_inner[idx+1];

            nc = std::distance(start_col,sentinel_col);

            //memory reserved
            row_c_idx.reserve(nc);
            val_c_idx.reserve(nc);

            //copy of the value: done in parallel since operating on std::vector
            std::copy(std::execution::par,
                      start_col,sentinel_col,std::back_inserter(row_c_idx));
            std::copy(std::execution::par,
                      start_val,sentinel_val,std::back_inserter(val_c_idx));
        }   
        else                                                            //uncompressed
        {   
            //problem if looking for the first col, since upper bound retrieve the first greater
            //and not greater_equal: cannot give {0,-1} as key, since these elements are size_t, it is necessary to
            //take account for the element (0,0)
            if (idx==0)
            { 
                key_type ub{0,idx+1};                                   //upper bound of the row
                auto start_r = m_mat_uncomp.cbegin();                   //beginning of the map
                auto sentinel_r = m_mat_uncomp.lower_bound(ub);         //sentinel for the row

                nc = std::distance(start_r,sentinel_r);                 //number of the elements in the row
                row_c_idx.reserve(nc);                                  //reserving memory
                val_c_idx.reserve(nc);

                //putting elements in the right vector
                std::for_each(start_r,sentinel_r,[&row_c_idx,&val_c_idx](auto p){row_c_idx.push_back(p.first[0]); val_c_idx.emplace_back(p.second);});
            }
            else
            {
                key_type lb{this->m()+1,idx-1};                         //lower bound of the row
                key_type ub{0,idx+1};                                   //upper bound of the row
                auto start_r = m_mat_uncomp.upper_bound(lb);
                auto sentinel_r = m_mat_uncomp.lower_bound(ub);         //sentinel for the row

                nc = std::distance(start_r,sentinel_r);                 //number of the elements in the row
                row_c_idx.reserve(nc);                                  //reserving memory for the storaging
                val_c_idx.reserve(nc);

                //putting elements in the right vector
                std::for_each(start_r,sentinel_r,[&row_c_idx,&val_c_idx](auto p){row_c_idx.push_back(p.first[0]); val_c_idx.emplace_back(p.second);});
            }
        }
    }   
    else                                                                //RowWise
    {   

        if (this->is_compressed())                                      //compressed
        {   
            //number of element in the col: can be done in parallel since it is std::vector
            nc = std::count(std::execution::par,                        
                            m_outer.cbegin(),m_outer.cend(),
                            idx);

            //memory reserving
            row_c_idx.reserve(nc);
            val_c_idx.reserve(nc);

            //variable for looping
            std::size_t k{0};
            //idea: -the value can be retrieved triavially since it is only matter of 
            //       accessing the position of m_value_comp in which m_outer has idx as element
            //      -the indeces of the rows can be retrieved as:
            //       the j-th elem of m_inner contains the number of elem up to the row (j-1)-th
            //       so, looping over all the rows, starting from the last one in which an elem is
            //       found, until the corrispettive elem in m_inner is bigger of the position
            //       of the elem in m_outer, that indicate which elems are in the idx-th col
            //Looping until all the elemnts have been found
            std::size_t check_point_col_rows{0};

            while (val_c_idx.size()<nc)
            {
                if (m_outer[k]==idx)
                {   
                    //pushing the value
                    val_c_idx.emplace_back(m_val_comp[k]);

                    //searching the first row available as
                    auto row_ptr = std::find_if(m_inner.cbegin()+check_point_col_rows,
                                                m_inner.cend(),
                                                [&k](auto elem){return elem>k;});
                    row_c_idx.push_back(std::distance(m_inner.cbegin(),row_ptr-1));
                    check_point_col_rows+=row_c_idx.back();
                }
                k++;
            }
        }
        else                                                            //uncompressed
        {
            //map to temporary store col idx-th
            MatrixUncompressed<T,O> temp_map;
            //copy only the element of the col idx-th
            std::ranges::copy_if(m_mat_uncomp,std::inserter(temp_map,temp_map.begin()),[&idx](auto p){return p.first[1]==idx;});
            
            //number of elements in the col idx-th
            nc = temp_map.size();

            //reserving memory
            row_c_idx.reserve(nc);
            val_c_idx.reserve(nc);

            //pushing back the rows and the values in the col idx-th
            std::ranges::for_each(temp_map,[&row_c_idx,&val_c_idx](auto p){row_c_idx.push_back(p.first[0]); val_c_idx.emplace_back(p.second);});

            //clearing memory
            temp_map.clear();      
        }
    }
    return std::make_pair(row_c_idx,val_c_idx); 
}