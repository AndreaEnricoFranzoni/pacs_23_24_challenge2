#ifndef HH_SPARSE_MATRIX_CHALLENGE2_HH
#define HH_SPARSE_MATRIX_CHALLENGE2_HH

#include <iostream>
#include <cassert>
#include <complex>
#include <algorithm>
#include <ranges>
#include <iterator>
#include <numeric>
#include <execution>
#include "matrix_types_def.hpp"


namespace algebra
{

/*!
*Template class for sparse matrix.
* @param T is the type stored: works with integral, floating and complex types
* @param O is the way of storing elements in the matrix, either RowWise (default) either ColumnWise
*/

template <class T, StorageOrder O = StorageOrder::RowWise>      //default: row-major
class Matrix{

public:

   /*!
   * Constructor may take number of rows, columns and format.
   * @param m number of rows
   * @param n number of cols
   * @param compressed if in compressed format (default is false, so can be filled)
   */
    Matrix(size_t m=0, size_t n=0, bool compressed = false) 
        :   m_m(m), m_n(n), m_nnz(0), m_compressed(compressed){}


    /*!
    * Getter for number of rows. 
    * @return the private member m_m: read-only
    */
    auto m() const {return m_m;};
    /*!
    * Setter for number of rows. 
    * @return the private member m_m, modifying it
    */
    auto &m() {return m_m;};


    /*!
    * Getter for number of cols.
    * @return the private member m_n: read-only
    */ 
    auto n() const {return m_n;};
    /*!
    * Setter for number of cols. 
    * @return the private member m_n, modifying it
    */
    auto &n() {return m_n;};

    
    /*!
    * Getter for number of non-zero elements.
    * @return the private member m_nnz, the number of non-zero elemts: only const version since cannot be modified if not adding or removing an element
    */
    auto nnz() const {return m_nnz;};

    
    /*!
    * Getter for the matrix in uncompressed state.
    * @return the private member m_mat_uncomp, that is the matrix in uncompressed state, so a std::map<key_type,StorageOrder::RowWise> or a std::map<key_type,StorageOrder::ColumnWise,CompOp_ColWise>;
    *         only const version since cannot be set
    */
    auto matUnc() const {return m_mat_uncomp;};

    
    /*!
    * Getter for the values of the matrix in compressed state.
    * @return the private member m_val_comp, that is the values of the matrix in uncompressed state, so a std::vector<T>;
    *         only const version since cannot be set
    */
    auto matCom() const {return m_val_comp;};

    
    /*!
    * Getter for the inner indeces.
    * @return the private member m_inner, that is a std::vector<std::size_t>; only const version allowed
    */
    auto innerCom() const {return m_inner;};

    
    /*!
    * Getter for the outer indices.
    * @return the private member m_outer, that is a std::vector<std::size_t>; only const version allowed
    */
    auto outerCom() const {return m_outer;};

    
    /*!
    * Getter for if it is compressed.
    * @return the private member m_compressed, that is a bool: read-only
    */
    bool is_compressed() const {return m_compressed;};
    /*!
    * Setter for if it is compressed. 
    * @return the private member m_compressed, modifying it
    */
    bool &is_compressed() {return m_compressed;};


    /*!
    *   Clearing the buffers used to store elements. Called to resize it, to compress and to uncompress the matrix.
    */
    inline void clear_buffer(){
        if (this->is_compressed())
        {
            this->m_mat_uncomp.clear();
        }
        else
        {
            this->m_val_comp.clear();
            this->m_outer.clear();
            this->m_inner.clear();
        }
    }


    /*!
    * Resizing the matrix, after having cleared the buffer(s) 
    * @param m_new: new number of rows
    * @param n_new: new number of cols
    */
    inline void resize(const size_t m_new, const size_t n_new){
        this->clear_buffer();
        this->m() = m_new;
        this->n() = n_new;
    }


    /*!
    * Checking if row idx-th is present.
    * @param idx: the index of the row whose presence has to be checked
    * @return true if there is at least a nnz element in that row, false otherwise;
    *         aborting if the param exceeds the dimension of the matrix
    */
    bool check_presence_row(const std::size_t &idx) const;


    /*!
    * Checking if col idx-th is present.
    * @param idx: the index of the col whose presence has to be checked
    * @return true if there is at least a nnz element in that col, false otherwise;
    *         aborting if the param exceeds the dimension of the matrix
    */
    bool check_presence_col(const std::size_t &idx) const;

    
    /*!
    * Checking if element(i,j) is present (non-null).
    * @param i: row's index of the element whose presence has to be checked
    * @param j: col's index of the element whose presence has to be checked
    * @return true if elem(i,j) is non-null, false otherwise;
    *         aborting if the params exceed the dimensions of the matrix
    */
    bool check_presence(const std::size_t &i, const std::size_t &j) const;

    
    /*!
    * Extracting row idx-th.
    * @param idx: index of the row that has to be extracted
    * @return a pair: first:  std::vector<std::size_t> : the indeces of the cols of row idx-th that contain nnz elements
    *                 second: std::vector<T>           : the values of the nnz elements in row idx-th
    */
    std::pair<std::vector<std::size_t>,std::vector<T>> get_row(const std::size_t &idx) const;

    
    /*!
    * Extracting col idx-th.
    * @param idx: index of the col that has to be extracted
    * @return a pair: first:  std::vector<std::size_t> : the indeces of the rows of col idx-th that contain nnz elements
    *                 second: std::vector<T>           : the values of the nnz elements in col idx-th
    */    
    std::pair<std::vector<std::size_t>,std::vector<T>> get_col(const std::size_t &idx) const;

    
    /*!
    * Compressing the matrix. Passing from uncompress to compress state, clearing the container for the uncompress state.
    */
    void compress();        


    /*!
    * Uncompressing the matrix. Passing from compress to uncompress state, clearing the containers for the compress state.
    */
    void uncompress();

    
    /*!
    * Reading the matrix from a .mtx/.txt file in which is stored as MatrixMarketFormat.
    * Assuming that the comments are at the beginning of the file, and that the file follows MMF standards.
    * The reading is done in uncompress state, and so the matrix is returned as uncompressed.
    * @param file_name: name of the file in which the matrix is written in MMF
    * @return true if the matrix has been read correctly, false otherwise
    */
    bool reader_mmf(const std::string &file_name);  


    /*!
    * Get element A(i,j). Const version, read-only.
    * @param i: number of row
    * @param j: number of col
    * @return the value in position (i,j) if non-zero, 0 casted to type T otherwise
    */
    T operator()(std::size_t i, std::size_t j) const;


    /*!
    * Upload element A(i,j) with a new value.
    * Adding a zero means removing an element.
    * Can update a null element only if uncompressed.
    * Can remove an element only if uncompressed.
    * @param i; number of row
    * @param j: number of col
    * @param newvalue: value to be inserted
    * @return if inserted, return the new value inserted. If not, raises an error since there is the attempt to modify a null element in compressed format
    */
    T & operator()(std::size_t i, std::size_t j, const T newvalue);


    /*!
    * Matrix-vector product. Matrix and vector are read only, so are passed as const lvalue ref.
    * Works with complex elements. 
    * Aborting if the dimensions are wrong
    * @param M: algebra::Matrix<T,O> of dimensions m x n
    * @param v: std::vector<T> of size n
    * @return a std::vector<T> of size m containg the matrix-vector product
    * 
    * @note exploitation of random access of std::vector once indeces of nnz elements of a row/col are known
    */
    friend
    std::vector<T>
    operator * (const Matrix<T,O> &M,const std::vector<T> &v){
        
        //checking if the size of the std::vector is correct
        assert(v.size()==M.n());

        //storing results
        std::vector<T> results;
        
        //ROWWISE: for each row: computing the scalar product only with the element
        //                        of v that corresponds at non null element in each row

        if constexpr(O==StorageOrder::RowWise)                              //ROWWISE
        {   
            results.reserve(M.m());                                         //reserving memory for results

            for (size_t i = 0; i < M.m(); ++i)                              //Looping on all the rows
            {
                T temp = static_cast<T>(0);                                 //initializing the scalar product between row i-th and v
                auto row_i = M.get_row(i);                                  //getting row i-th                      

                for (size_t j = 0; j < row_i.first.size(); ++j)             //loop on all the cols of row i-th that contain non null elements
                {                                                           
                    temp += row_i.second[j]*v[row_i.first[j]];              //scalar product: accessing only the elements of v that correspond to non null element in the row i-th 
                }
                results.emplace_back(temp);                                  //storing the result
            }
        }
        //COLWISE: for each i in {0,...,n-1}, for column i-th: computing the product of each non null element of the 
        //         col j-th with the element j-th of the vector. Then, summing up these vectors 

        else                                                                //COLWISE
        {   
            results.resize(M.m());                                          
            std::ranges::fill(results,static_cast<T>(0));                   //initializing results

            for (size_t i = 0; i < M.n(); ++i)                              //looping on the cols
            {   
                auto col_i = M.get_col(i);                                  //getting the col i-th
    
                for (size_t j = 0; j < col_i.first.size(); ++j)             //looping on all the rows that contain non zero elements
                {
                    results[col_i.first[j]] += col_i.second[j]*v[i];        //summing up the products of the col i-th with the element i-th of v
                }                
            }
        }
        return results;
    }


    /*!
    * Matrix-matrix product, being able to handle also a matrix-vector product if the dimensions are correct
    * The storaging of the result will be the same of the first factor.
    * Can know only on the fly if is doing a matrix-matrix product or actually a matrix-vector one: it will compile all the portion of the code.
    * Aborting if the dimensions are incorrect.
    * @param M: algebra::Matrix<T,O> of dimensions m1 x n1
    * @param v: algebra::Matrix<T,O1> of dimensions m2 x n2
    * @return an algebra::Matrix<T,O> of dimensions m1 x n2 as the result of the product between the two matrices
    * 
    * @note the scalar product exploits the constant call operator to relay on the same code structure of the std::vector version
    */
    template<StorageOrder O1>
    friend
    Matrix<T,O>
    operator * (const Matrix<T,O> &M, const Matrix<T,O1> &v){

        //checking dimensions
        assert(M.n()==v.m());

        //assessing if it is a product with a vector
        bool vector_prod;
        (v.n()==1) ? (vector_prod=true) : (vector_prod=false);

        //resulting matrix: same storage of M, uncompressed since it is necessary to fill it
        Matrix<T,O> results(M.m(),v.n());

        if (vector_prod==1){                                                    //MAT-VEC PROD
                
            auto rhs = v.get_col(0);                                            //getting the vector as a col
            
            if constexpr(O==StorageOrder::RowWise)                              //ROWWISE
            {   
                for (size_t i = 0; i < M.m(); ++i)                              //Looping on all the rows
                {
                    T temp = static_cast<T>(0);                                 //initializing the scalar product between row i-th and v
                    auto row_i = M.get_row(i);                                  //getting row i-th                      

                    for (size_t j = 0; j < row_i.first.size(); ++j)             //looping on all the cols of row i-th that contain non null elements
                    {                                                           
                        temp += row_i.second[j]*v(row_i.first[j],0);            //scalar product: I access using the const call operator only the elements of v that correspond to non null element in the row i-th 
                    }
                    results(i,0,temp);                                          //storing the result
                }
            }
            //COLWISE: for each i in {0,...,n-1}, for column i-th: computing the product of each non null element of the 
            //         col j-th with the element j-th of the vector. Then, summing up these vectors 

            else                                                                //COLWISE
            {   
                std::vector<T> temp;                                            //temporary std::vector to store results (since in this way is possible to updating all the rows), that will be copied in the Matrix<T,O>
                temp.resize(M.m());                                          
                std::ranges::fill(temp,static_cast<T>(0));                      //initializing results

                for (size_t i = 0; i < M.n(); ++i)                              //looping on the cols
                {   
                    auto col_i = M.get_col(i);                                  //get the col i-th
    
                    for (size_t j = 0; j < col_i.first.size(); ++j)             //looping on all the rows that contain non zero elements
                    {
                        temp[col_i.first[j]] += col_i.second[j]*v(i,0);         //summing up the products of the col i-th with the element i-th of v: accessing using the const call operator of Matrix<T,O>
                    }                
                }

                for (auto it = temp.begin(); it != temp.end(); ++it)
                {
                    results(std::distance(temp.begin(),it),0,*it);              //filling the algebra::Matrix<T,O> results
                }
                
                temp.clear();                                                   //clearing the temporary container
            }

        }
        else                                                                    //MAT-MAT PROD
        {
            //for every pair (i,j), checking if for row_i and col_j there is at least one nnz elem: if yes, the inner product of row and col is done:
            //row_i and col_j are taken using the getters functions, that would take care of storaging and compressed/uncompressed
            for (size_t j = 0; j < v.n(); ++j)
            {
                auto col_j = v.get_col(j);                                                                   //getting the col

                for (size_t i = 0; i < M.m(); ++i)
                {
                    auto row_i = M.get_row(i);                                                               //getting the row
                    
                    //how many elements in row_i and col_j are nnz and in the same position
                    auto n_idx = std::ranges::count_if(row_i.first,[&col_j](auto el){return(std::find(std::execution::par,col_j.first.cbegin(),col_j.first.cend(),el)!=col_j.first.cend());});
                    
                    if (n_idx!=0)                                                                            //if the resulting elem (i,j) is nnz
                    {
                        std::vector<std::size_t> idx;                                                        //storing the indeces of corrispective elements in row_i and col_j non simultaneously null              
                        idx.reserve(n_idx);                 

                        std::ranges::set_intersection(row_i.first,col_j.first,std::back_inserter(idx));      //taking the indeces of the nnz elems

                        auto res=std::transform_reduce(std::execution::par,                                  //inner prod of row_i with col_j: in parallel
                                                       idx.begin(),                                          //works also with complex values
                                                       idx.end(),
                                                       static_cast<T>(0),
                                                       std::plus{},
                                                       [&row_i,&col_j](auto el)
                                                       {    
                                                            //taking the right value in row and col to be multiplied and then aggregated with the sum
                                                            auto idx_r = std::distance(row_i.first.begin(),std::ranges::find(row_i.first,el));
                                                            auto idx_c = std::distance(col_j.first.begin(),std::ranges::find(col_j.first,el));

                                                            return row_i.second[idx_r]*col_j.second[idx_c];
                                                        });

                        results(i,j,res);                                                                    //storing the result in the resulting matrix
                        idx.clear();                                                                         //clearing the temporary container
                    }                    
                }
            }
        }
        return results;     
    }
    
    
    /*!
    * Stream operator to visualize an algebra::Matrix<T,O>.
    * If uncompressed: store the matrix in the stream object in order to visualize it as a classical table, filled with zeros.
    * If compressed: store the matrix in the stream object in order to visualize it as three arrays according to CSR/CSC format.
    * Good to visualize small matrices.
    * @param str: an lvalue ref to an std::ostream object
    * @param M: the algebra::Matrix<T,O> to be visualized after being store in the stream object
    * @return the ref to @param str with the matrix stored
    * 
    * @note taking out the comment bars from that lines of code allows to have the print of storage order and compressed/uncompressed. For the sake of semplicity, they are left out.
    */
    friend 
    std::ostream &
    operator << (std::ostream & str, Matrix<T,O> const &M)
    {   
        //(O==StorageOrder::RowWise) ? (str << "Row-wise storage" << '\n') : (str << "Col-wise storage" << '\n'); 
        
        if (!M.is_compressed())
        {   
            //str << "Matric uncompressed" << '\n';
            for (size_t i = 0; i < M.m_m; ++i)
            {
                for (size_t j = 0; j < M.m_n; ++j)
                {
                    str << M(i,j) << ' ';
                }
                str << '\n'; 
            }
        }   
        else
        {   
            //str << "Matric compressed" << '\n';
            str << "Values:" << '\n';
            for (size_t i = 0; i < M.m_nnz; ++i)
            {
                str << M.m_val_comp[i] << '\n';
            }
            str << "Inner" << '\n';
            for (size_t i = 0; i < M.m_inner.size(); ++i)
            {
                str << M.m_inner[i] << '\n';
            }
            str << "Outer:" << '\n';
            for (size_t i = 0; i < M.m_nnz; ++i)
            {
                str << M.m_outer[i] << '\n';
            }   
        }
        return str;  
    }
    

    /*!
    * Norm evaluation of the matrix.
    * Template function.
    * @param N: template parameter that indicates an enum for the requested norm
    * @return tag dispatching to the correct function to evaluate the norm as requested
    */
    template<NormType N>
    auto norm() const{ return norm(NT<N>{});};  
    

private:
    /*!Number of rows.*/
    size_t m_m;
    /*!Number of cols.*/                              
    size_t m_n;
    /*!Number of non-zero elements*/                              
    size_t m_nnz;
    /*!True if the matrix is compressed, false if uncompressed*/                            
    bool m_compressed;                       

    //Unompressed format containers
    /*!Uncompressed container. 
    *std::map<key_type,T> that has the comparison operator overloaded if storage is ColumnWise.
    *Size is nnz.*/
    MatrixUncompressed<T,O> m_mat_uncomp;    

    //Compressed format containers
    /*!std::vector<T> that contains the values of non-zero elements.
    *Size is nnz.*/
    std::vector<T> m_val_comp;    
    /*!std::vector<T> that contains inner indices.
    *If RowWise: size is m+1, elem i-th of the vector contains the number of nnz elements up to row i-th.
    *If ColumnWise: size is nnz, and contains, for each value of m_val_comp, their column.*/           
    std::vector<std::size_t> m_inner;        
    /*!std::vector<T> that contains outer indices.
    *If RowWise: size is nnz, and contains, for each value of m_val_comp, their row.
    *If ColumnWise: size is n+1, elem i-th of the vector contains the number of nnz elements up to col i-th.*/  
    std::vector<std::size_t> m_outer;        

    /*!
    *Evaluation of Norm One of the matrix.
    * @return the value of the Norm One
    */
    auto norm(NT<NormType::One>) const;
    /*!
    *Evaluation of Norm Infinity of the matrix.
    * @return the value of the Norm Infinity
    */      
    auto norm(NT<NormType::Infinity>) const;
    /*!
    *Evaluation of Frobenius Norm of the matrix.
    * @return the value of the Frobenius Norm
    */
    auto norm(NT<NormType::Frobenius>) const;//Norm Frobenius
};
}

#include "matrix_imp.hpp"
#include "matrix_get_row_col_imp.hpp"
#include "matrix_reader_imp.hpp"
#include "matrix_norm_imp.hpp"


#endif     //SPARSE_MATRIX_CHALLENGE2_HPP