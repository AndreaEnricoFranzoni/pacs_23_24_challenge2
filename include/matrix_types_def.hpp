#ifndef MATRIX_TYPES_CHALLENGE2_PACS_HPP
#define MATRIX_TYPES_CHALLENGE2_PACS_HPP

#include <vector>
#include <array>
#include <map>
#include <type_traits>
#include <concepts>

/*!
*Types and structures definitions to build sparse matrices.
*/


/*!
* To specify the matrix storage ordering.
*/
enum StorageOrder{
    RowWise = 0,
    ColumnWise = 1
};


/*!
* To specify the norm evaluation types.
*/
enum NormType{
    One = 0,
    Infinity = 1,
    Frobenius = 2
};


/*!
* Doing tag dispatching to relays on the correct function to evaluate the norm of the matrix.
* @param N: template parameter for the norm type to be mapped by std::integral_constant
*/
template <NormType N>
using NT = std::integral_constant<NormType,N>;


/*!
* Definition of the key for the std::map to contain uncompressed matrix: array of size 2 with std::size_t elements.
*/
using key_type = std::array<std::size_t,2>;


/*!
* Functor defining column-wise ordering to be passed as template parameter for columnwise storaging.
* @param lhs: const lvalue ref to an element of type key_type: left side of the comparison
* @param rhs: const lvalue ref to an element of type key_type: right side of the comparison
* @return true if the left side is less than the right one checking first the number of col, and, if tied, the number of the row, false otherwise
*/
struct CompOp_ColWise{
    bool operator() (const key_type &lhs, const key_type &rhs) const{
        return ((lhs[1]<rhs[1]) or ((lhs[1]==rhs[1]) and (lhs[0]<rhs[0]))); 
    };
};


/*!
* Template definition for sparse uncompressed matrix.
* @param T: the type of elements that are stored
* @param O: StorageOrder 
* @return using std::conditional: O can be interpreted as false if is 0 (RowWise): the map will have the defaulted lexicographical less operator;
*                                 O can be interpreted as true if is 1 (ColumnWise): the map will have the CompOp_ColWise functor as less operator
*/
template <class T,StorageOrder O>
using MatrixUncompressed = std::conditional<O,
                                            std::map<key_type, T, CompOp_ColWise>,
                                            std::map<key_type,T>>::type;


#endif ///MATRIX_TYPES_CHALLENGE2_PACS_HPP