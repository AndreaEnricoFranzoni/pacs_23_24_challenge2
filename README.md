# Content
The repository contains a class to handle sparse matrices, able to read them from a MMF(Matrix Market Format) file, to efficiently update their values dynamically and to efficiently perform matrix-vector product, matrix-matrix product, norm one, norm infinity and Frobenius norm due to its capability to compress/uncompress the format of the matrix.

# Data structure
The class is, under the namespace `algebra`, a template class: `Matrix<T,O>`.

The template parameter `T` refers to the type of the values stored: works for integral, floating and complex types, but the operations are suited only for scalar values.

The template parameter `O` referes to the storage order: `StorageOrder::RowWise` and `StorageOrder::ColumnWise`. Once a matrix has been constructed with a storage order, it is not possible to change it anymore.

Matrices can switch between uncompressed and compressed format.

Uncompressed format: COOmap. Data are stored using `std::map`, where key is `std::array<std::size_t,2>`, the first element of the array referring to the row where an element is stored, the second one to the row, and value is `T`. Lexicographical ordering of the `std::array` handles row-wise storage, while a functor to handle the column-wise storage is passed in the other case, using `std::conditional`.

Compressed format: CSR/CSC. Data are stored using three `std::vector`, storing values of type `T` for the vector concerning values storage, while of type `std::size_t` for the other two that handle indeces storaging.

# Requirements 
Code exploits on `std::execution::par` whenever it is possible. Linking `libttb` library is necessary.

The PACs utility `Chrono.hpp` is used. Linking `libpacs.so` is also required.

The file `Makefile` has to be modified, putting the proper `PACS_ROOT` at the beginning.

The file `Doxyfile` has to be modified at the line 2239: after `INCLUDE_PATH = ./include \ .\ `, is necessary to put `PACS_ROOT/include`.

# Running the code

Compiling, with the flag `-O3` for optimization already set:
~~~
make
~~~

Running the executable:
~~~
./main
~~~

Documentation: 
~~~
make doc
~~~
create a directory `doc`, in which there are two subfolders (`html` and `latex`): in `html`, opening the file `index.html`, will open on the browser the code documentation constructed with Doxygen.

Cleaning from object files and executables:
~~~
make distclean
~~~

Removing documentation:
~~~
rm -r doc
~~~

# Files content
-**`main.cpp`**: matrix in `lnsp_131.mtx` is read through MMF reader. It is stored both row-wise and column-wise. To test matrix-vector product, both an `std::vector<T>' and an `algebra::Matrix<T,O>` of the right dimensions are constructed, filled with increasing values starting from 0.
The main function, using `Chrono.hpp' test the perdormances of:
1. matrix-vector product with the vector being `std::vector`;
2. matrix-vector product with the vector being a 1-col `algebra::Matrix<T,O>`;
3. matrix-matrix product doing the square of the read matrix;
4. norm One of the read matrix;
5. norm Infinity of the read matrix;
6. Frobenius norm of the read matrix.

Performances are evaluated for both row-wise and column-wise storage, and in both cases for uncompressed and compressed format: elapsed time (microseconds) is lower in compressed format.
Switching off other processes allows the parallelization due to `std::execution::par` to work better.

7. Construction of a 2x2 matrix with complex coefficients and of a vector with complex coefficients as before. Matrix-vector product, in both its version, matrix-matrix product (square of the matrix) and all the three types of norm's evaluation are tested, displaying results to show that returned types are coherent. The test is done in all the storage order and formats, but only the row-wise uncompressed, for the sake of semplicity, is shown, while other cases are left commented.

-**folder `/include`**:
- `matrix.hpp`: declaration of the class under the nanmespace `algebra`, as well as the definitions of the friend `operator*` and the friend `operator <<`;
- `matrix_imp.hpp`: definition of `compress` and `uncompress` methods, as well as the ones of the call `operator()`, in both its const and non-const version;
- `matrix_types_def.hpp`: definition of the enumerators and types used, of the functor handling column-wise ordering;
matrix_types_def.hpp contains the definition of the enumerator, types and the functor to handle col-wise ordering;
- `matrix_reader_imp.hpp`: definition of the reader of MMF;
- `matrix_norm_imp.hpp`: definition of the three functions to evaluate the norm;
- `matrix_get_row_col_imp.hpp`: defintion of `check_presence_row`, `check_presence_col`, `check_presence`, `get_row` and `get_col` methods.
					
# Features
- Use of `std::conditional` to construct correctly the uncompressed matrix depending on its storaging.
- Use of `if constexpr(O==StrorageOrder::RowWise)` and `else` to compile only the code regarding the specific storage.
- Use of tag dispatching for evaluating the norm as requested.
- Expolitation of `STL` algorithms, both in their original (to relay on parallelization if possible) and constrained version. Getting a row or a col is done looking for the correct iterator describing that range of values if uncompressed, while looking for the correct indices of the stored values if compressed. If uncompressed, is trivial due to the overloading of the `less` operator. If compressed, is easy to retrain a row/col if row-wise/column-wise respectively, while the viceversa is not trivial (since the values will not be stored contiguously a priori). Further explanations in the comments of the code. 




Issues:
-friend operators are defined in the header file. Putting them in a .cpp did not work since it would not compile in case;
-resizing the matrix means clearing it;
-the operator*, if done with an std::vector, return a std::vector;
-the operator*, if done with an algebra::Matrix<T,O> has some problems:
			-has to check on the fly if it is a matrix-vector or a matrix-matrix product: I could not differentiate it since the arguments would have been the same: I am checking with an if every 				  time, but I have no idea about how to improve this, since I am compiling anyway all the code in that operation;
			-the matrix returned will have the same storage of the first factor: also in this case, I do not know how to pass this type of parameter, such as saying: these are my factors, but I want 				
			 it with a specific storage ordering
-in get_col and get_row, for uncompressed format, it is necessary to check manually if we are not extracting the row/col 0: working with upper_bound and lower_bound, is necessary to handle this case separately. But, also in this case, will compile both branches of the if.
-the non-const call operator returns a reference to a casted 0 if in a compress state matrix I try to add/remove an element. But returning 0 here is compulsory: I know that is an error doing it this way, but I do not know how to handle it differently.
