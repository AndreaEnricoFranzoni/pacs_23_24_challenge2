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
-`main.cpp`: matrix in `lnsp_131.mtx` is read through MMF reader. It is stored both row-wise and column-wise. To test matrix-vector product, both an `std::vector<T>' and an `algebra::Matrix<T,O>` of the right dimensions are constructed, filled with increasing values starting from 0.
The main function, using `Chrono.hpp' test the perdormances of:
1. matrix-vector product with the vector being `std::vector`;
2. matrix-vector product with the vector being a 1-col `algebra::Matrix<T,O>`;
3. matrix-matrix product doing the square of the read matrix;
4. norm One of the read matrix;
5. norm Infinity of the read matrix;
6. Frobenius norm of the read matrix.
Performances are evaluated for both row-wise and column-wise storage, and in both cases for uncompressed and compressed format: elapsed time is lower in compressed format.
Switching
The performances are evaluated for the following format of the read matrix: row-wise uncompressed, row-wise compressed, col-wise uncompressed, col-wise compressed. The time elapsed is lower when relaying on compressed format, for both the storage orders.
Switching off other processes and activating optimization allows the program to work at its best. 

Then, a small overview on matrix with complex coefficients is done: a 2x2 matrix is constructed with complex coefficients, as well as an std::vector with complex coefficients of size 2 and another algebra::Matrix 2x1 with complex coefficients.
The same operations as before are committed, but instead than evaluating the performance in this case simply results are displayed in order to check that results and types are coherent.
For the sake of simplicity, only the case of row-wise uncompressed is displayed, while the other are left commented. But obviously is possible to see the correct results by simply leaving out the comments.

In the folder /include:
matrix.hpp contains the declaration of the template class, under the namespace algebra, as well as the definition of the friend operator * and the friend operator <<.
matrix_imp.hpp contains the definitions of compress, uncompress methods and the call operator, in its const and non-const version.
matrix_types_def.hpp contains the definition of the enumerator and types used, the functor to handle col-wise ordering, as well as the use of std::conditional to use the correct less operator.
					In general, if constexpr is used to compile only the code lines related to a specific matrix, depending on its storing.
					For the norm, tag dispatching is exploited.
matrix_reader_imp.hpp contains the definition of the reader of MMF format.
matrix_norm_imp.hpp contains the definitions of the three different functions to evaluate the norm.
matrix_get_row_col_imp.hpp contains the definitions of the method for checking presence of a row, a col and and element, as well as the ones to get a row or a col.

In general, the code exploits STL algorithms, yet in their original version yet in their constrained version.
In particular, for getting a row or a col, the code tries to get iterators to the first and after-the-last one element of the row/col for the uncompressed format, while for the compressed one relays on understanding from outer and inner indices if the row/col is present and which are its extrema. 
For uncompressed format, the overload of the less operand makes trivial both cases.
While, for compressed format, is easy to retrain a row for row-wise and a col for col-wise storage, while the viceversa is not so trivial: a better explanation is done on the comment of the code.

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
