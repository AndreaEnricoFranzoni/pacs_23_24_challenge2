# pacs_23_24_challenge2

The repository contains a class to handle sparse matrices, from reading them from a MMF(Matrix Market Format) file to efficiently updating their values dynamically and efficiently performing matrix-vector product, matrix-matrix product and norm one, infinity and Frobenius due to its capability to compress/decompress the format of the matrix.

Matrices store elements of the template parameter: everything works for integral, floating and complex types.
Matrices can be stored row-wise or column-wise: it is a template parameter, and once a matrix is constructed with a storage ordering, it is not possible to change it anymore.
Uncompressed format: COOmap format is used. Key is a std::array<std::size_t,2>. Handling the different ways of storage ordering is done using std::conditional: for column-wise storing, the lexicographical ordering has to start from the second indices of the key, and its defined by a functor passed if column-wise is requested.
Compressed format: CSR/CSC format is exploited, using std::vector for all the three vectors.

Requirements: the code relays on std::execution::par whenever possible, since using std::vector that are sequential containers. It is necessary to link libttb library.
			 The pacs utility “Chrono.hpp” is required, so the link with libpacs.so is required.


Running the code:
The file Makefile has to be modified putting the proper PACS_ROOT, and the beginning.
make: compiling
./main: executing main function

Documentation:
In the file Doxyfile, is necessary to modify in the line 2239: after  INCLUDE_PATH = ./include \ .\  ,is necessary to put PACS_ROOT/include, with PACS_ROOT being the proper pacs root
make doc. will give back a directory doc: in the subfolder html, opening index.html will give back the Doxygen documentation for the code

make disctlean remove object files and executables, but not the doc folder: once is created, it has be removed manually.


Main function: firstly, the matrix contained in lnsp_131.mtx is read through the MMF reader. It is stored yet as RowWise yet as ColumnWise matrix. Then, yet a std::vector yet an algebra::Matrix<T,O> with 1-col, with the right dimensions, are constructed, filled with increasing values starting from 0 up to the number of cols of the matrix read.
The main tests the performances using “Chrono.hpp” doing :
-matrix-vector product with std::vector 
-matrix-vector product with algebra::Matrix<T,O>
-matrix-matrix product doing the square of the read matrix
-evaluating the norm One of the read matrix
-evaluating the norm Infinity of the read matrix
-evaluating the Frobenius norm of the read matrix

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
