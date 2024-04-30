#include "matrix.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <exception>
#include <type_traits>
#include <iomanip>

/*!
* Definition of reader_mmf, that allows to read a matrix from MatrixMarketFormat.
* 
* @note assuming that the file will have the standards of a MMF: comments, line with dimension and then all the stored elements. Everyhing will be read in this order.
*/
template<class T,StorageOrder O>
bool
algebra::Matrix<T,O>::reader_mmf(const std::string &file_name)
{
    std::ifstream file(file_name);

    //checking if is possible to open the file
    if(!file)
        throw std::runtime_error(std::string("Error in opening file ") + file_name);

    //Ignoring comments headers, shifting toward the beginning of the matrix definition
    while (file.peek() == '%') file.ignore(2048, '\n');

    size_t count = 0;
    bool sizes_read = false;
    int num_row{0}, num_col{0}, num_lines{0};

    //storing the lines of the document
    std::string buffer;

    //looping over the lines of the document
    while (std::getline(file,buffer))
    {   
        //if sizes not already read
        if(!sizes_read){

            //getting the line in which there are the dimensions
            std::istringstream line(buffer);
            line >> num_row >> num_col >> num_lines;

            //checking that dimensions are ok
            if (num_row<=0 or num_col<=0 or num_lines<0)
                throw std::runtime_error("Not possible to have negative or null dimensions or a negative number of nnz elements");
            
            //resizing the matrix and making it uncompressed in order to be able to store elements
            this->resize(num_row,num_col);
            this->is_compressed() = false;

            sizes_read = true;
        }
        //if is reading the actual elements
        else
        {   
            int i{0}, j{0};
            T value;
            
            //getting the line
            std::istringstream line(buffer);
            line >> i >> j >> value;

            //checking that the indices are ok
            if(i<=0 or i>this->m() or j<=0 or j>this->n())
                throw std::runtime_error("Wrong dimensions for an element's indices");

            //storing
            this->operator()(i-1,j-1,value);

            //counting how many NNZ is actually storing
            count++;
        }   
    }
    
    //checking that the number of stored elements is the one declared
    //(a little bit redudant since is storing elements dynamically)
    if(count!=num_lines)
    {
        std::cerr << count << "!=" << num_lines << "\n";
        throw std::runtime_error("Matrix file incorrectly formatted");
    }
    
    file.close();

    return true;
}