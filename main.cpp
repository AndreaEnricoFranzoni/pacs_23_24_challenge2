#include <vector>
#include <iostream>
#include "chrono.hpp"
#include "matrix.hpp"
#include <concepts>
#include <complex>

int main()
{

    int status(0);
    using namespace std;

    //types used
    using T = double;
    using V = complex<double>;

    //matrix in market format to be exploited
    std::string matrix_mmf = "lnsp_131.mtx";

    //storaging
    constexpr StorageOrder ro = StorageOrder::RowWise;
    constexpr StorageOrder co = StorageOrder::ColumnWise;

    //storing the matrix rw and cw
    algebra::Matrix<T,ro> M_rw;
    algebra::Matrix<T,co> M_cw;

    M_rw.reader_mmf(matrix_mmf);
    M_cw.reader_mmf(matrix_mmf);

    
    //vector to test the product
    std::size_t dim_vec = M_rw.n();

    //std::vector
    std::vector<T> v_vec;
    v_vec.reserve(dim_vec);

    //algebra::Matrix<T,O>
    algebra::Matrix<T,ro> v_mat(dim_vec,1);

    //filling them
    for (size_t i = 0; i < dim_vec; ++i)
    {
        v_vec.emplace_back(i);
        v_mat(i,0,i);
    }


    //CHECKING THE ROW WISE STORAGE
    cout << "ROWISE" << endl;
    cout << "Checking the performance in uncompressed format" << endl;

    cout << "Matrix*vector with std::vector" << endl;
    Timings::Chrono timer1;
    timer1.start();
    M_rw*v_vec;
    timer1.stop();
    cout << timer1 << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer2;
    timer2.start();
    M_rw*v_mat;
    timer2.stop();
    cout << timer2 << endl;

    cout << "Matrix*Matrix with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer3;
    timer3.start();
    M_rw*M_rw;
    timer3.stop();
    cout << timer3 << endl;

    cout << "Norm one" << endl;
    Timings::Chrono timer4;
    timer4.start();
    M_rw.norm<NormType::One>();
    timer4.stop();
    cout << timer4 << endl;
    
    cout << "Norm infinity" << endl;
    Timings::Chrono timer5;
    timer5.start();
    M_rw.norm<NormType::Infinity>();
    timer5.stop();
    cout << timer5 << endl;

    cout << "Norm Frobenius" << endl;
    Timings::Chrono timer6;
    timer6.start();
    M_rw.norm<NormType::Frobenius>();
    timer6.stop();
    cout << timer6 << endl;


    cout << "******" << endl;
    cout << "Checking the performance in compressed format" << endl;
    M_rw.compress();

    cout << "Matrix*vector with std::vector" << endl;
    Timings::Chrono timer7;
    timer7.start();
    M_rw*v_vec;
    timer7.stop();
    cout << timer7 << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer8;
    timer8.start();
    M_rw*v_mat;
    timer8.stop();
    cout << timer8 << endl;

    cout << "Matrix*Matrix with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer9;
    timer9.start();
    M_rw*M_rw;
    timer9.stop();
    cout << timer9 << endl;

    cout << "Norm one" << endl;
    Timings::Chrono timer10;
    timer10.start();
    M_rw.norm<NormType::One>();
    timer10.stop();
    cout << timer10 << endl;

    cout << "Norm infinity" << endl;
    Timings::Chrono timer11;
    timer11.start();
    M_rw.norm<NormType::Infinity>();
    timer11.stop();
    cout << timer11 << endl;

    cout << "Norm Frobenius" << endl;
    Timings::Chrono timer12;
    timer12.start();
    M_rw.norm<NormType::Frobenius>();
    timer12.stop();
    cout << timer12 << endl;


    cout << "******" << endl;
    cout << "******" << endl;

    //CHECKING THE COL WISE STORAGE
    cout << "COLWISE" << endl;
    cout << "Checking the performance in uncompressed format" << endl;

    cout << "Matrix*vector with std::vector" << endl;
    Timings::Chrono timer13;
    timer13.start();
    M_cw*v_vec;
    timer13.stop();
    cout << timer13 << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer14;
    timer14.start();
    M_cw*v_mat;
    timer14.stop();
    cout << timer14 << endl;

    cout << "Matrix*Matrix with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer15;
    timer15.start();
    M_cw*M_cw;
    timer15.stop();
    cout << timer15 << endl;

    cout << "Norm one" << endl;
    Timings::Chrono timer16;
    timer16.start();
    M_cw.norm<NormType::One>();
    timer16.stop();
    cout << timer16 << endl;

    cout << "Norm infinity" << endl;
    Timings::Chrono timer17;
    timer17.start();
    M_cw.norm<NormType::Infinity>();
    timer17.stop();
    cout << timer17 << endl;

    cout << "Norm Frobenius" << endl;
    Timings::Chrono timer18;
    timer18.start();
    M_cw.norm<NormType::Frobenius>();
    timer18.stop();
    cout << timer18 << endl;


    cout << "******" << endl;
    cout << "Checking the performance in compressed format" << endl;
    M_cw.compress();

    cout << "Matrix*vector with std::vector" << endl;
    Timings::Chrono timer19;
    timer19.start();
    M_cw*v_vec;
    timer19.stop();
    cout << timer19 << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer20;
    timer20.start();
    M_cw*v_mat;
    timer20.stop();
    cout << timer20 << endl;

    cout << "Matrix*Matrix with algebra::Matrix<T,O>" << endl;
    Timings::Chrono timer21;
    timer21.start();
    M_cw*M_cw;
    timer21.stop();
    cout << timer21 << endl;

    cout << "Norm one" << endl;
    Timings::Chrono timer22;
    timer22.start();
    M_cw.norm<NormType::One>();
    timer22.stop();
    cout << timer22 << endl;

    cout << "Norm infinity" << endl;
    Timings::Chrono timer23;
    timer23.start();
    M_cw.norm<NormType::Infinity>();
    timer23.stop();
    cout << timer23 << endl;

    cout << "Norm Frobenius" << endl;
    Timings::Chrono timer24;
    timer24.start();
    M_cw.norm<NormType::Frobenius>();
    timer24.stop();
    cout << timer24 << endl;


    cout << "******" << endl;
    cout << "******" << endl;

    //CHECKING HAVING COMPLEX VALUES
    cout << "COMPLEX VALUES" << endl;

    V a11 = 1.0 + 1.0i;
    V a12 = 2.0;
    V a21 = 1.0 + 3.0i;
    V a22 = 5.0 - 2.0i;

    V v1 = 2.0i;
    V v2 = 2.0 + 1.0i;


    algebra::Matrix<V,StorageOrder::RowWise> Mc_r(2,2);
    Mc_r(0,0,a11);
    Mc_r(0,1,a12);
    Mc_r(1,0,a21);
    Mc_r(1,1,a22);

    algebra::Matrix<V,StorageOrder::ColumnWise> Mc_c(2,2);
    Mc_c(0,0,a11);
    Mc_c(0,1,a12);
    Mc_c(1,0,a21);
    Mc_c(1,1,a22);


    std::vector<V> vc_vec;
    vc_vec.reserve(2);
    vc_vec.emplace_back(v1);
    vc_vec.emplace_back(v2);

    algebra::Matrix<V,StorageOrder::RowWise> vc_mat_r(2,1);
    vc_mat_r(0,0,v1);
    vc_mat_r(1,0,v2);

    algebra::Matrix<V,StorageOrder::ColumnWise> vc_mat_c(2,1);
    vc_mat_c(0,0,v1);
    vc_mat_c(1,0,v2);

    
    cout << "The matrix is " << endl;
    cout << Mc_r << endl;

    cout << "The vector is " << endl;
    cout << vc_mat_r << endl;

    
    cout << "ROW WISE UNCOMPRESSED" << endl;

    cout << "Matrix*vector with std::vector: the result is " << endl;
    std::vector<V> res_ru = Mc_r*vc_vec;
    cout << res_ru[0] << endl;
    cout << res_ru[1] << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_r*vc_mat_r << endl;
    //cout << Mc_r*vc_mat_c << endl;    //to check the product with the other storaging
    
    cout << "Square of the matrix with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_r*Mc_r << endl;
    //cout << Mc_r*Mc_c << endl;          //to check the product with the other storaging

    cout << "Norm one" << endl;
    cout << Mc_r.norm<NormType::One>() << endl;

    cout << "Norm infinity" << endl;
    cout << Mc_r.norm<NormType::Infinity>() << endl;

    cout << "Norm Frobenius" << endl;
    cout << Mc_r.norm<NormType::Frobenius>() << endl;


    /*    
    cout << "" << endl;
    cout << "ROW WISE COMPRESSED" << endl;
    Mc_r.compress();

    cout << "Matrix*vector with std::vector: the result is " << endl;
    std::vector<V> res_rc = Mc_r*vc_vec;
    cout << res_rc[0] << endl;
    cout << res_rc[1] << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_r*vc_mat_r << endl;
    //cout << Mc_r*vc_mat_c << endl;    //to check the product with the other storaging

    cout << "Square of the matrix with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_r*Mc_r << endl;
    //cout << Mc_r*Mc_c << endl;          //to check the product with the other storaging


    cout << "Norm one" << endl;
    cout << Mc_r.norm<NormType::One>() << endl;

    cout << "Norm infinity" << endl;
    cout << Mc_r.norm<NormType::Infinity>() << endl;

    cout << "Norm Frobenius" << endl;
    cout << Mc_r.norm<NormType::Frobenius>() << endl;
    */
    

    /*
    cout << "" << endl;
    cout << "COL WISE UNCOMPRESSED" << endl;

    cout << "Matrix*vector with std::vector: the result is " << endl;
    std::vector<V> res_cu = Mc_c*vc_vec;
    cout << res_cu[0] << endl;
    cout << res_cu[1] << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_c*vc_mat_c << endl;
    //cout << Mc_c*vc_mat_r << endl;      //to check the product with the other storaging                  
    
    cout << "Square of the matrix with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_c*Mc_c << endl;
    //cout << Mc_c*Mc_r << endl;          //to check the product with the other storaging

    cout << "Norm one" << endl;
    cout << Mc_c.norm<NormType::One>() << endl;

    cout << "Norm infinity" << endl;
    cout << Mc_c.norm<NormType::Infinity>() << endl;

    cout << "Norm Frobenius" << endl;
    cout << Mc_c.norm<NormType::Frobenius>() << endl;
    */


    /*
    cout << "" << endl;
    cout << "COL WISE COMPRESSED" << endl;
    Mc_c.compress();

    cout << "Matrix*vector with std::vector: the result is " << endl;
    std::vector<V> res_cc = Mc_c*vc_vec;
    cout << res_cc[0] << endl;
    cout << res_cc[1] << endl;

    cout << "Matrix*vector with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_c*vc_mat_c << endl;
    //cout << Mc_c*vc_mat_r << endl;      //to check the product with the other storaging                  
    
    cout << "Square of the matrix with algebra::Matrix<T,O>: the result is " << endl;
    cout << Mc_c*Mc_c << endl;
    //cout << Mc_c*Mc_r << endl;          //to check the product with the other storaging

    cout << "Norm one" << endl;
    cout << Mc_c.norm<NormType::One>() << endl;

    cout << "Norm infinity" << endl;
    cout << Mc_c.norm<NormType::Infinity>() << endl;

    cout << "Norm Frobenius" << endl;
    cout << Mc_c.norm<NormType::Frobenius>() << endl;
    */


   return status;
}

    

   



 


    

 
    

  


    

  
   