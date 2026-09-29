#include "./lib/qizingxy.h"

#include <iostream>
#include <complex>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>
using namespace std::literals::complex_literals;


int main() {
  Hamilitonian H(2, 1.0, 2.0, 0.5, 2, 3); 
  // H.print();
  return 0;
}

// #include <iostream>
// #include <xtensor/containers/xarray.hpp>
// #include <xtensor/io/xio.hpp>
// int main() {
//   xt::xarray<double> A = {{1, 2}, {3, 4}};
//   xt::xarray<double> B = {{5, 6}, {7, 8}};
//   Spin a(1, 2, 3);
//   std::cout << "S