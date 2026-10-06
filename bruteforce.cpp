#include "./lib/qizingxy.h"
#define N 2
#include <complex>
#include <fstream>
#include <iostream>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>
using namespace std::literals::complex_literals;

int main() {
  Hamilitonian H(N, 1.0, 2.0, 0.5, 2, 3);
  std::cout << H.rho(1) << std::endl;
  double w = 0;
  std::ofstream MyLog("log.txt");
  for (int i = 0; i < pow(N, 2); i++) {
    std::string bin = std::format("{:0{}b}", i, N);
    MyLog << bin << " " << H.rho(1)(i, i) << std::endl;
  }
  MyLog.close();
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
