#pragma once

// #include <xtensor/xutils.hpp>
#include <complex>
#include <math.h>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/containers/xtensor.hpp>
#include <xtensor/io/xio.hpp>

typedef xt::xarray<std::complex<double>> arcomp;
typedef xt::xarray<double> ardouble;
enum class Axis { X, Y, Z };

using namespace std::literals::complex_literals;

class Sigma {
private:
  arcomp mat;
  int i;
  Axis ax;
  int N;

public:
  Sigma(int ii, int NN, Axis axx);
  arcomp get_mat();
};

class Hamilitonian {
private:
  // H = -J * ((1 + gamma) * sigma^x_i+1 * sigma^x_i + (1-gamma) * sigma^y_i+1 *
  // sigma^y_i) - h sigma^z_i
  long unsigned int N;
  long unsigned int STATES;
  int height;
  int length;
  double h;
  double J;
  double gamma; // -1 <= gamma <= 1
  arcomp ham;
  ardouble evectors;
  ardouble energy_spectr;

public:
  Hamilitonian(int NN, double hh, double JJ, double ggama, int hheight,
               int llength);
  void change_grid(int new_length, int new_height);
  void change_amount(int newN);
  void magn_field(double nh, double nJ, double ngamma);
  void CalcHam();
  ardouble rho(double T);
  void print() const;
  ardouble espectra() const;
  ardouble eigenvectors() const;
};
