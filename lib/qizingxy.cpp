#include "./qizingxy.h"
#include <string>


// Hamiltonian realization
Hamilitonian::Hamilitonian(int NN, double hh, double JJ, double ggamma,
                           int hheight, int llength): N((long unsigned int)NN), h(hh), J(JJ), gamma(ggamma),
      length(llength), height(hheight) {

  ham = xt::zeros<std::complex<double>>(
      {static_cast<std::size_t>(N), static_cast<std::size_t>(N)});
  ham.fill(0);
  // change_grid(length, height);

    int i_max = N / length;
    int index = 0;
    int j_max = (N % length > 0 && i_max > 0) ? N % length : length - 2;
    std::cout << xt::imag(Sigma(0, 3, Axis::Y).get_mat()) << std::endl;    

};


void Hamilitonian::print() { std::cout << ham << std::endl; }

// void Hamilitonian::set(xt::xarray<std::complex<double>> arg)
// {
//   if (ham.shape() != arg.shape())
//   throw std::invalid_argument("Shapes aren't the same\n");
//   ham = arg;
// }

xt::xarray<double> Hamilitonian::eigenvalues() {
   return xt::real(xt::linalg::eigvals(ham));
}



// Simga realization

Sigma::Sigma(int ii, int NN, Axis axx): i(ii), ax(axx), N(NN) 
{ 
  mat = xt::zeros<std::complex<double>>(
      {static_cast<std::size_t>(pow(2, N)), static_cast<std::size_t>(pow(2, N))});
      switch(axx) {
        case Axis::X:
        for (int index = 0; index < pow(2, N); index++)
        {
          int reversed = index ^ (1 << i);
          std::cout << std::format("{:0{}b}", index, N) << " " << std::format("{:0{}b}", reversed, N) << std::endl;
          mat(index, reversed) = 1;
        }
        break;
        case Axis::Y:
        for (int index = 0; index < pow(2, N); index++)
        {
          int newState = index ^ (1 << i);
          float sign = ((index >> i) & 1) ? -1.0 : 1.0;
          std::complex<double> coeff(0, sign);
          mat(index, newState) = coeff;
        }
        break;
        case Axis::Z:
        for (int index = 0; index < pow(2, N); index++)
        {
          std::string bin = std::format("{:0{}b}", index, N);
          mat(index, index) = (char)bin[N - i] == '0' ? 1 : -1;
          std::cout << mat(index, index) << std::endl;
        }
      break;
    default:
      throw std::invalid_argument("Wrong axis");
  }
};

arcomp Sigma::get_mat() {
  return mat;
}