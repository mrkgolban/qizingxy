#include "./qizingxy.h"
#include <string>

// Hamiltonian realization
Hamilitonian::Hamilitonian(int NN, double hh, double JJ, double ggamma,
                           int hheight, int llength)
    : N((long unsigned int)NN), h(hh), J(JJ), gamma(ggamma), length(llength),
      height(hheight) {

  ham = xt::zeros<std::complex<double>>({static_cast<std::size_t>(pow(2, N)),
                                         static_cast<std::size_t>(pow(2, N))});
  ham.fill(0);
  // change_grid(length, height);

  int i_max = N / length;
  int j_max = (N % length > 0 && i_max > 0) ? N % length : length - 2;
  int index;
  int j = 0;
  double Jx = J * (1 + gamma);
  double Jy = J * (1 - gamma);
  for (int i = 0; i < height && i * length + j < N; i++) {
    for (; j < length && i * length + j < N; j++) {
      index = i * length + j;
      if (index + length < N) {
        ham -= Jy * xt::linalg::dot(Sigma(N - index - 1, N, Axis::Y).get_mat(),
                                    Sigma(index, N, Axis::Y).get_mat());
      }
      if ((j + 1 < length - 1 && i != height - 1) ||
          (i == height - 1 && index + 1 < N)) {
        ham -= Jx * xt::linalg::dot(Sigma(N - index - 1, N, Axis::X).get_mat(),
                                    Sigma(index, N, Axis::X).get_mat());
      }
      std::cout << index << " Z" << std::endl;
      std::cout << Sigma(index, N, Axis::Z).get_mat() << std::endl;
      ham -= h * Sigma(index, N, Axis::Z).get_mat();
    }
  }
};

void Hamilitonian::print() { std::cout << ham << std::endl; }

// void Hamilitonian::set(xt::xarray<std::complex<double>> arg)
// {
//   if (ham.shape() != arg.shape())
//   oooikk
//   throw std::invalid_argument("Shapes aren't the same\n");
//   ham = arg;
// }

xt::xarray<double> Hamilitonian::eigenvalues() {
  return xt::real(xt::linalg::eigvals(ham));
}

// Simga realization

Sigma::Sigma(int ii, int NN, Axis axx) : i(ii), ax(axx), N(NN) {
  mat = xt::zeros<std::complex<double>>({static_cast<std::size_t>(pow(2, N)),
                                         static_cast<std::size_t>(pow(2, N))});
  switch (axx) {
  case Axis::X:
    for (int index = 0; index < pow(2, N); index++) {
      int reversed = index ^ (1 << i);
      mat(index, reversed) = 1;
    }
    break;
  case Axis::Y:
    for (int index = 0; index < pow(2, N); index++) {
      int newState = index ^ (1 << i);
      float sign = ((index >> i) & 1) ? -1.0 : 1.0;
      std::complex<double> coeff(0, sign);
      mat(index, newState) = coeff;
    }
    break;
  case Axis::Z:
    for (int index = 0; index < pow(2, N); index++) {
      std::string bin = std::format("{:0{}b}", index, N);
      mat(index, index) = (char)bin[i] == '0' ? 1 : -1;
    }
    break;
  default:
    throw std::invalid_argument("Wrong axis");
  }
};

arcomp Sigma::get_mat() { return mat; }
