#include "./qizingxy.h"
#include <string>

// Hamiltonian realization
Hamilitonian::Hamilitonian(int NN, double hh, double JJ, double ggamma,
                           int hheight, int llength)
    : N((long unsigned int)NN), h(hh), J(JJ), gamma(ggamma), length(llength),
      height(hheight) {
  CalcHam();
};

void Hamilitonian::CalcHam() {
  STATES = pow(2, N);
  ham = xt::zeros<std::complex<double>>({static_cast<std::size_t>(pow(2, N)),
                                         static_cast<std::size_t>(pow(2, N))});
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
      ham -= h * Sigma(index, N, Axis::Z).get_mat();
    }
  }
  auto [es, ev] = xt::linalg::eigh(ham);
  energy_spectr = es;
  evectors = xt::real(ev);
}

void Hamilitonian::print() const { std::cout << ham << std::endl; }

ardouble Hamilitonian::espectra() const { return energy_spectr; }

ardouble Hamilitonian::eigenvectors() const { return evectors; }

void Hamilitonian::change_amount(int newN) {
  N = newN;
  CalcHam();
}

void Hamilitonian::change_grid(int new_length, int new_height) {
  length = new_length;
  height = new_height;
  CalcHam();
}

void Hamilitonian::magn_field(double nh, double nJ, double ngamma) {
  h = nh;
  J = nJ;
  gamma = ngamma;
  CalcHam();
}

ardouble Hamilitonian::rho(double T) {
  ardouble D = xt::zeros<double>({pow(2, N), pow(2, N)});
  long double Z = 0;
  for (int i = 0; i < pow(2, N); i++) {
    D(i, i) = std::exp((1 / T) * energy_spectr[i]);
    Z += D(i, i);
  }

  return (1 / Z) *
         xt::linalg::dot(xt::linalg::dot(evectors, D), xt::transpose(evectors));
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
