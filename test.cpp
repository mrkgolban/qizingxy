#include "./lib/qizingxy.h"

#include <iostream>
#include <complex>
#include <string>
#include <format>
#include <fstream>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xarray.hpp>
#include <xtensor/io/xio.hpp>
#include <xtensor/io/xcsv.hpp>
#include <xtensor/misc/xsort.hpp>

using namespace std::literals::complex_literals;

int main() {
    std::ifstream inans;
    std::ifstream intst;
    int N;
    Hamilitonian HH(2, 1.0, 2.0, 0.5, 2, 3);
    xt::xarray<double> a = {{1, 2}, {3, 4}, {5, 6}};
    try {
        HH.set(a);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "You've put wrong matrix" << std::endl;
        return -1;
    }



    for (int l = 0; l < 3; l++)
    {
        intst.open(std::format("./tests/test{}.csv", l + 1));
        inans.open(std::format("./tests/test{}_ans.csv", l + 1));
        auto data_test = xt::load_csv<double>(intst, ' ');
        auto data_ans = xt::load_csv<double>(inans, ' ');

        if (data_test.shape()[0] != data_test.shape()[1])
        {
            throw std::runtime_error("Wrong data");
        }
        Hamilitonian H(data_test.shape()[0], 1.0, 2.0, 0.5, 2, 3);
        H.set(data_test);
        auto result = xt::sort(H.eigenvalues());
        auto ans_result = xt::sort(data_ans);
    
        if (xt::amax(xt::abs(result - ans_result))() < 1e-3)
        {
            std::cout << std::format("Test {} completed", l + 1) << std::endl;
        }
        else 
        {
            std::cout << std::format("Test {} failed", l + 1) << std::endl;
        }

        intst.close();
        inans.close(); 
    }
  return 0;
}