#ifndef SPIRAL3_CONFIG_H
#define SPIRAL3_CONFIG_H

// user defined options
#define WITH_SQLITE3
// #define WITH_BOOST

#ifdef WITH_BOOST
#include <boost/math/constants/constants.hpp>
const double pi = boost::math::constants::pi<double>();
const double two_pi = 2.*boost::math::constants::pi<double>();
#else
const double pi = 3.141592653589793;
const double two_pi = 6.283185307179586;
#endif

// const int iter_num = 50;
const double eps1 = 1.e-4;
const double eps2 = 1.e-6;

// define EIGEN_USE_MKL_ALL

#endif //
