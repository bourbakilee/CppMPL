#ifndef SPIRAL3_H
#define SPIRAL3_H


#include "config.h"

#ifdef WITH_SQLITE3
#include <sqlite3.h>
#endif

#ifdef WITH_BOOST
#include <boost/math/constants/constants.hpp>
const double pi = boost::math::constants::pi<double>();
const double two_pi = 2.*boost::math::constants::pi<double>();
#else
const double pi = 3.141592653589793;
const double two_pi = 6.283185307179586;
#endif

#include <Eigen/Dense> 
using namespace Eigen;
// const int iter_num = 50;
const double eps1 = 1.e-4;
const double eps2 = 1.e-6;

namespace spiral3 {
	// p - [p0=k0, p1=k(sg/3), p2=k(2*sg/3), p3=k1, p4=sg]
	// r - [a,b,c,d,sg]
	// p[4], s - must be positive,non-negative
	inline double __a(VectorXd& p) { return p[0]; }
	inline double __b(VectorXd& p) { return -(11 * p[0] - 18 * p[1] + 9 * p[2] - 2 * p[3]) / (2 * p[4]); }
	inline double __c(VectorXd& p) { return 9 * (2 * p[0] - 5 * p[1] + 4 * p[2] - p[3]) / (2 * p[4] * p[4]); }
	inline double __d(VectorXd& p) { return -9 * (p[0] - 3 * p[1] + 3 * p[2] - p[3]) / (2 * p[4] * p[4] * p[4]); }
	inline double __k(double s, VectorXd& r) { return r[0] + s*(r[1] + s*(r[2] + s*r[3])); }
	inline double __theta(double s, VectorXd& r) { return s*(r[0] + s*(r[1] / 2 + s*(r[2] / 3 + s*r[3] / 4))); }
	void __thetas(VectorXd& thetas, VectorXd& ss, VectorXd& r);
	void __xy(double& x,double& y,double s, VectorXd& r, double ref_size = 8.0);
	// double __y(double s, double r[], double ref_size = 8.0);
	void __jacobian(Matrix3d& Jcb, VectorXd& p, VectorXd& r);
	void optimize(VectorXd& p, VectorXd& r, VectorXd& bd_con, int iter_num=50, double k_m = 0.2);

#ifdef WITH_SQLITE3
	void select_init_val(VectorXd& p, VectorXd& r, VectorXd& bd_con, sqlite3* db);
	void spiral3(VectorXd& r, VectorXd& q0, VectorXd& q1, sqlite3* db, double k_m=0.2);
	// 接口中的矩阵、向量可以通过Eigen::Map从纯C、C++数组转换而来
	void spiral3(double r[], double q0[], double q1[], sqlite3* db, double k_m = 0.2);
#else
	void spiral3(VectorXd& r, VectorXd& q0, VectorXd& q1, double k_m = 0.2);
	// 接口中的矩阵、向量可以通过Eigen::Map从纯C、C++数组转换而来
	void spiral3(double r[], double q0[], double q1[], double k_m = 0.2);
#endif

	void path(ArrayXXd& points, VectorXd& r, VectorXd& q0, VectorXd& q1, double length = -1, double ref_size = 0.2);
}

#endif