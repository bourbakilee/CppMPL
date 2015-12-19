#ifndef SPIRAL3_H
#define SPIRAL3_H

#include <sqlite3.h>
#include <armadillo>
using namespace arma;
#include <boost/math/constants/constants.hpp>
const double pi = boost::math::constants::pi<double>();
const double eps1 = 1.e-4;
const double eps2 = 1.e-6;
const double eps = 1.e-8;
const int iter_num = 50;

namespace spiral3 {
	// p - [p0=k0, p1=k(sg/3), p2=k(2*sg/3), p3=k1, p4=sg]
	// r - [a,b,c,d,sg]
	// p[4], s - must be positive,non-negative
	inline double __a(vec& p) { return p[0]; }
	inline double __b(vec& p) { return -(11 * p[0] - 18 * p[1] + 9 * p[2] - 2 * p[3]) / (2 * p[4]); }
	inline double __c(vec& p) { return 9 * (2 * p[0] - 5 * p[1] + 4 * p[2] - p[3]) / (2 * p[4] * p[4]); }
	inline double __d(vec& p) { return -9 * (p[0] - 3 * p[1] + 3 * p[2] - p[3]) / (2 * p[4] * p[4] * p[4]); }
	inline double __k(double s, vec& r) { return r[0] + s*(r[1] + s*(r[2] + s*r[3])); }
	inline double __theta(double s, vec& r) { return s*(r[0] + s*(r[1] / 2 + s*(r[2] / 3 + s*r[3] / 4))); }
	void __thetas(vec& thetas, vec& ss, vec& r);
	void __xy(double& x,double& y,double s, vec& r, double ref_size = 8.0);
	// double __y(double s, double r[], double ref_size = 8.0);
	void __jacobian(mat& Jcb, vec& p, vec& r);
	void optimize(vec& pp, vec& bd_con, double k_m=0.2);
	void select_init_val(vec& pp, vec& bd_con, sqlite3* db);
	void calc_path(vec& p, vec& r, vec& q0, vec& q1, sqlite3* db, double k_m=0.2);

}

#endif