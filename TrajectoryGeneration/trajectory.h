#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <limits>
// #include "spiral3.h"
// define EIGEN_USE_MKL_ALL
#include <Eigen/Dense> 

#include "environment.h"

namespace trajectory {
	using namespace Eigen;
	/*
	struct Configuration {
		// double s,l; int i,j;
		double x;
		double y;
		double theta;
		double k;
		// construct functions
		Configuration() :Configuration(0., 0., 0., 0.) {}
		Configuration(double x, double y, double t, double k) :x(x), y(y), theta(t), k(k) {}
		// member functions
		void set(double x, double y, double t, double k) { this->x = x; this->y = y; this->theta = t; this->k = k; }
	};
	*/
	// infinity, use isinf() function to determine if a number is infinity
	const double inf = std::numeric_limits<double>::infinity();
	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	const double cost_weights[] = { 5., 10., -0.1, 10., 0.1, 0.1, 50., 5, 40., -4. };
	//kinematic_limits - { k_m, dk_m, v_max, v_min, a_max, a_min, ac_m }
	const double kinematic_limits[] = { 0.2, 0.15, 20., 0., 2.1, -6., 6.};
	// u: u0,u1,u2. tg
	// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	void velocity(double u[], double v0, double a0, double vg, double sg);
	void trajectory(ArrayXXd& traj, double r[], double u[], double ref_length=0., double ref_time=0.);

	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	using Eval_Res = std::pair<double, bool>;
	Eval_Res eval_traj(ArrayXXd& traj, const environment::Vehicle& vehicle, const environment::CostMap& cost_map, const double *weights = cost_weights, const double* k_limits = kinematic_limits, environment::Road* road = nullptr, bool truncate=true);

}

#endif
