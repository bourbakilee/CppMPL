#ifndef TRAJECTORY_NN_H
#define TRAJECTORY_NN_H

#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "InitialValueGuess_types.h"

#include "rt_nonfinite.h"
#include "InitialValueGuess.h"
#include "InitialValueGuess_terminate.h"
#include "InitialValueGuess_initialize.h"

#include <limits>
#include <Eigen/Dense>
using namespace Eigen;

namespace TrajectoryNN
{
	inline double __a(VectorXd& p) { return p[0]; }
	inline double __b(VectorXd& p) { return -(11 * p[0] - 18 * p[1] + 9 * p[2] - 2 * p[3]) / (2 * p[4]); }
	inline double __c(VectorXd& p) { return 9 * (2 * p[0] - 5 * p[1] + 4 * p[2] - p[3]) / (2 * p[4] * p[4]); }
	inline double __d(VectorXd& p) { return -9 * (p[0] - 3 * p[1] + 3 * p[2] - p[3]) / (2 * p[4] * p[4] * p[4]); }
	inline double __k(double s, VectorXd& r) { return r[0] + s*(r[1] + s*(r[2] + s*r[3])); }
	inline double __theta(double s, VectorXd& r) { return s*(r[0] + s*(r[1] / 2 + s*(r[2] / 3 + s*r[3] / 4))); }
	void __thetas(VectorXd& thetas, VectorXd& ss, VectorXd& r);
	void __xy(double& x, double& y, double s, VectorXd& r, double ref_size = 8.0);
	void __jacobian(Matrix3d& Jcb, VectorXd& p, VectorXd& r);
	void optimize(VectorXd& p, VectorXd& r, VectorXd& bd_con, int iter_num = 20, double k_m = 0.2);

	void spiral3(VectorXd& r, VectorXd& q0, VectorXd& q1, double k_m = 0.2);
	void spiral3(double r[], double q0[], double q1[], double k_m = 0.2);

	void path(ArrayXXd& points, double r[], double q0[], double q1[], double length = -1, double ref_size = 0.2);


	void velocity(double u[], double v0, double a0, double vg, double sg);
	// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	void trajectory(ArrayXXd& traj, double r[], double u[], double ref_length = 0., double ref_time = 0.);

	const double inf = std::numeric_limits<double>::infinity();
	//weights - (k, dk, v, a, a_c, offset, j)
	const double cost_weights[] = { 5., 1., 1., 10., 10., 0.1, 1. };
	//kinematic_limits - { k_m, dk_m, v_max, v_min, a_max, a_min, ac_m }
	const double kinematic_limits[] = { 0.2, 0.15, 20., 0., 2.1, -6.1, 9. };
	// eval
	double traj_eval(ArrayXXd& traj, double state_g[], const double weights[] = cost_weights, const double limits[] = kinematic_limits);

	//
	struct Traj
	{
		Traj(double state_i[], double state_g[], double a_i=0., const double weights[] = cost_weights, const double limits[] = kinematic_limits);
		~Traj();
		double r[5];
		double u[4];
		double cost;
		ArrayXXd points;
		void interp(double state[], double time);
	};


	Traj opt_traj(double state_i[], double state_g[], double a_i = 0., const double weights[] = cost_weights, const double limits[] = kinematic_limits);
}

#endif
