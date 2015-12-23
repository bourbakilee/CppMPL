#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "spiral3.h"

namespace trajectory {
	struct Configuration {
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
	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	const double cost_weights[] = { 10.,10.,1.,0.1,0.1,10.,0.1,0.1,10.,1. };
	//kinematic_limits - { k_m, dk_m, v_max, v_min, a_max, a_min }
	const double kinematic_limits[] = { 0.2, 0.1, 20., 0., 2., -6. };
	// u: u0,u1,u2. tg
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	void velocity(double u[], double v0, double a0, double vg, double sg);
	void trajectory(ArrayXXd& traj, double r[], double u[], double ref_length=0., double ref_time=0.);

	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	double eval_traj(ArrayXXd& traj, double weights[]);
}

#endif