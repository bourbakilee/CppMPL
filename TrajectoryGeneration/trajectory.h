#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "spiral3.h"

namespace trajectory {
	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	const double cost_weights[] = { 10.,10.,1.,0.1,0.1,10.,0.1,0.1,10.,1. };
	const double kinematic_limits[] = { 0. };
	// u: u0,u1,u2. tg
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	void velocity(double u[], double v0, double a0, double vg, double sg);
	void trajectory(ArrayXXd& traj, double r[], double u[], double ref_length=0., double ref_time=0.);

	//weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	double eval_traj(ArrayXXd& traj, double weights[]);
}

#endif