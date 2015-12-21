#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "spiral3.h"

namespace trajectory {
	// u: u0,u1,u2. tg
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	void velocity(double u[], double v0, double a0, double vg, double sg);
	void trajectory(ArrayXXd& traj, ArrayXXd& points, double r[], double u[]);
	double eval_traj(ArrayXXd& traj);
}

#endif