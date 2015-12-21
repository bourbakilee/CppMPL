#include "trajectory.h"
#include <cmath>

namespace trajectory {
	// u: u0,u1,u2. tg
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	void velocity(double u[], double v0, double a0, double vg, double sg) 
	{
		/*
		# v(t) = u0 + u1*t + u2*t**2
		# return: q0~q2, tg
		u0, u1, u2, tg = None, None, None, None
		delta = (2*v0+vg)**2 + 6*a0*sg
		if delta >= 0:
			u0 = v0
			u1 = a0
			tg = 3*sg/(2*v0+vg) if np.abs(a0)<1.e-6 else (np.sqrt(delta)-2*v0-vg)/a0
			u2 = (vg - v0 - a0*tg)/tg**2
		return (u0,u1,u2,tg)
		*/
		u[0] = v0;
		u[1] = a0;
		double delta = (2 * v0 + vg)*(2 * v0 + vg) + 6 * a0*sg;
		if (delta > eps2)
		{
			if (std::abs(a0) < eps2)
				u[3] = 3 * sg / (2 * v0 + vg);
			else
				u[3] = (std::sqrt(delta) - 2 * v0 - vg) / a0;
		}
		else if (std::abs(delta) < eps2)
		{
			u[3] = (-2 * v0 - vg) / a0;
		}
		else
		{
			u[3] = -1.;
		}
		u[2] = (vg - v0 - a0*u[3]) / (u[3] * u[3]);
	}

	
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	// r: a, b, c, d, sg
	// u: u0,u1,u2. tg
	void trajectory(ArrayXXd& traj, ArrayXXd& points, double r[],double u[])
	{
		int N = points.rows();
		traj.resize(N, 11);
		traj.block(0, 1, N, 5) = points; // 1-5 col
		//0 col: t
		ArrayXXd times = VectorXd::LinSpaced(N, 0., u[3]);
		ArrayXXd lengths = times*(u[0] + times*(u[1] / 2. + times*u[2] / 3.));
	
		traj(0, 0) = 0., traj(N - 1, 0) = u[3];
		int j = 1;
		for (int i = 1; i < N-1; i++)
		{
			while(traj(i, 1) >= lengths(j, 0))
				j++;
			traj(i, 0) = times(j - 1, 0) + (traj(i, 1) - lengths(j - 1, 0))*(times(j, 0) - times(j - 1, 0)) / (lengths(j, 0) - lengths(j - 1, 0));
		}
		// 7 col: v
		traj.col(7) = u[0] + traj.col(0)*(u[1] + traj.col(0)*u[2]);
		// 6 col: dk/dt = v*dk/ds
		traj.col(6) = traj.col(7)*(r[1] + traj.col(1)*(2 * r[2] + 3 * r[3] * traj.col(1)));
		// 8 col: a
		traj.col(8) = u[1] + 2 * u[2] * traj.col(0);
		// 9 col: j
		traj.col(9) = 2 * u[2]*ArrayXXd::Ones(N,1);
	    // 10 col: al
		traj.col(10) = traj.col(7)*traj.col(7)*traj.col(5);
	}


	double eval_traj(ArrayXXd& traj)
	{
	}

}