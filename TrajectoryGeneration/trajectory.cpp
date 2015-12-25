#include <cmath>

#include "config.h"
#include "trajectory.h"

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

	
	// TODO: considering add road and obstacle(static and dynamic) items - 2015.12.22
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	// r: a, b, c, d, sg
	// u: u0,u1,u2. tg
	// parameter traj must be pre-computed by spiral3::path procedure, and colums 1-5 must be filled
	void trajectory(ArrayXXd& traj, double r[],double u[], double ref_length, double ref_time)
	{
		int N = traj.rows();
		// traj.resize(N, 11);
		// traj.block(0, 1, N, 5) = points; // 1-5 col

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
		/*
		// 9 col: j
		//traj.col(9) = 2 * u[2]*ArrayXXd::Ones(N,1);
	    // 10 col: al
		//traj.col(10) = traj.col(7)*traj.col(7)*traj.col(5);
		*/
		// ref_time -> absolute time
		traj.col(0) += ref_time;
		// ref_length -> absolute length
		traj.col(1) += ref_length;
	}

	// Eigen3 <-> OpenCV
	// I wouldn't necessarily recommend to use cv2eigen() and eigen2cv() though
	// use Eigen::Map to just map the memory (no cost in copying anything) and cv::Mat(void*, ...) to map the data back
	// MatrixXd matrix; double* arrayd = matrix.data();


	// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	// weights - (k, dk, v, a, a_c, offset, env, j, t, s)
	// kinematic_limits - { k_m, dk_m, v_max, v_min, a_max, a_min }
	double eval_traj(ArrayXXd& traj, const double *weights, const double* k_limits, environment::Vehicle* vehicle, environment::CostMap* cost_map, environment::Road* road)
	{
		int N = traj.rows();
		ArrayXXd cost(N, 1);
		ArrayXXd cost_matrix(N,7); // k, dk, v, a, a_c, off_set(x,y), env(t,x,y,theta)
		cost_matrix.col(0) = weights[0] * traj.col(5).abs(); // |k|
		cost_matrix.col(1) = weights[1] * traj.col(6).abs(); //|dk|
		cost_matrix.col(2) = weights[2] * traj.col(7)*traj.col(7); // v^2
		cost_matrix.col(3) = weights[3] * traj.col(8)*traj.col(8); // a^2
		cost_matrix.col(4) = weights[4] / weights[2] * traj.col(5).abs()*cost_matrix.col(2); // v^2*k
		if (road != nullptr)
		{
			ArrayXXd sl(N, 2);
			road->traj2sl(traj, sl);
			cost_matrix.col(5) = weights[5] * sl.col(1); // road
		}
		else
		{
			cost_matrix.col(5) *= 0.;
		}
		if (cost_map != nullptr)
		{
			if (vehicle == nullptr)
				*vehicle = environment::Vehicle();
			// ArrayXXd cost(N, 1);
			cost_map->query(vehicle, traj, cost);
			cost_matrix.col(6) = weights[6] * cost; // environment
		}
		else
		{
			cost_matrix.col(6) *= 0.;
		}
		// 
		cost = cost_matrix.rowwise().sum();
		//
		// check kinematics limits
		ArrayXXd flag(N, 4);
		flag.col(0) = (traj.col(5).abs() > k_limits[0]).cast<double>(); // k
		flag.col(1) = (traj.col(6).abs() > k_limits[1]).cast<double>(); // dk
		flag.col(2) = (traj.col(7) > k_limits[2] || traj.col(7) < k_limits[3]).cast<double>(); // v
		flag.col(4) = (traj.col(8) > k_limits[4] || traj.col(8) < k_limits[5]).cast<double>(); //a
		ArrayXXd total_flag = flag.rowwise().sum();
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			if (total_flag(i, 0) > 0)
				cost(i, 0) = inf;
		}
		// truncate the trajectory
		double total_cost = 0.;
		int M = 0;
		while (!std::isinf(cost(M, 0)) && M < N)
		{
			total_cost += cost(M, 0);
			M += 1;
		}
		if (M > 0)
		{
			ArrayXXd tmp1 = traj.block(0, 0, M, 9);
			traj.resize(M, 9);
			traj = tmp1;
		}
		else
		{
			ArrayXXd tmp2 = traj.row(0);
			traj.resize(1, 9);
			traj = tmp2;
		}
		return total_cost;
	}
}
