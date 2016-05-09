#include "TrajectoryNN2.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
// #include <iostream>


const double pi = 3.141592653589793;
const double two_pi = 6.283185307179586;
const double eps1 = 1.e-4;
const double eps2 = 1.e-6;


// -pi - pi
inline double mod2pi(double theta)
{
	double v = std::fmod(theta, two_pi);
	if (v < -pi)
		v += two_pi;
	else if (v > pi)
		v -= two_pi;
	return v;
}

// thetas与ss大小一致
inline void TrajectoryNN2::__thetas(VectorXd& thetas, VectorXd& ss, VectorXd& r)
{
	for (int i = 0; i < ss.rows(); i++)
	{
		thetas[i] = TrajectoryNN2::__theta(ss[i], r);
	}
}

inline void TrajectoryNN2::__xy(double& x, double& y, double s, VectorXd& r, double ref_size)
{
	if (s < 0) { s = r[4]; }
	int N = (int)std::ceil(s / ref_size);
	VectorXd s_list = VectorXd::LinSpaced(N + 1, 0., s);
	//s_list= VectorXd::LinSpaced(N + 1,0., s);
	double h = s / (ref_size * N);
	VectorXd thetas(9), fs(9), gs(9);
	x = y = 0.;
	for (int i = 0; i < N; i++) {
		VectorXd ss = VectorXd::LinSpaced(9, s_list[i], s_list[i + 1]);
		TrajectoryNN2::__thetas(thetas, ss, r);
		fs = thetas.array().cos();
		gs = thetas.array().sin();
		x += h / 3 * (fs[0] + 2 * (fs[2] + fs[4] + fs[6]) + 4 * (fs[1] + fs[3] + fs[5] + fs[7]) + fs[8]);
		y += h / 3 * (gs[0] + 2 * (gs[2] + gs[4] + gs[6]) + 4 * (gs[1] + gs[3] + gs[5] + gs[7]) + gs[8]);
	}
}


void TrajectoryNN2::__jacobian(Matrix3d& Jcb, VectorXd& p, VectorXd& r)
{
	VectorXd ss = VectorXd::LinSpaced(9, 0., p[4]), thetas(9);
	TrajectoryNN2::__thetas(thetas, ss, r);
	VectorXd cos_t = thetas.array().cos();
	VectorXd sin_t = thetas.array().sin();
	VectorXd c_p_1(9);
	c_p_1 << 0., 1851 / 8192., 363 / 1024., 9963 / 8192., 51 / 64., 14475 / 8192., 891 / 1024., 13083 / 8192., 3 / 8.;
	c_p_1 *= p[4] * p[4] / 24;
	VectorXd c_p_2(9);
	c_p_2 << 0., 795 / 8192., 123 / 1024., 2187 / 8192., 3 / 64., -2325 / 8192., -405 / 1024., -10437 / 8192., -3 / 8.;
	c_p_2 *= p[4] * p[4] / 24;
	VectorXd c_s_1(9);
	c_s_1 << 0., (2871 * p[0] + 1851 * p[1] - 795 * p[2] + 169 * p[3]) / 8192.,
		(247 * p[0] + 363 * p[1] - 123 * p[2] + 25 * p[3]) / 1024.,
		(4071 * p[0] + 9963 * p[1] - 2187 * p[2] + 441 * p[3]) / 8192.,
		(15 * p[0] + 51 * p[1] - 3 * p[2] + p[3]) / 64.,
		(3655 * p[0] + 14475 * p[1] + 2325 * p[2] + 25 * p[3]) / 8192.,
		(231 * p[0] + 891 * p[1] + 405 * p[2] + 9 * p[3]) / 1024.,
		(3927 * p[0] + 13083 * p[1] + 10437 * p[2] + 1225 * p[3]) / 8192.,
		(p[0] + 3 * p[1] + 3 * p[2] + p[3]) / 8.;
	c_s_1 *= p[4] / 24;
	VectorXd c_s_2(9);
	c_s_2 << 1 / 24., 1 / 6., 1 / 12., 1 / 6., 1 / 12., 1 / 6., 1 / 12., 1 / 6., 1 / 24.;
	Jcb << -c_p_1.dot(sin_t), c_p_2.dot(sin_t), c_s_2.dot(cos_t) - c_s_1.dot(sin_t),
		c_p_1.dot(cos_t), -c_p_2.dot(cos_t), c_s_2.dot(sin_t) + c_s_1.dot(cos_t),
		0.375 * p[4], 0.375 * p[4], c_s_1[8];
}

// bd_con: boudary conditions [k0,x1,y1,theta1,k1]
// pp: initial value of parameters [p1,p2,sg]
void TrajectoryNN2::optimize(VectorXd& p, VectorXd& r, VectorXd& bd_con, int iter_num, double k_m)
{
	Matrix3d Jcb;
	VectorXd pp(3);
	pp << p[1], p[2], p[4];
	pp[0] = std::max(std::min(pp[0], k_m), -k_m);
	pp[1] = std::max(std::min(pp[1], k_m), -k_m);
	pp[2] = std::max(std::min(pp[2], 1000.), 1.);
	/*
	if (pp[2] <= 0)
	{
		pp[0] = (2 * bd_con[0] + bd_con[4]) / 3.;
		pp[1] = (bd_con[0] + 2 * bd_con[4]) / 3.;
		// pp[2] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::min(std::abs(bd_con[3]), two_pi - std::abs(bd_con[3]));
		pp[2] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::abs(bd_con[3]);
	}
	*/
	VectorXd q_g(3), q_p(3);
	q_g << bd_con[1], bd_con[2], bd_con[3]; // x1,y1,theta1
											// VectorXd p{ bd_con[0], pp[0], pp[1], bd_con[4], pp[2] };
											// VectorXd r{ TrajectoryNN2::__a(p),TrajectoryNN2::__b(p), TrajectoryNN2::__c(p), TrajectoryNN2::__d(p), p[4] };
	VectorXd dq(3);
	dq << 1., 1., 1.;
	int times = 0;
	double x_p = 0., y_p = 0., theta_p = 0.;
	while ((std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2) && (times < iter_num)) {
		times += 1;
		//std::cout << "iteration times: " << times << std::endl;
		TrajectoryNN2::__jacobian(Jcb, p, r);
		theta_p = TrajectoryNN2::__theta(p[4], r);
		TrajectoryNN2::__xy(x_p, y_p, p[4], r);
		q_p << x_p, y_p, theta_p;
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		dq = q_g - q_p;
		// dq[2] = dq[2];
		//
		pp += Jcb.colPivHouseholderQr().solve(dq);
		//pp[0] = p[1] = sgn(pp[0])*std::min(k_m, std::abs(pp[0]));
		pp[0] = p[1] = std::max(std::min(pp[0], k_m), -k_m);
		//pp[1] = p[2] = sgn(pp[1])*std::min(k_m, std::abs(pp[1]));
		pp[1] = p[2] = std::max(std::min(pp[1], k_m), -k_m);
		pp[2] = r[4] = p[4] = std::max(std::min(pp[2], 1000.), 1.);
		r[0] = TrajectoryNN2::__a(p);
		r[1] = TrajectoryNN2::__b(p);
		r[2] = TrajectoryNN2::__c(p);
		r[3] = TrajectoryNN2::__d(p);
	}
	if(times == iter_num) //(std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2)
	{
		r[4] = -1.; //通过检查r[4]来判断结果是否有效
	}
}

void TrajectoryNN2::spiral3(VectorXd & r, VectorXd & q0, VectorXd & q1, double k_m)
{
	//bd_con
	double cc = cos(q0[2]);
	double ss = sin(q0[2]);
	double theta_r = std::fmod(q1[2] - q0[2], two_pi);
	if (theta_r > pi)
		theta_r -= two_pi;
	else if (theta_r < -pi)
		theta_r += two_pi;
	VectorXd bd_con(5);
	bd_con << q0[3], (q1[0] - q0[0])*cc + (q1[1] - q0[1])*ss, -(q1[0] - q0[0])*ss + (q1[1] - q0[1])*cc, theta_r, q1[3];

	//initilize p and r
	double pp[3] = { 0.,0.,0. };
	InitialValueGuess2(bd_con.data(), pp);  // before spiral3 be called, InitialValueGuess_initialize() must be called; after InitialValueGuess_terminate();
	VectorXd p(5);
	p << q0[3], pp[0], pp[1], q1[3], pp[2];
	r[0] = TrajectoryNN2::__a(p);
	r[1] = TrajectoryNN2::__b(p);
	r[2] = TrajectoryNN2::__c(p);
	r[3] = TrajectoryNN2::__d(p);
	r[4] = std::max(p[4], 0.);

	// optimize
	TrajectoryNN2::optimize(p, r, bd_con);
}

void TrajectoryNN2::spiral3(double r[], double q0[], double q1[], double k_m)
{
	VectorXd rr(5), qs(4), qg(4);
	//rr << r[0], r[1], r[2], r[3], r[4];
	qs << q0[0], q0[1], q0[2], q0[3];
	qg << q1[0], q1[1], q1[2], q1[3];
	TrajectoryNN2::spiral3(rr, qs, qg, k_m);
	r[0] = rr[0];
	r[1] = rr[1];
	r[2] = rr[2];
	r[3] = rr[3];
	r[4] = rr[4];
}


// points - {t, s, x, y, theta, k, dk, v, a}, here only 1-5 cols are used
void TrajectoryNN2::path(ArrayXXd& points, double r[], double q0[], double q1[], double length, double ref_size)
{
	if (length < 0) { length = r[4]; }
	int N = (int)std::ceil(length / ref_size);
	double delta = length / N;
	points.resize(N + 1, 9); // t, s, x, y, theta, k, dk, v, a
	points.col(1) = VectorXd::LinSpaced(N + 1, 0., length); // s
	points.col(5) = r[0] + points.col(1)*(r[1] + points.col(1)*(r[2] + points.col(1)*r[3])); // k
	points.col(4) = q0[2] + points.col(1)*(r[0] + points.col(1)*(r[1] / 2 + points.col(1)*(r[2] / 3 + points.col(1)*r[3] / 4))); // theta
	ArrayXXd cos_t = points.col(4).cos();
	ArrayXXd sin_t = points.col(4).sin();
	ArrayXXd d_x = (cos_t.block(0, 0, N, 1) + cos_t.block(1, 0, N, 1))*delta / 2.;
	ArrayXXd d_y = (sin_t.block(0, 0, N, 1) + sin_t.block(1, 0, N, 1))*delta / 2.;
	points(0, 2) = q0[0];
	points(0, 3) = q0[1];
	for (int i = 1; i <= N; i++)
	{
		points(i, 2) = points(i - 1, 2) + d_x(i - 1);
		points(i, 3) = points(i - 1, 3) + d_y(i - 1);
	}
	points(N, 2) = q1[0];
	points(N, 3) = q1[1];
	points(N, 4) = q1[2];
	points(N, 5) = q1[3];
}

void TrajectoryNN2::velocity(double u[], double v0, double a0, double vg, double sg)
{
	u[0] = v0;
	u[1] = a0;
	u[2] = u[3] = -1.;
	double v_tmp = 2 * v0 + vg;
	double a0_abs = std::abs(a0);
	double delta = v_tmp*v_tmp + 6 * a0*sg;
	double delta_abs = std::abs(delta);
	if (delta > eps1)
	{
		if (a0_abs > eps1)
		{
			u[3] = (std::sqrt(delta) - v_tmp) / a0;
		}
		else if (a0_abs <= eps1 && v_tmp > eps1)
		{
			u[3] = 3 * sg / v_tmp;
		}
	}
	else if (delta_abs < eps1 && a0_abs > eps1)
	{
		u[3] = (-2 * v0 - vg) / a0;
	}
	if (u[3] > 0.)
	{
		u[2] = (vg - v0 - a0*u[3]) / (u[3] * u[3]);
	}
}


// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
// r: a, b, c, d, sg
// u: u0,u1,u2. tg
// parameter traj must be pre-computed by spiral3::path procedure, and colums 1-5 must be filled
void TrajectoryNN2::trajectory(ArrayXXd& traj, double r[], double u[], double ref_length, double ref_time)
{
	int N = traj.rows();
	// traj.resize(N, 11);
	// traj.block(0, 1, N, 5) = points; // 1-5 col

	//0 col: t
	ArrayXXd times = VectorXd::LinSpaced(N, 0., u[3]);
	ArrayXXd lengths = times*(u[0] + times*(u[1] / 2. + times*u[2] / 3.));

	traj(0, 0) = 0., traj(N - 1, 0) = u[3];
	int j = 1;
	for (int i = 1; i < N - 1; i++)
	{
		while (traj(i, 1) >= lengths(j, 0) && j < N - 1)
			j++;
		traj(i, 0) = times(j - 1, 0) + (traj(i, 1) - lengths(j - 1, 0))*(times(j, 0) - times(j - 1, 0)) / (lengths(j, 0) - lengths(j - 1, 0));
	}
	// 7 col: v
	traj.col(7) = u[0] + traj.col(0)*(u[1] + traj.col(0)*u[2]);
	// 6 col: dk/dt = v*dk/ds
	traj.col(6) = traj.col(7)*(r[1] + traj.col(1)*(2 * r[2] + 3 * r[3] * traj.col(1)));
	// 8 col: a
	traj.col(8) = u[1] + 2 * u[2] * traj.col(0);

	// ref_time -> absolute time
	traj.col(0) += ref_time;
	// ref_length -> absolute length
	traj.col(1) += ref_length;
}

double TrajectoryNN2::traj_eval(ArrayXXd & traj, double state_g[], const double weights[], const double limits[])
{
	int N = traj.rows();
	ArrayXXd cost_matrix = ArrayXXd::Zero(N, 6); // (k, dk, v, a, a_c, offset)
	cost_matrix.col(0) = weights[0] * traj.col(5).abs(); // |k|
	// cost_matrix.col(1) = weights[1] * traj.col(6).abs(); //|dk|
	// cost_matrix.col(2) = weights[2] * traj.col(7).abs(); // |v|
	cost_matrix.col(3) = weights[3] * traj.col(8)*traj.col(8); // |a|
	cost_matrix.col(4) = weights[4] * traj.col(5)* traj.col(5)*traj.col(7)*traj.col(7)*traj.col(7)*traj.col(7); // a_c = |v^2*k|
	cost_matrix.col(5) = weights[5] * ((state_g[0] - traj.col(2))*sin(state_g[2]) - (state_g[1] - traj.col(3))*cos(state_g[2])).abs(); //offset
	// ArrayXXd lane_offset = ArrayXXd::Zero(N, 1);

	// check kinematics limits
	ArrayXXd flag(N, 5);
	flag.col(0) = (traj.col(5).abs() > limits[0]).cast<double>(); // k
	flag.col(1) = (traj.col(6).abs() > limits[1]).cast<double>(); // dk
	flag.col(2) = (traj.col(7) > limits[2] || traj.col(7) < limits[3]).cast<double>(); // v
	flag.col(3) = (traj.col(8) > limits[4] || traj.col(8) < limits[5]).cast<double>(); //a
	flag.col(4) = (cost_matrix.col(4).abs() > weights[4] * limits[6] * limits[6]).cast<double>(); // a_c

	if(flag.sum()>0.9)
	{
		return inf;
	}
	double jerk = (traj(1, 8) - traj(0, 8)) / (traj(1, 0) - traj(0, 0));
	double length = traj(N - 1, 1) - traj(0, 1);
	double delta_s = traj(1, 1) - traj(0, 1);
	return cost_matrix.sum()*delta_s + weights[6] * jerk*jerk*length;
}



TrajectoryNN2::Traj::Traj(double state_i[], double state_g[], double a_i, const double weights[], const double limits[])
{
	//
	InitialValueGuess2_initialize();
	//
	// this->points.resize(1000, 9);
	double q0[4] = { state_i[0], state_i[1], state_i[2], state_i[3] };
	double q1[4] = { state_g[0], state_g[1], state_g[2], state_g[3] };
	TrajectoryNN2::spiral3(this->r, q0, q1);
	// std::cout << this->r[4] << std::endl;

	if (this->r[4] > 0)
	{
		TrajectoryNN2::path(this->points, this->r, q0, q1);

		TrajectoryNN2::velocity(this->u, state_i[4], a_i, state_g[4], this->r[4]);

		TrajectoryNN2::trajectory(this->points, this->r, this->u);

		this->cost = TrajectoryNN2::traj_eval(this->points, state_g, weights, limits);
	}
	else {
		this->points.resize(1, 9);
		this->cost = inf;
	}
}

TrajectoryNN2::Traj::~Traj()
{
	InitialValueGuess2_terminate();
}

void TrajectoryNN2::Traj::interp(double state[], double time)
{
	int N = this->points.rows() - 1;
	if (time<this->points(0, 0) || time > this->points(N, 0))
	{
		state[0] = state[1] = state[2] = state[3] = state[4] = 0.;
	}
	else
	{
		double s = time*(this->u[0] + time*(this->u[1] / 2. + time*this->u[2] / 3.));
		double delta_s = this->points(1, 1) - this->points(0, 1);
		int n = (int)std::floor(s / delta_s);
		ArrayXXd state_t = (this->points.row(n + 1)*(s - this->points(n, 1)) + this->points.row(n)*(this->points(n + 1, 1) - s)) / delta_s;
		state[0] = state_t(0,2);
		state[1] = state_t(0,3);
		state[2] = state_t(0,4);
		state[3] = state_t(0,5);
		state[4] = state_t(0,7);
	}
}




TrajectoryNN2::Traj TrajectoryNN2::opt_traj(double state_i[], double state_g[], double a_i, const double weights[], const double limits[])
{
	using TrajPtr = std::shared_ptr<TrajectoryNN2::Traj>;
	double cos_theta = std::cos(state_g[2]);
	double sin_theta = std::sin(state_g[2]);
	// double s_step, l_step;
	double state_t[5] = { state_g[0], state_g[1], state_g[2], state_g[3], state_g[4] };
	std::vector<TrajPtr> vec_traj;

	/*
#pragma omp parallel for
	for (int i = 0; i < 5;i++)
	{
		for (int j = 0; j < 5;j++)
		{
			s_step = i - 2.;
			l_step = (j - 2)*0.3;
			state_t[0] = state_g[0] + s_step*cos_theta + l_step*sin_theta;
			state_t[1] = state_g[1] - s_step*sin_theta + l_step*cos_theta;
			vec_traj.push_back(std::make_shared<TrajectoryNN2::Traj>(state_i, state_t, a_i));
		}
	}
	*/
	for (auto s_step : {-2.,-1.,0.,1.,2.})
	{
		for (auto l_step : {-0.6,-0.3,0.,0.3,0.6})
		{
			state_t[0] = state_g[0] + s_step*cos_theta + l_step*sin_theta;
			state_t[1] = state_g[1] - s_step*sin_theta + l_step*cos_theta;
			vec_traj.push_back(std::make_shared<TrajectoryNN2::Traj>(state_i, state_t, a_i, weights, limits));
		}
	}
	

	auto opt_traj_ptr_iter = std::min_element(std::begin(vec_traj), std::end(vec_traj), [](TrajPtr const & t1, TrajPtr const & t2) {return t1->cost < t2->cost; });
	return **opt_traj_ptr_iter;
}

TrajectoryNN2::Planner::Planner()
{
	this->start_time = -1.;
	this->end_time = -1.;

	this->start_planning = false;
	this->end_planning = true;

	this->r[0] = this->r[1] = this->r[2] = this->r[3] = this->r[4] = -1.;
	this->u[0] = this->u[1] = this->u[2] = this->u[3] = -1.;

	this->start_state[0] = this->start_state[1] = this->start_state[2] = this->start_state[3] = this->start_state[4] = -1.;
	this->goal_state[0] = this->goal_state[1] = this->goal_state[2] = this->goal_state[3] = this->goal_state[4] = -1.;

	InitialValueGuess2_initialize();
}

TrajectoryNN2::Planner::~Planner()
{
	InitialValueGuess2_terminate();
}

void TrajectoryNN2::Planner::update(double state_i[], double state_g[], double start_time, double a_i, const double weights[], const double limits[])
{
	this->start_planning = true;
	this->end_planning = false;
	this->start_state[0] = state_i[0];
	this->start_state[1] = state_i[1];
	this->start_state[2] = state_i[2];
	this->start_state[3] = state_i[3];
	this->start_state[4] = state_i[4];
	this->goal_state[0] = state_g[0];
	this->goal_state[1] = state_g[1];
	this->goal_state[2] = state_g[2];
	this->goal_state[3] = state_g[3];
	this->goal_state[4] = state_g[4];
	TrajectoryNN2::spiral3(this->r, state_i, state_g);
	if (this->r[4] > 0)
	{
		TrajectoryNN2::velocity(this->u, state_i[4], a_i, state_g[4], this->r[4]);
		if (u[3] > 0)
		{
			this->start_time = start_time;
			this->end_time = this->start_time + u[3];
		}
	}
	else
	{
		this->start_planning = false;
		this->end_planning = true;
	}
}

void TrajectoryNN2::Planner::reset()
{
	this->start_time = -1.;
	this->end_time = -1.;

	this->start_planning = false;
	this->end_planning = true;

	this->r[0] = this->r[1] = this->r[2] = this->r[3] = this->r[4] = -1.;
	this->u[0] = this->u[1] = this->u[2] = this->u[3] = -1.;

	this->start_state[0] = this->start_state[1] = this->start_state[2] = this->start_state[3] = this->start_state[4] = -1.;
	this->goal_state[0] = this->goal_state[1] = this->goal_state[2] = this->goal_state[3] = this->goal_state[4] = -1.;
}

void TrajectoryNN2::Planner::interp(double time, double state[])
{
	// state - [x,y,theta,k,v]
	state[4] = this->u[0] + this->u[1] * time + this->u[2] * time*time; // v
	double s = this->u[0] * time + this->u[1] * time*time / 2. + this->u[2] * time*time*time / 3.; // s
	VectorXd rr(5);
	rr << this->r[0], this->r[1], this->r[2], this->r[3], this->r[4];
	state[3] = __k(s, rr); // k
	state[2] = this->start_state[2] + __theta(s, rr); // theta
	double x = -1., y = -1.;
	__xy(x, y, s, rr);
	state[0] = this->start_state[0] + x*cos(this->start_state[2]) - y*sin(this->start_state[2]);
	state[1] = this->start_state[1] + x*sin(this->start_state[2]) + y*cos(this->start_state[2]);
}
