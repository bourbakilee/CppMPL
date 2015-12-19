#include "spiral3.h"
#include <cmath>
#include <algorithm>

inline int sgn(double t)
{
	return (t > 0) - (t < 0);
}

// thetas与ss大小一致
inline void spiral3::__thetas(vec& thetas, vec& ss, vec& r)
{
	for (int i = 0; i < ss.n_rows; i++)
	{
		thetas[i] = spiral3::__theta(ss[i], r);
	}
}

inline void spiral3::__xy(double& x, double& y, double s, vec& r, double ref_size = 8.0)
{
	int N = std::ceil(s / ref_size);
	vec s_list = linspace(0., s, N + 1);
	double h = s / (8 * N);
	vec thetas(9), fs(9), gs(9);
	x = y = 0.;
	for (int i = 0; i < N; i++) {
		vec ss = linspace(s_list[i], s_list[i + 1], 9);
		spiral3::__thetas(thetas, ss, r);
		fs = cos(thetas);
		gs = sin(thetas);
		x += h / 3 * (fs[0] + 2 * (fs[2] + fs[4] + fs[6]) + 4 * (fs[1] + fs[3] + fs[5] + fs[7]) + fs[8]);
		y += h / 3 * (gs[0] + 2 * (gs[2] + gs[4] + gs[6]) + 4 * (gs[1] + gs[3] + gs[5] + gs[7]) + gs[8]);
	}
}

/*
inline double spiral3::__y(double s, double r[], double ref_step_size = 8.0)
{

}
*/

void spiral3::__jacobian(mat& Jcb, vec& p, vec& r)
{
	vec ss = linspace(0., p[4], 9), thetas(9);
	spiral3::__thetas(thetas, ss, r);
	vec cos_t = cos(thetas);
	vec sin_t = sin(thetas);
	vec c_p_1{ 0, 1851 / 8192, 363 / 1024, 9963 / 8192, 51 / 64, 14475 / 8192, 891 / 1024, 13083 / 8192, 3 / 8 };
	c_p_1 *= p[4] * p[4] / 24;
	vec c_p_2{ 0, 795 / 8192, 123 / 1024, 2187 / 8192, 3 / 64, -2325 / 8192, -405 / 1024, -10437 / 8192, -3 / 8 };
	c_p_2 *= p[4] * p[4] / 24;
	vec c_s_1{ 0, (2871 * p[0] + 1851 * p[1] - 795 * p[2] + 169 * p[3]) / 8192, 
		(247 * p[0] + 363 * p[1] - 123 * p[2] + 25 * p[3]) / 1024, 
		(4071 * p[0] + 9963 * p[1] - 2187 * p[2] + 441 * p[3]) / 8192, 
		(15 * p[0] + 51 * p[1] - 3 * p[2] + p[3]) / 64, 
		(3655 * p[0] + 14475 * p[1] + 2325 * p[2] + 25 * p[3]) / 8192, 
		(231 * p[0] + 891 * p[1] + 405 * p[2] + 9 * p[3]) / 1024, 
		(3927 * p[0] + 13083 * p[1] + 10437 * p[2] + 1225 * p[3]) / 8192, 
		(p[0] + 3 * p[1] + 3 * p[2] + p[3]) / 8 };
	c_s_1 *=  p[4] / 24;
	vec c_s_2{ 1 / 24, 1 / 6, 1 / 12, 1 / 6, 1 / 12, 1 / 6, 1 / 12, 1 / 6, 1 / 24 };
	Jcb[0, 0] = -dot(c_p_1, sin_t);
	Jcb[0, 1] = dot(c_p_2, sin_t);
	Jcb[0, 2] = dot(c_s_2, cos_t) - dot(c_s_1, sin_t);
	Jcb[1, 0] = dot(c_p_1, cos_t);
	Jcb[1, 1] = -dot(c_p_2, cos_t);
	Jcb[1, 2] = dot(c_s_2, sin_t) + dot(c_s_1, cos_t);
	Jcb[2, 0] = 0.375 / p[4];
	Jcb[2, 1] = Jcb[2, 0];
	Jcb[2, 2] = c_s_1[8];
}

// bd_con: boudary conditions [k0,x1,y1,theta1,k1]
// pp: initial value of parameters [p1,p2,sg]
void spiral3::optimize(vec& pp, vec& bd_con, double k_m = 0.2)
{
	mat Jcb(3, 3);
	if (pp[2] <= 0)
	{
		pp[0] = (2 * bd_con[0] + bd_con[4]) / 3;
		pp[1] = (bd_con[0] + 2 * bd_con[4]) / 3;
		pp[2] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::min(std::abs(bd_con[3]), 2 * pi - std::abs(bd_con[3]));
	}
	vec q_g{ bd_con[1], bd_con[2],bd_con[3] }; // x1,y1,theta1
	vec p{ bd_con[0], pp[0], pp[1], bd_con[4], pp[2] };
	vec r{ spiral3::__a(p),spiral3::__b(p), spiral3::__c(p), spiral3::__d(p), p[4] };
	vec dq{ 1.,1.,1. };
	int times = 0;
	double x_p = 0., y_p = 0., theta_p = 0.;
	while (std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2 || times < iter_num) {
		times += 1;
		spiral3::__jacobian(Jcb, p, r);
		theta_p = spiral3::__theta(p[4], r);
		spiral3::__xy(x_p, y_p, p[4], r);
		dq = q_g - vec({ x_p,y_p,theta_p });
		pp += solve(Jcb, dq);
		pp[0] = p[1] = sgn(pp[0])*std::min(k_m, std::abs(pp[0]));
		pp[1] = p[2] = sgn(pp[1])*std::min(k_m, std::abs(pp[1]));
		pp[2] = r[4] = p[4] = std::max(pp[2], 1.);
		r[0] = spiral3::__a(p);
		r[1] = spiral3::__b(p);
		r[2] = spiral3::__c(p);
		r[3] = spiral3::__d(p);
	}
	if (std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2)
	{
		pp[0] = pp[1] = 0.;
		pp[2] = -1.;
	}
}

void spiral3::select_init_val(vec& pp, vec& bd_con, sqlite3* db)
{

}

void spiral3::calc_path(vec& p, vec& r, vec& q0, vec& q1, sqlite3* db, double k_m = 0.2)
{
	
}