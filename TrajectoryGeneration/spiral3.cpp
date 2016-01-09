#include "spiral3.h"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace Eigen;
/*
inline int sgn(double t)
{
	return (t > 0) - (t < 0);
}
*/

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
inline void spiral3::__thetas(VectorXd& thetas, VectorXd& ss, VectorXd& r)
{
	for (int i = 0; i < ss.rows(); i++)
	{
		thetas[i] = spiral3::__theta(ss[i], r);
	}
}

inline void spiral3::__xy(double& x, double& y, double s, VectorXd& r, double ref_size)
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
		spiral3::__thetas(thetas, ss, r);
		fs = thetas.array().cos();
		gs = thetas.array().sin();
		x += h / 3 * (fs[0] + 2 * (fs[2] + fs[4] + fs[6]) + 4 * (fs[1] + fs[3] + fs[5] + fs[7]) + fs[8]);
		y += h / 3 * (gs[0] + 2 * (gs[2] + gs[4] + gs[6]) + 4 * (gs[1] + gs[3] + gs[5] + gs[7]) + gs[8]);
	}
}

/*
inline double spiral3::__y(double s, double r[], double ref_step_size = 8.0)
{

}
*/

void spiral3::__jacobian(Matrix3d& Jcb, VectorXd& p, VectorXd& r)
{
	VectorXd ss = VectorXd::LinSpaced(9,0., p[4]), thetas(9);
	spiral3::__thetas(thetas, ss, r);
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
	c_s_1 *=  p[4] / 24;
	VectorXd c_s_2(9);
	c_s_2 << 1 / 24., 1 / 6., 1 / 12., 1 / 6., 1 / 12., 1 / 6., 1 / 12., 1 / 6., 1 / 24.;
	Jcb << -c_p_1.dot(sin_t), c_p_2.dot(sin_t), c_s_2.dot(cos_t) - c_s_1.dot(sin_t),
		c_p_1.dot(cos_t), -c_p_2.dot(cos_t), c_s_2.dot(sin_t) + c_s_1.dot(cos_t),
		0.375 * p[4], 0.375 * p[4], c_s_1[8];
	//std::cout << Jcb << '\n';
	/*
	Jcb[0, 0] = -arma::dot(c_p_1, sin_t);
	Jcb[0, 1] = arma::dot(c_p_2, sin_t);
	Jcb[0, 2] = arma::dot(c_s_2, cos_t) - arma::dot(c_s_1, sin_t);
	Jcb[1, 0] = arma::dot(c_p_1, cos_t);
	Jcb[1, 1] = -arma::dot(c_p_2, cos_t);
	Jcb[1, 2] = arma::dot(c_s_2, sin_t) + arma::dot(c_s_1, cos_t);
	Jcb[2, 0] = 0.375 / p[4];
	Jcb[2, 1] = Jcb[2, 0];
	Jcb[2, 2] = c_s_1[8];
	*/
}

// bd_con: boudary conditions [k0,x1,y1,theta1,k1]
// pp: initial value of parameters [p1,p2,sg]
void spiral3::optimize(VectorXd& p, VectorXd& r, VectorXd& bd_con, int iter_num, double k_m)
{
	Matrix3d Jcb;
	VectorXd pp(3);
	pp << p[1], p[2], p[4];
	if (pp[2] <= 0)
	{
		pp[0] = (2 * bd_con[0] + bd_con[4]) / 3.;
		pp[1] = (bd_con[0] + 2 * bd_con[4]) / 3.;
		// pp[2] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::min(std::abs(bd_con[3]), two_pi - std::abs(bd_con[3]));
		pp[2] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::abs(bd_con[3]);
	}
	pp[0] = std::max(std::min(pp[0], k_m), -k_m);
	pp[1] = std::max(std::min(pp[1], k_m), -k_m);
	pp[2] = std::max(std::min(pp[2], 1000.), 1.);
	VectorXd q_g(3), q_p(3);
	q_g << bd_con[1], bd_con[2], bd_con[3]; // x1,y1,theta1
	// VectorXd p{ bd_con[0], pp[0], pp[1], bd_con[4], pp[2] };
	// VectorXd r{ spiral3::__a(p),spiral3::__b(p), spiral3::__c(p), spiral3::__d(p), p[4] };
	VectorXd dq(3);
	dq << 1., 1., 1.;
	int times = 0;
	double x_p = 0., y_p = 0., theta_p = 0.;
	while( (std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2 )&&( times < iter_num )) {
		times += 1;
		//std::cout << "iteration times: " << times << std::endl;
		spiral3::__jacobian(Jcb, p, r);
		theta_p = spiral3::__theta(p[4], r);
		spiral3::__xy(x_p, y_p, p[4], r);
		q_p << x_p, y_p, mod2pi(theta_p);
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		dq = q_g - q_p;
		dq[2] = mod2pi(dq[2]);
		//
		pp += Jcb.colPivHouseholderQr().solve(dq);
		//pp[0] = p[1] = sgn(pp[0])*std::min(k_m, std::abs(pp[0]));
		pp[0] = p[1] = std::max(std::min(pp[0], k_m), -k_m);
		//pp[1] = p[2] = sgn(pp[1])*std::min(k_m, std::abs(pp[1]));
		pp[1] = p[2] = std::max(std::min(pp[1], k_m), -k_m);
		pp[2] = r[4] = p[4] = std::max(std::min(pp[2],1000.), 1.);
		r[0] = spiral3::__a(p);
		r[1] = spiral3::__b(p);
		r[2] = spiral3::__c(p);
		r[3] = spiral3::__d(p);
	}
	if (std::abs(dq[0]) > eps1 || std::abs(dq[1]) > eps1 || std::abs(dq[2]) > eps2)
	{
		r[4] = -1.; //通过检查r[4]来判断结果是否有效
	}
}

#ifdef WITH_SQLITE3

//callback function for sqlite3_exec
static int callback(void *data, int argc, char **argv, char **azColName)
{
	double *pp = (double*)data;
	pp[0] = strtod(argv[0], nullptr);
	pp[1] = strtod(argv[1], nullptr);
	pp[2] = strtod(argv[2], nullptr);
	return 0;
}

//sql_construct
inline std::string sql_construct(int i, int j, int k, int l, int m)
{
	std::stringstream sql;
	sql << "select p1,p2,sg from InitialGuessTable where k0=" << i << " and x1=" << j << " and y1=" << k << " and theta1=" << l << " and k1=" << m;
	return sql.str();
}

void spiral3::select_init_val(VectorXd& p, VectorXd& r, VectorXd& bd_con, sqlite3* db)
{
	int i = std::llround(bd_con[0] * 40);
	int j = std::llround(bd_con[1] * 16 / 49);
	int k = std::llround(bd_con[2] * 0.16);
	int l = std::llround(bd_con[3] * 16 / pi);
	int m = std::llround(bd_con[4] * 40);
	//
	p[0] = bd_con[0];
	p[3] = bd_con[4];
	//
	char* errMsg = nullptr;
	double pp[3]; //p1,p2,sg
	std::string sql = sql_construct(i, j, k, l, m);
	int rc = sqlite3_exec(db, sql.c_str(), callback, pp, &errMsg);
	if (rc == SQLITE_OK && pp[2]>0)
	{
		p[1] = pp[0];
		p[2] = pp[1];
		p[4] = r[4] = pp[2];
		r[0] = spiral3::__a(p);
		r[1] = spiral3::__b(p);
		r[2] = spiral3::__c(p);
		r[3] = spiral3::__d(p);
		//std::cout << "Initial values: " << std::endl;
		//std::cout << p << '\n';
		//std::cout << r << '\n';
	}
	else if(rc == SQLITE_OK && pp[2]<0)
	{
		p[1] = p[2] = r[0] = r[1] = r[2] = r[3] = 0.;
		p[4] = r[4] = -1.;
	}
	else
	{
		std::cout << errMsg;
	}
}

// q - (x,y,theta,k)
void spiral3::spiral3(VectorXd& r, VectorXd& q0, VectorXd& q1, sqlite3* db, double k_m)
{
	VectorXd p(5);
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
	if (std::abs(theta_r) > pi / 2 || bd_con[1] < 0 || bd_con[1]>50 || std::abs(bd_con[2])>50 || std::abs(bd_con[0])>0.2 || std::abs(bd_con[4])>0.2)
	{
		p[1] = (2 * bd_con[0] + bd_con[4]) / 3;
		p[2] = (bd_con[0] + 2 * bd_con[4]) / 3;
		p[4] = std::sqrt(bd_con[1] * bd_con[1] + bd_con[2] * bd_con[2]) + 10. * std::min(std::abs(bd_con[3]), two_pi - std::abs(bd_con[3]));
		p[0] = q0[3];
		p[3] = q1[3];
		if (p[4] > 0) {
			r[0] = spiral3::__a(p);
			r[1] = spiral3::__b(p);
			r[2] = spiral3::__c(p);
			r[3] = spiral3::__d(p);
		}
		else {
			r[0] = r[1] = r[2] = r[3] = 0.;
		}
		r[4] = p[4];
	}
	else
	{
		spiral3::select_init_val(p, r, bd_con, db);
	}
	
	// optimize p and r
	spiral3::optimize(p, r, bd_con);
}

void spiral3::spiral3(double r[], double q0[], double q1[], sqlite3* db, double k_m)
{
	VectorXd rr(5), qs(4), qg(4);
	//rr << r[0], r[1], r[2], r[3], r[4];
	qs << q0[0], q0[1], q0[2], q0[3];
	qg << q1[0], q1[1], q1[2], q1[3];
	spiral3::spiral3(rr, qs, qg, db, k_m);
	r[0] = rr[0];
	r[1] = rr[1];
	r[2] = rr[2];
	r[3] = rr[3];
	r[4] = rr[4];
}

#else

void spiral3::spiral3(VectorXd& r, VectorXd& q0, VectorXd& q1, double k_m)
{
	// p <-> r
	VectorXd p(5);
	p << 0., 0., 0., 0., -1.;
	//bd_con
	double cc = cos(q0[2]);
	double ss = sin(q0[2]);
	double theta_r = mod2pi(q1[2] - q0[2]);
	VectorXd bd_con(5);
	bd_con << q0[3], (q1[0] - q0[0])*cc + (q1[1] - q0[1])*ss, -(q1[0] - q0[0])*ss + (q1[1] - q0[1])*cc, theta_r, q1[3];

	spiral3::optimize(p, r, bd_con);
}

void spiral3::spiral3(double r[], double q0[], double q1[], double k_m)
{
	VectorXd rr(5), qs(4), qg(4);
	qs << q0[0], q0[1], q0[2], q0[3];
	qg << q1[0], q1[1], q1[2], q1[3];
	spiral3::spiral3(rr, qs, qg, k_m);
	r[0] = rr[0];
	r[1] = rr[1];
	r[2] = rr[2];
	r[3] = rr[3];
	r[4] = rr[4];
}

#endif



// points - {t, s, x, y, theta, k, dk, v, a}, here only 1-5 cols are used
void spiral3::path(ArrayXXd& points, VectorXd& r, VectorXd& q0, VectorXd& q1, double length, double ref_size)
{
	if (length < 0) { length = r[4]; }
	int N = (int)std::ceil(length / ref_size);
	double delta = length / N;
	points.resize(N+1, 9); // t, s, x, y, theta, k, dk, v, a
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