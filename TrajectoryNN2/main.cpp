#include "TrajectoryNN2.h"
#include<iostream>
#include <chrono>
#include <random>
using namespace std;
using namespace TrajectoryNN2;

const double pi = 3.141592653589793;

int test_1()
{
	double state_i[5] = { 0., 0., 0., 0.,1. };
	double state_g[5] = { 0., 10., pi, 0,1. };
	double a_i = 0.5;


	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	Traj traj = opt_traj(state_i, state_g, a_i);
	
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_s = end - start;
	std::cout << "elapsed time: " << elapsed_s.count() << "s\n";


	cout << "r:\n" << traj.r[0] << "\t" << traj.r[1] << "\t" << traj.r[2] << "\t" << traj.r[3] << "\t" << traj.r[4] << endl;
	cout << "u:\n" << traj.u[0] << "\t" << traj.u[1] << "\t" << traj.u[2] << "\t" << traj.u[3] << endl;
	cout << "cost:\n" << traj.cost << endl;
	cout << "traj:\n" << traj.points << endl;

	Traj traj1(state_i, state_g, a_i);
	cout << "r:\n" << traj1.r[0] << "\t" << traj1.r[1] << "\t" << traj1.r[2] << "\t" << traj1.r[3] << "\t" << traj1.r[4] << endl;
	cout << "u:\n" << traj1.u[0] << "\t" << traj1.u[1] << "\t" << traj1.u[2] << "\t" << traj1.u[3] << endl;
	cout << "cost:\n" << traj1.cost << endl;
	cout << "traj:\n" << traj1.points << endl;

	return 0;
}

int test_2()
{
	InitialValueGuess2_initialize();
	
	int success = 0;
	double r[5], u[4];
	double state_i[5] = { 0., 0., 0., 0.01, 2. };
	double a_i = 0.2;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (double x : {10., 20., 30., 40., 50.})
	{
		for (double y : {-20., -10., -5., 0., 5., 10., 20.})
		{
			for (double theta : {-6.*pi / 18, -5.*pi / 18, -4.*pi / 18, -3.*pi / 18, -2.*pi / 18, -1.*pi / 18, 0., 1.*pi / 18, 2.*pi / 18, 3.*pi / 18, 4.*pi / 18, 5.*pi / 18, 6.*pi / 18, })
			{
				for (double k : {-0.02, -0.015, -0.01, -0.005, 0., 0.005, 0.01, 0.015, 0.02})
				{
					double state_g[5] = { x,y,theta,k,3. };
					double q0[4] = { state_i[0], state_i[1], state_i[2], state_i[3] };
					double q1[4] = { state_g[0], state_g[1], state_g[2], state_g[3] };
					TrajectoryNN2::spiral3(r, q0, q1);
					if (r[4] > 0)
					{
						
						TrajectoryNN2::velocity(u, state_i[4], a_i, state_g[4], r[4]);
						if (u[3] > 0)
						{
							success++;
						}
					}					
				}
			}
		}
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_s = end - start;
	std::cout << "elapsed time: " << elapsed_s.count() << "s\n";
	std::cout << "average time: " << elapsed_s.count() / (5 * 7 * 13 * 9) << "s\n";
	std::cout << "Success: " << success << "\n";
	std::cout << "Success Ratio: " << success / (5 * 7 * 13 * 9) << "\n";

	InitialValueGuess2_terminate();
	return 0;
}

int test_3()
{
	InitialValueGuess2_initialize();

	int success = 0;
	int path_success = 0;
	double r[5], u[4];

	std::random_device rd;
	std::mt19937 gen(rd());
	//std::uniform_int_distribution<> dis_ki(-20, 20), dis_vi(0, 800), dis_xg(10, 500), dis_yg(-500, 500), dis_tg(-90, 90), dis_kg(-20, 20), dis_vg(0, 800), dis_ai(-300, 150);
	std::uniform_int_distribution<> dis_ki(-5, 5), dis_vi(0, 800), dis_xg(10, 500), dis_yg(-80, 80), dis_tg(-30, 30), dis_kg(-5, 5), dis_vg(0, 800), dis_ai(-300, 150);


	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (int i = 0; i < 100000;i++)
	{
		double q0[4] = { 0.,0.,0.,dis_ki(gen) / 100. };
		double q1[4] = { dis_xg(gen) / 10., dis_yg(gen) / 10., dis_tg(gen) * pi / 180. };
		double v_i = dis_vi(gen) / 100.;
		double a_i = dis_ai(gen) / 100.;
		double v_g = dis_vg(gen) / 100.;

		TrajectoryNN2::spiral3(r, q0, q1);
		if (r[4] > 0)
		{
			path_success++;
			TrajectoryNN2::velocity(u, v_i, a_i, v_g, r[4]);
			if (u[3] > 0)
			{
				success++;
			}
		}
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_s = end - start;
	std::cout << "elapsed time: " << elapsed_s.count() << "s\n";
	std::cout << "average time: " << elapsed_s.count() / 100000. << "s\n";
	std::cout << "Success: " << success << "\n";
	std::cout << "Success Ratio: " << success / 100000. << "\n";
	std::cout << "Path Success: " << path_success << "\n";
	std::cout << "Path Success Ratio: " << path_success / 100000. << "\n";

	InitialValueGuess2_terminate();
	return 0;
}


int main()
{
	//test_1();
	test_3();
}