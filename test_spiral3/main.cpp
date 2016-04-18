//#define TEST_WITHOUT_SQLITE3
#define TEST_WITH_SQLITE3
#include <iostream>
#include <string>
#include <chrono>
#include <random>
// using namespace std;

#ifdef TEST_WITH_SQLITE3
#ifndef WITH_SQLITE3
#define WITH_SQLITE3
#include <spiral3.h>
#include <trajectory.h>
//#include <chrono>
using namespace Eigen;
void test_with_sqlite3()
{
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		std::cout << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
		return;
	}
	else {
		std::cout << "Opened database successfully" << std::endl;
	}
	//
	VectorXd r(5);
	VectorXd q0(4);
	q0 << 10., 50., 0., 0.;
	VectorXd q1(4);
	// q1 << 100., 40., pi / 6, -0.01;
	q1 << 20, 51, 0, 0;

	double u[4];
	double v0 = 5., a0 = 0.5, vg = 6.;

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	spiral3::spiral3( r, q0, q1, db);
	
	//
	ArrayXXd points(10,9);


	spiral3::path(points, r, q0, q1);
	//std::cout << points;
	//traj
	
	trajectory::velocity(u, v0, a0, vg, r[4]);
	
	// double r1[] = { r[0],r[1],r[2],r[3],r[4] };
	//ArrayXXd traj(10, 5);
	trajectory::trajectory(points, r.data(), u);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_s = end - start;
	std::cout << "elapsed time: " << elapsed_s.count() << "s\n";

	std::cout <<"r: \n"<< r << '\n';
	std::cout << "u: \n"<<u[0] << std::endl << u[1] << std::endl << u[2] << std::endl << u[3] << std::endl;
	std::cout << points;
	//
	sqlite3_close(db);
}


void test_3()
{
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		std::cout << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
		return;
	}
	else {
		std::cout << "Opened database successfully" << std::endl;
	}

	int success = 0;
	int path_success = 0;
	double r[5], u[4];

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis_ki(-20, 20), dis_vi(0, 800), dis_xg(10, 500), dis_yg(-500, 500), dis_tg(-90, 90), dis_kg(-20, 20), dis_vg(0, 800), dis_ai(-300, 150);
	//std::uniform_int_distribution<> dis_ki(-10, 10), dis_vi(0, 800), dis_xg(10, 500), dis_yg(-100, 100), dis_tg(-45, 45), dis_kg(-10, 10), dis_vg(0, 800), dis_ai(-300, 150);


	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (int i = 0; i < 100000; i++)
	{
		double q0[4] = { 0.,0.,0.,dis_ki(gen) / 100. };
		double q1[4] = { dis_xg(gen) / 10., dis_yg(gen) / 10., dis_tg(gen) * pi / 180. };
		double v_i = dis_vi(gen) / 100.;
		double a_i = dis_ai(gen) / 100.;
		double v_g = dis_vg(gen) / 100.;

		spiral3::spiral3(r, q0, q1, db);
		if (r[4] > 0)
		{
			path_success++;
			trajectory::velocity(u, v_i, a_i, v_g, r[4]);
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

	sqlite3_close(db);
}

#endif
#endif

#ifdef TEST_WITHOUT_SQLITE3
#ifdef WITH_SQLITE3
#undef WITH_SQLITE3
#endif
#include <trajectory.h>

void test_no_sql()
{
	VectorXd r(5);
	VectorXd q0(4);
	q0 << 20., 20., -pi / 6, -0.01;
	VectorXd q1(4);
	// q1 << 100., 40., pi / 6, -0.01;
	q1 << 70., 30., pi / 6, -0.01;
	spiral3::spiral3(r, q0, q1);
	std::cout << r << '\n';
	//
	ArrayXXd points(1, 9);
	spiral3::path(points, r, q0, q1);

	double u[4];
	double v0 = 7., a0 = 0.5, vg = 12.;
	trajectory::velocity(u, v0, a0, vg, r[4]);
	std::cout << u[0] << std::endl << u[1] << std::endl << u[2] << std::endl << u[3] << std::endl;
	double r1[] = { r[0],r[1],r[2],r[3],r[4] };
	trajectory::trajectory(points, r1, u);
	std::cout << points;
}


#endif





int main()
{
#ifdef TEST_WITH_SQLITE3
	//test_with_sqlite3();
	test_3();
#endif

#ifdef TEST_WITHOUT_SQLITE3
	test_no_sql();
#endif

	return 0;
}