//#define TEST_WITHOUT_SQLITE3
#define TEST_WITH_SQLITE3
#include <iostream>
#include <string>
// using namespace std;

#ifdef TEST_WITH_SQLITE3
#ifndef WITH_SQLITE3
#define WITH_SQLITE3
#include <spiral3.h>
#include <trajectory.h>
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
	spiral3::spiral3( r, q0, q1, db);
	std::cout << r << '\n';
	//
	ArrayXXd points(10,9);
	spiral3::path(points, r, q0, q1);
	//std::cout << points;
	//traj
	double u[4];
	double v0=5., a0=0.5, vg=6.;
	trajectory::velocity(u, v0, a0, vg, r[4]);
	std::cout << u[0] << std::endl << u[1] << std::endl << u[2] << std::endl << u[3] << std::endl;
	double r1[] = { r[0],r[1],r[2],r[3],r[4] };
	//ArrayXXd traj(10, 5);
	trajectory::trajectory(points, r1, u);
	std::cout << points;
	//
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
	test_with_sqlite3();
#endif

#ifdef TEST_WITHOUT_SQLITE3
	test_no_sql();
#endif

	return 0;
}