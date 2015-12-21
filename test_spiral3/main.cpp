#include <spiral3.h>
#include <trajectory.h>
#include <iostream>
#include <string>
// using namespace std;

int main()
{
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		std::cout << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
		return 0;
	}
	else {
		std::cout << "Opened database successfully" << std::endl;
	}
	//
	VectorXd r(5);
	VectorXd q0(4);
	q0 << 20., 20., -pi/6, -0.01;
	VectorXd q1(4);
	// q1 << 100., 40., pi / 6, -0.01;
	q1 << 70., 30., pi / 6, -0.01;
	spiral3::spiral3( r, q0, q1, db);
	std::cout << r << '\n';
	//
	ArrayXXd points(10,5);
	spiral3::path(points, r, q0, q1);
	//std::cout << points;
	//traj
	double u[4];
	double v0=7., a0=0.5, vg=12.;
	trajectory::velocity(u, v0, a0, vg, r[4]);
	double r1[] = { r[0],r[1],r[2],r[3],r[4] };
	ArrayXXd traj(10, 5);
	trajectory::trajectory(traj, points, r1, u);
	std::cout << traj;
	//
	sqlite3_close(db);
	return 0;
}