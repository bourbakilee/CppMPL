#include <spiral3.h>
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
	VectorXd p(5), r(5);
	VectorXd q0(4);
	q0 << 0., 0., 0., 0.01;
	VectorXd q1(4);
	q1 << 100., 40., pi / 6, -0.01;
	spiral3::calc_path(p, r, q0, q1, db);
	std::cout << p << '\n';
	std::cout << r << '\n';

	sqlite3_close(db);
	return 0;
}