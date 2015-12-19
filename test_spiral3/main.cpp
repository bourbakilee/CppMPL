#include <spiral3.h>
#include <iostream>
#include <string>
using namespace std;

int main()
{
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		cout << "Can't open database: " << sqlite3_errmsg(db) << endl;
		sqlite3_close(db);
		return 0;
	}
	else {
		cout << "Opened database successfully" << endl;
	}
	//
	arma::vec p(5), r(5);
	arma::vec q0{ 0.,0.,0.,0.01 };
	arma::vec q1{ 100.,40.,pi / 6,-0.01 };
	spiral3::calc_path(p, r, q0, q1, db);
	cout << p << '\n';
	cout << r << '\n';

	sqlite3_close(db);
	return 0;
}