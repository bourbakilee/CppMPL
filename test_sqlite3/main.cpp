#include<sqlite3.h>
#include<iostream>
#include<cstdlib>
#include<string>
#include <sstream>
using namespace std;

static int callback(void *data, int argc, char **argv, char **azColName)
{
	double *pp = (double*)data;
	pp[0] = strtod(argv[0], nullptr);
	pp[1] = strtod(argv[1], nullptr);
	pp[2] = strtod(argv[2], nullptr);
	return 0;
}

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
	char* errMsg;
	char* sql = "select p1,p2,sg from InitialGuessTable where k0=0 and x1=10 and y1=2 and theta1=3 and k1=-2";
	double pp[3];
	rc = sqlite3_exec(db, sql, callback, &pp, &errMsg);
	cout << pp[0] << endl << pp[1] << endl << pp[2] << endl;
	//
	sqlite3_close(db);

	stringstream test;
	test << "test" << 1;
	cout << test.str() << endl;
	return 0;
}