#include <iostream>
#include <Eigen/Dense> 
using namespace Eigen;
using namespace std;

void test(VectorXd &vec)
{
	cout << vec << '\n';
	for (int i = 0; i < vec.rows(); i++)
	{
		vec[i] = i / 10.;
	}
	cout << vec << '\n';
}

void test(MatrixXd &mat)
{
	cout << mat << '\n';
	for (int i = 0; i < mat.rows(); i++)
		for (int j = 0; j < mat.cols(); j++)
			mat(i, j) = (i + j) / 10.;
	cout << mat << '\n';
	mat << 0.7, 0.8, 0.9,1.,2.,3.;
	cout << mat << '\n';
}

int main()
{
	VectorXd q(4);
	q << 1, 2, 3, 4;
	ArrayXXd t;
	t.resize(1, 5);
	t(0, 0) = 0;
	t(0, 1) = q(0);
	t(0, 2) = q(1);
	t(0, 3) = q(2);
	t(0, 4) = q(3);
	cout << t;


	/*
	VectorXd vec(4);
	vec << 1 , 2 , 3 , 4;
	test(vec);

	MatrixXd mat(2, 3);
	mat << 1, 2, 3, 4, 5, 6;
	test(mat);
	//
	VectorXd v2 = VectorXd::LinSpaced(10, 0., 12.);
	test(v2);
	v2 *= 2.;
	cout << v2 << '\n';
	
	Matrix3f A;
	Vector3f b;
	A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
	b << 3, 3, 4;
	cout << "Here is the matrix A:\n" << A << endl;
	cout << "Here is the vector b:\n" << b << endl;
	Vector3f x = A.colPivHouseholderQr().solve(b);
	cout << "The solution is:\n" << x << endl;
	cout << b.array().cos() << endl;
	*/

	//
	return 0;
}