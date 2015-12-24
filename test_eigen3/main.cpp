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

void test(double t[])
{
	t[0] = 0.;
	t[1] = 0.;
}

int main()
{
	ArrayXXd t1(2, 3);
	t1 << 1.1, 2.2, 3.3, 4.4, 5.5, 6.6;
	Array<unsigned int, Dynamic, 3> t2 = (t1 / 0.1).cast<unsigned int>();
	cout << t2 << endl;

	/*
	ArrayXXd t1(3, 2);
	t1 << 1, 2, 3, 4, 5, 6;
	cout << t1 << endl;
	ArrayXXd t2(3, 2);
	t2 << 7, 8, 9, 0, 1, 2;
	cout << t2 << endl;
	cout << t1.block(0, 0, 1, 2).transpose() + t2.block(0, 0, 2, 1) << endl;

	double *q = t1.data();
	cout << *q << endl;
	q[0] = 90., q[2] = 80; // matrix and array elements are in colum stored in memory
	cout << *q << endl;
	cout << t1 << endl;
	*/

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