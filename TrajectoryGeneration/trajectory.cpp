#include "trajectory.h"
#include <cmath>

namespace trajectory {
	// u: u0,u1,u2. tg
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a,j,al)]
	void velocity(double u[], double v0, double a0, double vg, double sg) 
	{
		/*
		# v(t) = u0 + u1*t + u2*t**2
		# return: q0~q2, tg
		u0, u1, u2, tg = None, None, None, None
		delta = (2*v0+vg)**2 + 6*a0*sg
		if delta >= 0:
			u0 = v0
			u1 = a0
			tg = 3*sg/(2*v0+vg) if np.abs(a0)<1.e-6 else (np.sqrt(delta)-2*v0-vg)/a0
			u2 = (vg - v0 - a0*tg)/tg**2
		return (u0,u1,u2,tg)
		*/
		u[0] = v0;
		u[1] = a0;
		double delta = (2 * v0 + vg)*(2 * v0 + vg) + 6 * a0*sg;
		if (delta > eps2)
		{
			if (std::abs(a0) < eps2)
				u[3] = 3 * sg / (2 * v0 + vg);
			else
				u[3] = (std::sqrt(delta) - 2 * v0 - vg) / a0;
		}
		else if (std::abs(delta) < eps2)
		{
			u[3] = (-2 * v0 - vg) / a0;
		}
		else
		{
			u[3] = -1.;
		}
		u[2] = (vg - v0 - a0*u[3]) / (u[3] * u[3]);
	}

	
	// TODO: considering add road and obstacle(static and dynamic) items - 2015.12.22
	// array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	// r: a, b, c, d, sg
	// u: u0,u1,u2. tg
	// parameter traj must be pre-computed by spiral3::path procedure, and colums 1-5 must be filled
	void trajectory(ArrayXXd& traj, double r[],double u[], double ref_length, double ref_time)
	{
		int N = traj.rows();
		// traj.resize(N, 11);
		// traj.block(0, 1, N, 5) = points; // 1-5 col

		//0 col: t
		ArrayXXd times = VectorXd::LinSpaced(N, 0., u[3]);
		ArrayXXd lengths = times*(u[0] + times*(u[1] / 2. + times*u[2] / 3.));
	
		traj(0, 0) = 0., traj(N - 1, 0) = u[3];
		int j = 1;
		for (int i = 1; i < N-1; i++)
		{
			while(traj(i, 1) >= lengths(j, 0))
				j++;
			traj(i, 0) = times(j - 1, 0) + (traj(i, 1) - lengths(j - 1, 0))*(times(j, 0) - times(j - 1, 0)) / (lengths(j, 0) - lengths(j - 1, 0));
		}
		// 7 col: v
		traj.col(7) = u[0] + traj.col(0)*(u[1] + traj.col(0)*u[2]);
		// 6 col: dk/dt = v*dk/ds
		traj.col(6) = traj.col(7)*(r[1] + traj.col(1)*(2 * r[2] + 3 * r[3] * traj.col(1)));
		// 8 col: a
		traj.col(8) = u[1] + 2 * u[2] * traj.col(0);
		/*
		// 9 col: j
		//traj.col(9) = 2 * u[2]*ArrayXXd::Ones(N,1);
	    // 10 col: al
		//traj.col(10) = traj.col(7)*traj.col(7)*traj.col(5);
		*/
		// ref_time -> absolute time
		traj.col(0) += ref_time;
		// ref_length -> absolute length
		traj.col(1) += ref_length;
	}

	// Eigen3 <-> OpenCV
	// I wouldn't necessarily recommend to use cv2eigen() and eigen2cv() though
	// use Eigen::Map to just map the memory (no cost in copying anything) and cv::Mat(void*, ...) to map the data back
	// MatrixXd matrix; double* arrayd = matrix.data();


	// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	double eval_traj(ArrayXXd& traj, double weights[])
	{
		int N = traj.rows();
		ArrayXXd cost_matrix(N,7); // k, dk, v, a, a_c, off_set(x,y), env(t,x,y,theta)
		cost_matrix.col(0) = weights[0] * traj.col(5).abs(); // |k|
		cost_matrix.col(1) = weights[1] * traj.col(6).abs(); //|dk|
		cost_matrix.col(2) = weights[2] * traj.col(7)*traj.col(7); // v^2
		cost_matrix.col(3) = weights[3] * traj.col(8)*traj.col(8); // a^2
		cost_matrix.col(4) = weights[4] / weights[2] * traj.col(5).abs()*cost_matrix.col(2); // v^2*k
	}
}

/*
// Eigen3 <-> OpenCV : Example
//
#define EIGEN_RUNTIME_NO_MALLOC // Define this symbol to enable runtime tests for allocations
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>
#include "opencv2/core/eigen.hpp"
#include "opencv2/opencv.hpp"
using namespace Eigen;
using namespace cv;
using namespace std;

void EnergyFilter(Mat& src,Mat& dst,double alpha)
{
int n_pixels=src.rows*src.cols;
// Image to row-vector
Mat m=src.reshape(1,n_pixels).clone();
// To double
m.convertTo(m,CV_64FC1);

// Eigen vectors
VectorXd I(n_pixels);
VectorXd u(n_pixels);

// convert image from openCV to Eigen
cv2eigen(m,I);

//
SparseMatrix<double> A(n_pixels,n_pixels);

// Fill sparse martix using triplets
typedef Eigen::Triplet<double> T;
std::vector<T> tripletList;

// Filter parameter (smoothing factor)
//double alpha=-0.1;

// Set values
for(int i=0;i<n_pixels;i++)
{
tripletList.push_back(T(i,i,1+4*alpha));
if((i+1) < n_pixels){tripletList.push_back(T(i,i+1,-alpha));} // +1
if((i-1) >= 0){tripletList.push_back(T(i,i-1,-alpha));} // -1
if((i+src.cols) < n_pixels){tripletList.push_back(T(i,i+src.cols,-alpha));} // +3
if((i-src.cols) >= 0){tripletList.push_back(T(i,i-src.cols,-alpha));} // -3
}

// Boundary values of main diag
tripletList.push_back(T(0,0,1+2*alpha));
for(int i=1;i<src.cols;i++)
{
tripletList.push_back(T(i,i,1+3*alpha));
}

//
tripletList.push_back(T(n_pixels-1,n_pixels-1,1+2*alpha));
for(int i=1;i<src.cols;i++)
{
tripletList.push_back(T(i,n_pixels-i-1,1+3*alpha));
}

// Init sparse matrix
A.setFromTriplets(tripletList.begin(),tripletList.end());

tripletList.clear();
// Solver init
ConjugateGradient<SparseMatrix<double> > cg;
cg.compute(A);
// Solve linear systyem
u = cg.solve(I);
std::cout << "#iterations:     " << cg.iterations() << std::endl;
std::cout << "estimated error: " << cg.error()      << std::endl;
// Get the solution
dst=Mat(n_pixels,1,CV_64FC1);
eigen2cv(u,dst);
dst=dst.reshape(1,src.rows);
dst.convertTo(dst,CV_8UC1);
}


int main(int argc, char* argv[])
{
namedWindow("image");
namedWindow("result");
Mat img=imread("d:\\ImagesForTest\\lena.jpg",1);
imshow("image",img);
waitKey(10);
Mat res;
vector<Mat> ch;
cv::split(img,ch);

for(int i=0;i<3;i++)
{
EnergyFilter(ch[i],res,3);
res.copyTo(ch[i]);
}

cv::merge(ch,res);
// show the resilt
imshow("result",res);
waitKey(0);
return 0;
}

*/