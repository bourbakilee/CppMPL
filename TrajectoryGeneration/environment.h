#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <Eigen\Dense>
using namespace Eigen;

namespace environment 
{
	// no need to use class
	// easy to be converted to class
	struct Road
	{
		// members
		// array of points on center line: (s, x, y, theta, k)
		ArrayXXd center_line;
		// ID: right -> left, 0 -> lane_num - 1
		unsigned int lane_num; // default: 3
		unsigned int current_lane;
		unsigned int target_lane;
		double lane_width;
		double length;
		double width; 
		//
		unsigned int longitudinal_grid_num;
		unsigned int lane_grid_num;
		unsigned int lateral_grid_num;
		double grid_length;
		double grid_width;
		
		// construct functions
		Road() {}
		Road(ArrayXXd& center_line, double ref_grid_length = 1.5, double ref_grid_width = 0.5, double lane_width = 3.5, unsigned int lane_num = 3);

		// member functions
		void set_current_lane(double q0[]);
		void set_target_lane(double q1[]);
		void sl2xy(double q[], double s, double l);
		void ij2xy(double q[], double s, double l);
		void xy2sl(double sl[], double x, double y);
		void xy2ij(double ij[], double x, double y);
		void traj2sl(ArrayXXd& traj, ArrayXXd& sl);
	};

	struct CostMap {
		bool isdynamic;
		double resolution;
		unsigned int cols;
		unsigned int rows;
		std::vector<ArrayXXd> data; // vector<ArrayXXd> data

		CostMap() :isdynamic(false),resolution(0.2), cols(500), rows(500), data{ ArrayXXd::Ones(500, 500) } {}
		CostMap(ArrayXXd&map, double resolution = 0.2) :isdynamic(false),resolution(0.2), cols(map.cols()), rows(map.rows()), data{ map } {}
	};



}


#endif
