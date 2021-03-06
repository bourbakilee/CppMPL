#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>
#include <cmath>
#include <fstream>
// define EIGEN_USE_MKL_ALL
#include <Eigen\Dense>


namespace environment 
{
	using namespace Eigen;

	struct Vehicle {
		// members
		double wheel_base;
		double front_offset;
		double rear_offset;
		double width;
		double length;
		double cover_radius;
		double cover_distance;
		// distances to rear center
		double geometric_center;
		// double front_cover_center;
		// double rear_cover_center;

		// construct functions
		Vehicle():Vehicle(2.94, 0.988,1.28,2.029) {}
		Vehicle(const Vehicle& veh);
		Vehicle(double wb, double fo, double ro, double wt);

		// member functions
		// when rear_center_traj has only one row, it represent the state of vehicle
		// cover_points - [(x1,y1,x2,y2,x3,y3)]
		// rear_center_traj - [(t,s,x,y,theta....)]
		void cover_centers(ArrayXXd& cover_points, const ArrayXXd& rear_center_traj) const;
		// return coordinates of vertexes of vehicle body
		// void vertexes(ArrayXXd& vertex_list, ArrayXXd& rear_center_traj)
	};


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
		double current_lane_center_line_offset;
		double target_lane_center_line_offset;
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
		Road():Road(ArrayXXd::Ones(1,5),1.5,0.5,3.5,3) {}
		Road(const Road& road);
		// the first column of road center_line array must be equidistantly
		Road(const ArrayXXd& center_line, double ref_grid_length = 1., double ref_grid_width = 0.5, double lane_width = 3.5, unsigned int lane_num = 3);

		// member functions
		void set_center_line(const ArrayXXd& line, double ref_grid_length=1.5);
		void set_current_lane(double q0[]);
		void set_current_lane(unsigned int i);
		void set_target_lane(double q1[]);
		void set_target_lane(unsigned int i);
		void sl2xy(double q[], double s, double l);
		void ij2xy(double q[], double i, double j);
		void xy2sl(double sl[], double x, double y);
		void xy2ij(unsigned int ij[], double x, double y);
		// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
		void traj2sl(const ArrayXXd& traj, ArrayXXd& sl);
		// extend the road length
		// void extend_road(const ArrayXXd& line);
		// coid extend_road(const Road& road);
	};

	struct CostMap {
		bool isdynamic;
		double start_time;
		double end_time;
		double resolution;
		unsigned int cols;
		unsigned int rows;
		unsigned int num;
		std::vector<ArrayXXd> data; // vector<ArrayXXd> data
		double delta_t;

		// construct functions
		CostMap() :isdynamic(false),start_time(0.), end_time(0.),resolution(0.2), cols(500), rows(500), data{ ArrayXXd::Ones(500, 500) }, delta_t(0.) {}
		CostMap(const CostMap& cm);
		CostMap(ArrayXXd&map, double resolution = 0.2) :isdynamic(false), start_time(0.), end_time(0.), resolution(resolution), cols(map.cols()), rows(map.rows()), data{ map }, delta_t(0.) {}
		CostMap(const std::vector<ArrayXXd>& maps, double start_time = 0., double end_time = 0., double resolution = 0.2);

		// member functions
		// query the cost of a vehicle
		// double query(Vehicle& vehicle);
		void query(const Vehicle& vehicle, const ArrayXXd& traj, ArrayXXd& cost) const;
	};

	void read_map(const char* filename, ArrayXXd& map, int rows = 500, int cols = 500);

	void read_road_center_line(const char* filename, ArrayXXd& center_line);

}


#endif
