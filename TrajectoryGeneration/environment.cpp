#include "environment.h"
#include <functional>
#include <vector>

#include <iostream>

namespace environment{

	Vehicle::Vehicle(const Vehicle& veh)
	{
		this->cover_distance = veh.cover_distance;
		this->cover_radius = veh.cover_radius;
		this->front_offset = veh.front_offset;
		this->geometric_center = veh.geometric_center;
		this->length = veh.length;
		this->rear_offset = veh.rear_offset;
		this->wheel_base = veh.wheel_base;
		this->width = veh.wheel_base;
	}


	Vehicle::Vehicle(double wb, double fo, double ro, double wt) : wheel_base(wb), front_offset(fo), rear_offset(ro), width(wt), length(wb + fo + ro)
	{
		this->cover_radius = std::sqrt(this->length*this->length / 9. + this->width*this->width)/2.;
		this->cover_distance = 2.*this->length / 3.;
		this->geometric_center = (wb + fo - ro) / 2; // distances to rear center
	}


	// when rear_center_traj has only one row, it represent the state of vehicle
	// cover_points - [(x1,y1,x2,y2,x3,y3)] - front, center, rear
	// rear_center_traj - [(t,s,x,y,theta,....)]
	void Vehicle::cover_centers(ArrayXXd & cover_points, const ArrayXXd & rear_center_traj) const
	{
		int N = rear_center_traj.rows();
		if (cover_points.rows() != N)
		{
			cover_points.resize(N, 6);
		}
		ArrayXXd cos_t = rear_center_traj.col(4).cos();
		ArrayXXd sin_t = rear_center_traj.col(4).sin();
		cover_points.col(2) = rear_center_traj.col(2) + this->geometric_center*cos_t; // x2
		cover_points.col(3) = rear_center_traj.col(3) + this->geometric_center*sin_t; // y2
		cover_points.col(0) = cover_points.col(2) + this->cover_distance*cos_t; // x1
		cover_points.col(1) = cover_points.col(3) + this->cover_distance*sin_t; // y1
		cover_points.col(4) = cover_points.col(2) - this->cover_distance*cos_t; // x3
		cover_points.col(5) = cover_points.col(3) - this->cover_distance*sin_t; // y3
	}


	Road::Road(const Road& road)
	{
		this->center_line = road.center_line;
		this->current_lane = road.current_lane;
		this->current_lane_center_line_offset = road.current_lane_center_line_offset;
		this->grid_length = road.grid_length;
		this->grid_width = road.grid_width;
		this->lane_grid_num = road.lane_grid_num;
		this->lane_num = road.lane_num;
		this->lane_width = road.lane_width;
		this->lateral_grid_num = road.lateral_grid_num;
		this->length = road.length;
		this->longitudinal_grid_num = road.longitudinal_grid_num;
		this->target_lane = road.target_lane;
		this->target_lane_center_line_offset = road.target_lane_center_line_offset;
		this->width = road.width;
	}

	Road::Road(const ArrayXXd& center_line, double ref_grid_length, double ref_grid_width, double lane_width, unsigned int lane_num)
	{
		// center_line is array: (s,x,y,\theta,k), s=linspae(0.,road_length,N)
		this->center_line = center_line;
		this->lane_num = lane_num;
		this->current_lane = (lane_num - 1) / 2;
		this->target_lane = this->current_lane;
		this->current_lane_center_line_offset = (this->current_lane - this->lane_num / 2. + 0.5)*this->lane_width;
		this->target_lane_center_line_offset = (this->target_lane - this->lane_num / 2. + 0.5)*this->lane_width;
		this->length = center_line(center_line.rows() -1, 0);
		this->lane_width = lane_width;
		this->width = lane_width*lane_num;

		this->longitudinal_grid_num = (int)std::ceil(this->length / ref_grid_length);
		this->lane_grid_num = 2*std::llround(0.5*this->lane_width / ref_grid_width);
		this->lateral_grid_num = this->lane_grid_num*this->lane_num;
		this->grid_length = this->length / this->longitudinal_grid_num;
		this->grid_width = this->lane_width / this->lane_grid_num;
	}

	void Road::set_center_line(const ArrayXXd& line, double ref_grid_length)
	{
		this->center_line = center_line;
		this->length = center_line(center_line.rows() - 1, 0);
		this->longitudinal_grid_num = (int)std::ceil(this->length / ref_grid_length);
		this->grid_length = this->length / this->longitudinal_grid_num;
	}

	//
	void Road::set_current_lane(double q0[])
	{
		double sl[] = { -1.,-1. };
		this->xy2sl(sl, q0[0], q0[1]);
		this->current_lane = std::floor((this->lane_num / 2. + sl[1] / this->lane_width));
		this->current_lane_center_line_offset = (this->current_lane - this->lane_num + 0.5)*this->lane_width;
		/*
		if (this->lane_num % 2)
		{
			this->current_lane = (this->lane_num - 1) / 2 + std::floor(sl[1] / this->lane_width + 0.5);
		}
		else
		{
			this->current_lane = this->lane_num / 2 + std::floor(sl[1] / this->lane_width);
		}
		*/
	}

	void Road::set_current_lane(unsigned int i)
	{
		if (0 <= i && i < this->lane_num)
		{
			this->current_lane = i;
			this->current_lane_center_line_offset = (this->current_lane - this->lane_num + 0.5)*this->lane_width;
		}
	}

	//
	void Road::set_target_lane(double q1[])
	{
		double sl[] = { -1.,-1. };
		this->xy2sl(sl, q1[0], q1[1]);
		this->target_lane = std::floor((this->lane_num / 2. + sl[1] / this->lane_width));
		this->target_lane_center_line_offset = (this->target_lane - this->lane_num + 0.5)*this->lane_width;
		/*
		if (this->lane_num % 2)
		{
			this->target_lane = (this->lane_num - 1) / 2 + std::floor(sl[1] / this->lane_width + 0.5);
		}
		else
		{
			this->target_lane = this->lane_num / 2 + std::floor(sl[1] / this->lane_width);
		}
		*/
	}

	void Road::set_target_lane(unsigned int i)
	{
		if (0 <= i && i < this->lane_num)
		{
			this->target_lane = i;
			this->target_lane_center_line_offset = (this->target_lane - this->lane_num + 0.5)*this->lane_width;
		}
	}

	// q - (x,y,\theta,k)
	void Road::sl2xy(double q[], double s, double l)
	{
		if (0. <= s && s <= this->length)
		{
			// linear interpolate
			// delta_s = center_line(1,0)
			int i = std::min(std::max(std::floor(s / this->center_line(1, 0)), 0.), this->longitudinal_grid_num - 2.);
			double delta1 = s - center_line(i, 0);
			double delta2 = center_line(i + 1, 0) - s;
			ArrayXXd station = (this->center_line.block(i, 1, 1, 4)*delta2 + this->center_line.block(i + 1, 1, 1, 4)*delta1) / this->center_line(1, 0);
			q[0] = station(0, 0) - l*sin(station(0, 2));
			q[1] = station(0,1) + l*cos(station(0,2));
			q[2] = station(0,2);
			if (std::abs(station(0, 3)) < 1.e-8)
				q[3] = 0.;
			else
			    q[3] = 1. / (1. / station(0,3) - l);
		}
		else 
		{
			q[0] = q[1] = -1.;
		}
	}

	void Road::ij2xy(double q[], double i, double j)
	{
		if (0 <= i && i <= this->longitudinal_grid_num)
		{
			this->sl2xy(q, i*this->grid_length, j*this->grid_width);
		}
		else
		{
			q[0] = q[1] = -1.;
		}
	}

	void Road::xy2sl(double sl[], double x, double y)
	{
		std::function<double(double)> f = [this, x, y](double s) {
			// linear interpolate
			// delta_s = center_line(1,0)
			int i = std::min(std::max(std::floor(s / this->center_line(1, 0)), 0.), this->longitudinal_grid_num-2.);
			double delta1 = s - this->center_line(i, 0);
			double delta2 = this->center_line(i + 1, 0) - s;
			ArrayXXd station = (this->center_line.block(i, 1, 1, 3)*delta2 + this->center_line.block(i + 1, 1, 1, 3)*delta1) / this->center_line(1, 0);
			return (x - station(0,0))*cos(station(0,2)) + (y - station(0,1))*sin(station(0,2));
		};
		double a = 0., b = this->length, fs;

		do {
			sl[0] = (a + b) / 2;
			fs = f(sl[0]);
			if (f(a)*fs < 0)
			{
				b = sl[0];
			}
			else if (f(b)*fs < 0)
			{
				a = sl[0];
			}
			else
			{
				sl[0] = sl[1] = -1;
				return;
			}
		} while (std::abs(fs) > 1.e-3);
		// linear interpolate
		// delta_s = center_line(1,0)
		int i = std::min(std::max(std::floor(sl[0] / this->center_line(1, 0)), 0.), this->longitudinal_grid_num - 2.);
		double delta1 = sl[0] - this->center_line(i, 0);
		double delta2 = this->center_line(i + 1, 0) - sl[0];
		ArrayXXd station = (this->center_line.block(i, 1, 1, 3)*delta2 + this->center_line.block(i + 1, 1, 1, 3)*delta1) / this->center_line(1, 0);
		/*
		if abs(tmp[2]) < 1.e-4:
            l0 = (y-tmp[1])/np.cos(tmp[2])
        else:
            l0 = (tmp[0]-x)/np.sin(tmp[2])
		*/
		double sin_t = std::sin(station(0,2));
		if (std::abs(sin_t) < 1.e-4)
		{
			sl[1] = (y - station(0,1)) / cos(station(0, 2));
		}
		else
		{
			sl[1] = (station(0,0) - x) / sin_t;
		}
	}

	void Road::xy2ij(unsigned int ij[], double x, double y)
	{
		ij[0] = ij[1] = -1;
		double sl[] = { -1.,-1. };
		this->xy2sl(sl, x, y);
		if (sl[0] >= 0)
		{
			ij[0] = std::llround(sl[0] / this->grid_length);
			ij[1] = std::llround(sl[1] / this->grid_width);
		}
	}

	// traj - array of points on trajectory - [(t,s,x,y,theta,k,dk,v,a)]
	// sl - array of coordinates in SL frame
	void Road::traj2sl(const ArrayXXd & traj, ArrayXXd & sl)
	{
		int N = traj.rows();
		// std::cout << "traj2sl - Rows of traj: " << N << std::endl;
		if (sl.rows() != N || sl.cols()<2)
		{
			sl.resize(N, 2);
		}
		// std::cout << "xy2sl - Size of SL array: " << sl.rows() << " " << sl.cols() << std::endl;
		double tmp[] = { -1.,-1. };

#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{

			//  std::cout << "xy2sl - i: " << i << " of " << N << std::endl;
			//  std::cout << "xy2sl - xy: " << traj(i, 2) << " " << traj(i, 3) << std::endl;

			this->xy2sl(tmp, traj(i, 2), traj(i, 3)); 
			sl.row(i) << tmp[0], tmp[1];
			
			// std::cout << "xy2sl - SL: " << sl.row(i) << std::endl;

		}
		// std::cout << "traj2sl - sl coordinates of traj: \n";
		// std::cout << sl << std::endl;
	}


	CostMap::CostMap(const CostMap& cm)
	{
		this->cols = cm.cols;
		this->data = cm.data;
		this->delta_t = cm.delta_t;
		this->end_time = cm.end_time;
		this->isdynamic = cm.isdynamic;
		this->num = cm.num;
		this->resolution = cm.resolution;
		this->rows = cm.rows;
		this->start_time = cm.start_time;
	}



	CostMap::CostMap(const std::vector<ArrayXXd>& maps, double start_time, double end_time, double resolution): data(maps),start_time(start_time),end_time(end_time),resolution(resolution)
	{
		if (this->data.size() == 1)
		{
			this->start_time = this->end_time = this->delta_t = 0.;
			this->isdynamic = false;
		}
		else
		{
			this->delta_t = (this->end_time - this->start_time) / (this->data.size()-1);
			this->isdynamic = true;
		}
		this->cols = maps[0].cols();
		this->rows = maps[0].rows();
		this->num = maps.size();
	}

	/*
	double CostMap::query(Vehicle & vehicle)
	{

		return 0.0;
	}
	*/

	// cover_points - [(x1,y1,x2,y2,x3,y3)]
	// traj - [(t,s,x,y,theta....)]
	// (x,y) <-> (col, row)
	void CostMap::query(const Vehicle& vehicle, const ArrayXXd & traj, ArrayXXd& cost) const
	{
		int N = traj.rows();
		if (cost.rows() != N && cost.cols() != 1)
		{
			cost.resize(N, 1);
		}

		ArrayXXd cover_points(N, 6);
		vehicle.cover_centers(cover_points, traj);
		// covert the global coordinates to index of costmap
		Array<unsigned int, Dynamic, 6> cover_indexes = (cover_points / this->resolution).cast<unsigned int>();
		if (this->isdynamic)
		{
#pragma omp parallel for
			for (int i = 0; i < N; i++)
			{
				if (this->start_time <= traj(i, 0) && traj(i, 0) <= this->end_time)
				{
					int t_index = (int)std::floor(traj(i, 0) / this->delta_t);
					cost(i, 0) = (((t_index + 1)*this->delta_t - traj(i, 0))
						*(1.5*this->data[t_index](cover_indexes(i, 1), cover_indexes(i, 0)) 
							+ this->data[t_index](cover_indexes(i, 3), cover_indexes(i, 2)) 
							+ 0.5*this->data[t_index](cover_indexes(i, 5), cover_indexes(i, 4))) 
						+ (traj(i, 0) - t_index*this->delta_t))
						*(1.5*this->data[t_index + 1](cover_indexes(i, 1), cover_indexes(i, 0)) 
							+ this->data[t_index + 1](cover_indexes(i, 3), cover_indexes(i, 2)) 
							+ 0.5*this->data[t_index + 1](cover_indexes(i, 5), cover_indexes(i, 4))) 
						/ this->delta_t;
				}
			}
		}
		else
		{
#pragma omp parallel for
			for (int i = 0; i < N; i++)
			{
					// int t_index = (int)std::floor(traj(i, 0) / this->delta_t);
					cost(i, 0) = 1.5*this->data[0](cover_indexes(i, 1), cover_indexes(i, 0))
							+ this->data[0](cover_indexes(i, 3), cover_indexes(i, 2))
							+ 0.5*this->data[0](cover_indexes(i, 5), cover_indexes(i, 4));
			}
		}
	}


	void read_map(const char* filename, ArrayXXd& map, const int rows, const int cols)
	{
		if (map.cols() != cols || map.rows() != rows)
		{
			map.resize(rows, cols);
		}
		std::ifstream input(filename);
		for (int i = 0; i < rows * cols; i++)
		{
			// if input number is negative, it should be modified to inf
			input >> map(i / cols, i % rows);
		}

		ArrayXXd flag = (map < 0.).cast<double>();
#pragma omp parallel for
		for (int i = 0; i < rows; i++)
		{
#pragma omp parallel for
			for (int j = 0; j < cols; j++)
			{
				if (flag(i, j) > 0.5)
				{
					map(i, j) =  std::numeric_limits<double>::infinity();
				}
			}
		}
	}


	void read_road_center_line(const char* filename, ArrayXXd& center_line)
	{
		std::vector<double> vec_tmp;
		std::ifstream input(filename);
		double tmp;
		while (!(input >> tmp).fail()) 
		{
			vec_tmp.push_back(tmp);
		}
		center_line.resize(vec_tmp.size() / 5, 5);
#pragma omp parallel for
		for (int i = 0; i < center_line.rows(); i++)
		{
#pragma omp parallel for
			for (int j = 0; j < center_line.cols(); j++)
			{
				center_line(i, j) = vec_tmp[i*center_line.cols() + j] < 0 ? std::numeric_limits<double>::infinity() : vec_tmp[i*center_line.cols() + j];
			}
		}
	}

// end of namespace environment
}