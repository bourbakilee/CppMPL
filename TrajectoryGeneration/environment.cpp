#include "environment.h"
#include <functional>

namespace environment{
	Road::Road(const ArrayXXd& center_line, double ref_grid_length, double ref_grid_width, double lane_width, unsigned int lane_num)
	{
		// center_line is array: (s,x,y,\theta,k), s=linspae(0.,road_length,N)
		this->center_line = center_line;
		this->lane_num = lane_num;
		this->current_lane = (lane_num - 1) / 2;
		this->target_lane = this->current_lane;
		this->length = center_line(center_line.rows()-1, 0);
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
		if (0 <= i < this->lane_num)
		{
			this->current_lane = i;
		}
	}

	//
	void Road::set_target_lane(double q1[])
	{
		double sl[] = { -1.,-1. };
		this->xy2sl(sl, q1[0], q1[1]);
		this->target_lane = std::floor((this->lane_num / 2. + sl[1] / this->lane_width));
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
		if (0 <= i < this->lane_num)
		{
			this->target_lane = i;
		}
	}

	// q - (x,y,\theta,k)
	void Road::sl2xy(double q[], double s, double l)
	{
		if (0. <= s <= this->length)
		{
			// linear interpolate
			// delta_s = center_line(1,0)
			int i = std::floor(s / this->center_line(1, 0));
			double delta1 = s - center_line(i, 0);
			double delta2 = center_line(i + 1, 0) - s;
			ArrayXXd station = (this->center_line.block(i, 1, 1, 4)*delta2 + this->center_line.block(i + 1, 1, 1, 4)*delta1) / this->center_line(1, 0);
			q[0] = station(0, 0) - l*sin(station(0, 2));
			q[1] = station(0,1) + l*cos(station(0,2));
			q[2] = station(0,2);
			q[3] = 1 / (1 / station(0,3) - l);
		}
		else 
		{
			q[0] = q[1] = -1.;
		}
	}

	void Road::ij2xy(double q[], double i, double j)
	{
		if (0 <= i <= this->longitudinal_grid_num)
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
			int i = std::floor(s / this->center_line(1, 0));
			double delta1 = s - center_line(i, 0);
			double delta2 = center_line(i + 1, 0) - s;
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
		int i = std::floor(sl[0] / this->center_line(1, 0));
		double delta1 = sl[0] - center_line(i, 0);
		double delta2 = center_line(i + 1, 0) - sl[0];
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
			sl[1] = (y - station(0,1)) / sin_t;
		}
		else
		{
			sl[1] = (station(0,0) - x) / cos(station(0,0));
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
	void Road::traj2sl(ArrayXXd & traj, ArrayXXd & sl)
	{
		int N = traj.rows();
		sl.resize(N, 2);
		double tmp[] = { -1.,-1. };

#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			this->xy2sl(tmp, traj(i, 2), traj(i, 3));
			sl.row(i) << tmp[0], tmp[1];
		}
	}

	CostMap::CostMap(std::vector<ArrayXXd>& maps, double start_time, double end_time, double resolution): data(maps),start_time(start_time),end_time(end_time),resolution(resolution)
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
	}

	double CostMap::query(Vehicle & vehicle)
	{

		return 0.0;
	}

	Vehicle::Vehicle(double wb, double fo, double ro, double wt) : wheel_base(wb), front_offset(fo), rear_offset(ro), width(wt), length(wb + fo + ro)
	{
	}

	void Vehicle::cover_centers(ArrayXXd & cover_points, ArrayXXd & rear_center_traj)
	{
	}

}