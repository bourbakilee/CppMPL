#include "SearchGraph.h"

namespace SearchGraph
{
	// default construct function
	State::State()
	{
		this->time = 0.;
		this->length = 0.;
		this->x = 0.;
		this->y = 0.;
		this->theta = 0.;
		this->k = 0.;
		this->dk = 0.;
		this->v = 0.;
		this->a = 0.;
		this->r_s = 0.;
		this->r_l = 0.;
		this->r_i = 0;
		this->r_j = 0;
		this->cost = inf;
		this->priority = inf;
	}

	State::State(int i, int j, double velocity, double acc, Road * road, State * goal)
	{
		this->time = 0.;
		this->length = 0.;
		this->dk = 0.;
		this->v = velocity;
		this->a = acc;
		this->cost = 0.;
		//
		double q[4];
		road->ij2xy(q, i, j);
		this->x = q[0];
		this->y = q[1];
		this->theta = q[2];
		this->k = q[3];
		//
		this->r_i = i;
		this->r_j = j;
		this->r_s = i*road->grid_length;
		this->r_l = j*road->grid_width;
		//
		this->priority = this->cost + this->heuristic(goal);
	}

	State::State(double r_s, double r_l, double velocity, double acc, Road * road, State * goal)
	{
		this->time = 0.;
		this->length = 0.;
		this->dk = 0.;
		this->v = velocity;
		this->a = acc;
		this->cost = 0.;
		//
		double q[4];
		road->sl2xy(q, r_s, r_l);
		this->x = q[0];
		this->y = q[1];
		this->theta = q[2];
		this->k = q[3];
		//
		this->r_i = std::llround(r_s / road->grid_length);
		this->r_j = std::llround(r_l / road->grid_width);
		this->r_s = r_s;
		this->r_l = r_l;
		//
		this->priority = this->cost + this->heuristic(goal);
	}

	State::State(double x, double y, double theta, double k, double v, double acc, State * goal, Road * road)
	{
		this->time = 0.;
		this->length = 0.;
		this->x = x;
		this->y = y;
		this->theta = theta;
		this->k = k;
		this->dk = 0.;
		this->v = v;
		this->a = acc;
		this->cost = 0.;
		//
		if (road != nullptr)
		{
			double sl[2];
			road->xy2sl(sl, x, y);
			this->r_s = sl[0];
			this->r_l = sl[1];
			this->r_i = std::llround(sl[0] / road->grid_length);
			this->r_j = std::llround(sl[1] / road->grid_width);
		}
		else
		{
			this->r_s = 0.;
			this->r_l = 0.;
			this->r_i = 0;
			this->r_j = 0;
		}
		this->priority = this->cost + this->heuristic(goal);
	}


	// prev and this will be all modified: dk, ag ....
	// traj has been evaluated, and has more than 1 row
	// traj - array of points on trajectory - [(t, s, x, y, theta, k, dk, v, a)]
	// ref_rows - rows of traj before eval
	State::State(State*prev, ArrayXXd & traj, int ref_rows, double traj_cost, Road * road, State * goal)
	{
		int N = traj.rows();
		// prev
		prev->dk = traj(0, 6);
		// this
		this->cost = prev->cost + traj_cost;
		this->time = traj(N, 0);
		this->length = traj(N, 1);
		this->x = traj(N, 2);
		this->y = traj(N, 3);
		this->theta = traj(N, 4);
		this->k = traj(N, 5);
		this->dk = traj(N, 6);
		this->v = traj(N, 7);
		this->a = traj(N, 8);
		if (road != nullptr)
		{
			//this->r_s
		}

	}

	double State::heuristic(State * goal)
	{
		return 0.0;
	}

	void State::extend()
	{
	}

	double heuristic(State * s1, State * s2)
	{
		return s1->heuristic(s2);
	}

	bool connect(State * s1, State * s2, ArrayXXd & traj)
	{
		return false;
	}

}
