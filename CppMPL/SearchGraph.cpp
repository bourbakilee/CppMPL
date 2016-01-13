#include "SearchGraph.h"
#include <cmath>
#include <exception>

namespace SearchGraph
{

	State::State(double r_s, double r_l, Road * road, double velocity, double acc, double cost, double dk)
	{
		this->time = 0.;
		this->t_i = 0;
		this->length = 0.;
		
		this->r_s = r_s;
		this->r_l = r_l;
		this->r_i = std::llround(r_s / road->grid_length);
		this->r_j = std::llround(r_l / road->grid_width);

		double q[4];
		road->sl2xy(q, r_s, r_l);
		this->x = q[0];
		this->y = q[1];
		this->theta = q[2];
		this->k = q[3];

		this->v = velocity;
		this->v_i = std::llround(velocity / 2);
		this->a = acc;
		this->dk = dk;

		this->cost = cost;
		this->heuristic = 0.;
		this->priority = this->cost + this->heuristic;

		this->parent = nullptr;
		this->reach = false;
		this->extend = false;
	}

	State::State(double x, double y, double theta, double k, Road * road, double velocity, double acc, double cost, double dk)
	{
		this->time = 0.;
		this->t_i = 0;
		this->length = 0.;

		this->x = x;
		this->y = y;
		this->theta = theta;
		this->k = k;

		this->v = velocity;
		this->v_i = std::llround(velocity / 2);
		this->a = acc;
		this->dk = dk;

		double sl[2];
		road->xy2sl(sl, x, y);
		this->r_s = sl[0];
		this->r_l = sl[1];
		this->r_i = std::llround(sl[0] / road->grid_length);
		this->r_j = std::llround(sl[1] / road->grid_width);

		this->cost = cost;
		this->heuristic = 0.;
		this->priority = this->cost + this->heuristic;

		this->parent = nullptr;
		this->reach = false;
		this->extend = false;

	}




	void State::successors(VecSuccs & outs, State_Dict& state_dict, Road * road, StatePtr goal, const std::vector<double>& as, const std::vector<double>& vs, const std::vector<double>& ts, const double* p_lims)
	{
		outs.clear();
		if ((this->v + goal->v)*2.5 > goal->r_s - this->r_s)
		{
			outs.push_back(goal);
		}
#pragma omp parallel for
		for (double n1 : as)
		{
#pragma omp parallel for
			for (double n2 : vs)
			{
#pragma omp parallel for
				for (double n3 : ts)
				{
					double v = std::min(std::max(this->v + n1*n3, p_lims[3]), p_lims[2]);
					double l = this->r_l + n2*n3;
					double s = this->r_s + (this->v + v) / 2 * n3;
					int r_i = std::llround(s / road->grid_length);
					int r_j = std::llround(l / road->grid_width);
					int v_i = std::llround(v / 2.);
					
					State_Index si = std::make_tuple(r_i, r_j, v_i);
					StatePtr state = nullptr;
					try
					{
						state = state_dict.at(si);
					}
					catch (const std::out_of_range& e)
					{
						if (this->r_i <= r_i <= goal->r_i && std::abs(r_j) * 2 < road->lateral_grid_num)
						{
							state = std::make_shared<State>(State(r_i, r_j, road, v_i*2., n1));
						}
					}
					if (state != nullptr)
					{
						outs.push_back(state);
					}
				}
			}
		}
	}

	// the current state is successfully connected by traj from pareant state
	bool State::update(StatePtr current, StatePtr parent, double cost, ArrayXXd & traj, Traj_Dict & traj_dict, const HeuristicMap& hm, StatePtr goal, const Vehicle& veh, double delta_t)
	{
		bool res = false;

		if (!parent->extend)
		{
			parent->extend = true;
		}
		double c = cost + parent->cost;
		if (current->cost > c)
		{
			int N = traj.rows() - 1;
			current->time = traj(N, 0);
			current->t_i = std::llround(traj(N, 0) / delta_t);
			current->length = traj(N, 1);
			if (!current->reach)
			{
				current->reach = true;
		    }
			if (current->parent != nullptr)
			{
				traj_dict.erase(Traj_Index(current->parent, current));
			}
			current->parent = parent;
			traj_dict[Traj_Index(current->parent, current)] = ArrayXXd(traj);
			current->cost = c;
			current->heuristic = hm.query(current, veh, goal);
			current->priority = current->cost + current->heuristic;
			res = true;
		}

		return res;
	}

	void State::post_process(StatePtr current, StatePtr successor, Eval_Res & res, ArrayXXd & traj, 
		PQ & pq, State_Dict & state_dict, Traj_Dict & traj_dict, StatePtr goal, 
		const Vehicle & veh, Road * road, const CostMap & cost_map, const HeuristicMap & hm, 
		sqlite3 * db, const double * weights)
	{
		if (!res.second)  // the trajectory is not truncated
		{
			if (State::update(successor, current, res.first, traj, traj_dict, hm, goal, veh))
			{
				pq.push(successor);
				if (successor != goal)
				{
					State_Index si = std::make_tuple(successor->r_i, successor->r_j, successor->v_i);
					auto search = state_dict.find(si);
					if (search == state_dict.end())
					{
						state_dict[si] = successor;
					}
				}
			}
		}
		else  // the trajectory is truncated
		{
			successor = std::make_shared<State>(traj, road);
			State_Index si = std::make_tuple(successor->r_i, successor->r_j, successor->v_i);
			auto search = state_dict.find(si);
			if (search == state_dict.end())
			{
				state_dict[si] = successor;
			}
			else
			{
				StatePtr state = state_dict[si];
				if (State::distance(successor, state) < 1.e-4)  // successor is close to state
				{
					successor = state;
				}
				else  // bias successor to state
				{
					if (connect(current, state, traj, db))  // successfully compute the trajectory from current to state
					{
						res = eval_traj(traj, veh, cost_map, weights, kinematic_limits, road, false);
						if (!std::isinf(res.first)) // successfully connect current to state
						{
							successor = state;
						}
						else
						{
							successor = nullptr;
						}
					}
					else
					{
						successor = nullptr;
					}
				}
			}
			if (successor != nullptr)
			{
				if (State::update(successor, current, res.first, traj, traj_dict, hm, goal, veh))
				{
					pq.push(successor);
				}
			}
		}
	}


	double HeuristicMap::query(StatePtr current, const Vehicle & veh, StatePtr goal) const
	{
		double heuristic = 0.;
		ArrayXXd cover_points = ArrayXXd::Zero(1, 6);
		ArrayXXd rear_center(1, 9);
		rear_center << current->time, current->length, current->x, current->y, current->theta, current->k, current->dk, current->v, current->a;
		veh.cover_centers(cover_points, rear_center);
		Array<unsigned int, 1, 6> cover_index = (cover_points/this->resolution).cast<unsigned int>();

		if (this->isdynamic)
		{
			if (this->start_time <= current->time <= this->end_time)
			{
				int t_index = (int)std::floor(current->time / this->delta_t);
				heuristic = (((t_index + 1)*this->delta_t - current->time)
					*(1.5*this->data[t_index](cover_index(0, 1), cover_index(0, 0))
						+ this->data[t_index](cover_index(0, 3), cover_index(0, 2))
						+ 0.5*this->data[t_index](cover_index(0, 5), cover_index(0, 4)))
					+ (current->time - t_index*this->delta_t))
					*(1.5*this->data[t_index + 1](cover_index(0, 1), cover_index(0, 0))
						+ this->data[t_index + 1](cover_index(0, 3), cover_index(0, 2))
						+ 0.5*this->data[t_index + 1](cover_index(0, 5), cover_index(0, 4)))
					/ this->delta_t;
			}
		}
		else
		{
			heuristic = 1.5*this->data[0](cover_index(0, 1), cover_index(0, 0))
				+ this->data[0](cover_index(0, 3), cover_index(0, 2))
				+ 0.5*this->data[0](cover_index(0, 5), cover_index(0, 4));
		}
		return heuristic*5.;
	}


#ifdef WITH_SQLITE3
	bool connect(StatePtr s1, StatePtr s2, ArrayXXd & traj, sqlite3 * db)
	{
		return false;
	}
#endif

}
