#include "SearchGraph.h"
#include <cmath>
#include <exception>
#include <algorithm> 
#include <iostream>
// #include <cassert>

namespace SearchGraph
{

	HeuristicMap::HeuristicMap(const HeuristicMap& hm)
	{
		this->cols = hm.cols;
		this->data = hm.data;
		this->delta_t = hm.delta_t;
		this->end_time = hm.end_time;
		this->isdynamic = hm.isdynamic;
		this->num = hm.num;
		this->resolution = hm.resolution;
		this->rows = hm.rows;
		this->start_time = hm.start_time;
	}



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




	void State::successors(VecSuccs & outs, State_Dict& state_dict, Road * road, StatePtr goal, const std::vector<double>& as, const std::vector<double>& vs, const std::vector<double>& ts, const double* p_lims) const
	{
		outs.clear();
		if ((this->v + goal->v)*2.5 > goal->r_s - this->r_s)
		{
			outs.push_back(goal);
		}
// #pragma omp parallel for
		// for (const double& n1 : as)
		for (int n1 = 0; n1 < as.size();n1++)
		{
// #pragma omp parallel for
			// for (const double& n2 : vs)
			for (int n2 = 0; n2 < vs.size();n2++)
			{
// #pragma omp parallel for
				// for (const double& n3 : ts)
				for (int n3 = 0; n3 < ts.size();n3++)
				{
					double v = std::min(std::max(this->v + as[n1]*ts[n3], p_lims[3]), p_lims[2]);
					double l = this->r_l + vs[n2]*ts[n3];
					double s = this->r_s + (this->v + v) / 2 * ts[n3];
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
						if ((this->r_i < r_i && r_i < goal->r_i) && (std::abs(r_j) * 2 < road->lateral_grid_num))
						{
							state = std::make_shared<State>(State(r_i, r_j, road, v_i*2., as[n1]));
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



	// the current state is successfully connected by traj from pareant state.
	// the time, length, parent, cost, heuristic, reach info of current state must be updated.
	// the extend info of parent must be updated.
	bool State::update(StatePtr successor, StatePtr parent, double cost, const ArrayXXd & traj, Traj_Dict & traj_dict, const HeuristicMap& hm, StatePtr goal, const Vehicle& veh, double delta_t)
	{
		bool res = false;

		if (!parent->extend)
		{
			parent->extend = true;
		}
		double c = cost + parent->cost;
		if (successor->cost > c)
		{
			int N = traj.rows() - 1;
            // assert(N>0);
	        successor->time = traj(N, 0);
			successor->t_i = std::llround(traj(N, 0) / delta_t);
			successor->length = traj(N, 1);
			if (!successor->reach)
			{
				successor->reach = true;
		  }
			if (successor->parent != nullptr)
			{
				traj_dict.erase(std::make_tuple(successor->parent, successor));
			}
			successor->parent = parent;
			traj_dict[std::make_tuple(successor->parent, successor)] = ArrayXXd(traj);
			successor->cost = c;
			successor->heuristic = hm.query(successor, veh, goal);
			successor->priority = successor->cost + successor->heuristic;
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
				if (State::distance(successor, state) < 1.e-3)  // successor is close to state
				{
					successor = state;
				}
				else  // bias successor to state
				{
					if (connect(current, state, traj, db))  // successfully compute the trajectory from current to state
					{
						res = eval_traj(traj, veh, cost_map, road, false, weights, kinematic_limits);
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
			if (this->start_time <= current->time && current->time <= this->end_time)
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
		return heuristic*8.;
	}


#ifdef WITH_SQLITE3
	bool connect(StatePtr s1, StatePtr s2, ArrayXXd & traj, sqlite3 * db)
	{
		bool res = false;
		VectorXd r(5);
		VectorXd q0(4);
		q0 << s1->x, s1->y, s1->theta, s1->k;
		VectorXd q1(4);
		q1 << s2->x, s2->y, s2->theta, s2->k;
		spiral3::spiral3(r, q0, q1, db);
		if (r[4] > 0)
		{
			spiral3::path(traj, r, q0, q1);
			double u[4];
			trajectory::velocity(u, s1->v, s1->a, s2->v, r[4]);
			if (u[3] > 0)
			{
				double r1[] = { r[0],r[1],r[2],r[3],r[4] };
				trajectory::trajectory(traj, r1, u, s1->length, s1->time);
				res = true;
			}
		}

		return res;
	}
#endif


	bool Astar(PQ& pq, StatePtr goal, State_Dict & state_dict, Traj_Dict & traj_dict, Road * road, const Vehicle & veh, CostMap & cost_map, HeuristicMap & hm, sqlite3 * db)
	{
		// bool res = false;

		StatePtr current;
		VecSuccs vec_succ;
		ArrayXXd traj = ArrayXXd::Zero(1, 9);
		//  Eval_Res eval_res;

		int i = 0;

		while (!goal->extend && !pq.empty())
		{
			int j = 0;
			current = pq.top();
			//std::cout << "Size of Priority Queue: " << pq.size() << std::endl;
			i++;
			// std::cout << "Start of " << i << "th point. " << current->time<<"\t" << current->length << "\t" << current->x << "\t" << current->y << "\t" << current->theta << "\t" << current->k << "\t" << current->dk << "\t" << current->v << "\t" << current->a << std::endl;
			// std::cout << "SL: " << current->r_s << "\t" << current->r_l << std::endl;
			// std::cout << "IJ: " << current->r_i << "\t" << current->r_j << std::endl;
			pq.pop();
			current->successors(vec_succ, state_dict, road, goal);
			current->extend = true;
			for (auto successor : vec_succ)
			{
				j++;
				// std::cout << j << "th successor: " << successor->time << "\t" << successor->length << "\t" << successor->x << "\t" << successor->y << "\t" << successor->theta << "\t" << successor->k << "\t" << successor->dk << "\t" << successor->v << "\t" << successor->a << std::endl;
				// std::cout << "SL: " << successor->r_s << "\t" << successor->r_l << std::endl;
				// std::cout << "IJ: " << successor->r_i << "\t" << successor->r_j << std::endl;
				if (connect(current, successor, traj, db))
				{
					Eval_Res eval_res;
					// std::cout << "Rows of traj: " << traj.rows() << std::endl;
					if (successor == goal)
					{
						eval_res = eval_traj(traj, veh, cost_map, road, false);
					}
					else
					{
						eval_res = eval_traj(traj, veh, cost_map, road, true);
					}
					if (!std::isinf(eval_res.first) && traj.rows() > 2)
					{
						State::post_process(current, successor, eval_res, traj, pq, state_dict, traj_dict, goal, veh, road, cost_map, hm, db);
					}
				}
				//  std::cout << "End of " << j << "th trajectory. " << traj.row((traj.rows() - 1)) << std::endl;
			}
			// std::cout <<"End of "<< i << "th point. " << traj.row((traj.rows() - 1)) << std::endl;
		}
		//  std::cout << i << "th point. " << traj.row((traj.rows() - 1)) << std::endl;
		//  std::cout << "Size of State Dict: " << state_dict.size() << std::endl;
		//  std::cout << "Size of Trajectory Dict: " << traj_dict.size() << std::endl;
		//  std::cout << "Cost of Goal State: " << goal->cost << std::endl;
		return goal->extend;
	}
	
	Simulation::Simulation(MODE mode, const char * road, const char * costmap, const char * heuristicmap, const char * dbfile)
	{
		this->mode = mode;
		read_road_center_line(road, this->traj);
		this->road = Road(this->traj);
		std::cout << "Road Length: " << this->road.length << std::endl;
		read_map(costmap, this->traj);
		this->cost_map = CostMap(this->traj);
		std::cout << "Cost Map Size: " << this->cost_map.rows << " X " << this->cost_map.cols << std::endl;
		read_map(heuristicmap, this->traj);
		this->heuristic_map = HeuristicMap(std::vector<ArrayXXd>{this->traj}, this->cost_map);
		std::cout << "Heuristic Map Size: " << this->heuristic_map.rows << " X " << this->heuristic_map.cols << std::endl;
		sqlite3_open(dbfile, &db);
		std::cout << "DataBase" << std::endl;
		this->veh = Vehicle();
		std::cout << "Vehicle Length: " << this->veh.length << std::endl;

		this->start = nullptr;
		this->goal = nullptr;
		this->res = false;
		this->traj.resize(1, 9);

		this->state_list = std::vector<StatePtr>();
		// this->size_list = std::vector<int>();
		this->traj_list = std::vector<ArrayXXd*>();

		this->pq = PQ();
		this->state_dict = State_Dict();
		this->traj_dict = Traj_Dict();
	}

	Simulation::~Simulation()
	{
		sqlite3_close(this->db);
	}

	bool Simulation::set_boundary_conditions(double x1, double y1, double theta1, double k1, double x2, double y2, double theta2, double k2, double v1, double v2, double a1)
	{
		return false;
	}

	bool Simulation::set_boundary_conditions(double r_s1, double r_l1, double r_s2, double r_l2, double v1, double v2, double a1)
	{
		bool res = false;
		if (this->mode == MODE::OnRoadDynamic || this->mode == MODE::OnRoadStatic)
		{
			this->start = std::make_shared<State>(r_s1, r_l1, &(this->road), v1, a1, 0.);
			this->goal = std::make_shared<State>(r_s2, r_l2, &(this->road), v2);
			this->state_dict[std::make_tuple(this->start->r_i, this->start->r_j, this->start->v_i)] = this->start;
			this->state_dict[std::make_tuple(this->goal->r_i, this->goal->r_j, this->goal->v_i)] = this->goal;
			std::cout << "StateDict Size: " << this->state_dict.size() << std::endl;
			std::cout << "Start(SL): " << this->start->r_s << "\t" << this->start->r_l << std::endl;
			std::cout << "Goal(SL): " << this->goal->r_s << "\t" << this->goal->r_l << std::endl;
			this->pq.push(this->start);
			std::cout << "PQ Size: " << this->pq.size() << std::endl;
			res = true;
		}
		return res;
	}

	bool Simulation::set_boundary_conditions(int r_i1, int r_j1, int r_i2, int r_j2, int v_i1, int v_i2, double a1)
	{
		bool res = false;
		if (this->mode == MODE::OnRoadDynamic || this->mode == MODE::OnRoadStatic)
		{
			this->start = std::make_shared<State>(r_i1, r_j1, &(this->road), v_i1*2., a1, 0.);
			this->goal = std::make_shared<State>(r_i2, r_j2, &(this->road), v_i2*2.);
			this->state_dict[std::make_tuple(this->start->r_i, this->start->r_j, this->start->v_i)] = this->start;
			this->state_dict[std::make_tuple(this->goal->r_i, this->goal->r_j, this->goal->v_i)] = this->goal;
			std::cout << "StateDict Size: " << this->state_dict.size() << std::endl;
			std::cout << "Start(SL): " << this->start->r_s << "\t" << this->start->r_l << std::endl;
			std::cout << "Goal(SL): " << this->goal->r_s << "\t" << this->goal->r_l << std::endl;
			this->pq.push(this->start);
			std::cout << "PQ Size: " << this->pq.size() << std::endl;
			res = true;
		}
		return res;
	}

	bool Simulation::run()
	{
		std::cout << "Start Astar...... \n";
		this->res = Astar(this->pq, this->goal, this->state_dict, this->traj_dict, &(this->road), this->veh, this->cost_map, this->heuristic_map, this->db);
		std::cout << "Planning Result: " << this->res << std::endl;
		if (res)
		{
			// this->state_list, this->res, this-> traj
			StatePtr state = this->goal;
			this->state_list.push_back(state);
			int N = 0;
			while (state->parent != nullptr)
			{
				this->traj_list.push_back(&(this->traj_dict[std::make_tuple(state->parent, state)]));
				N += this->traj_list.back()->rows();
				this->state_list.push_back(state->parent);
				state = state->parent;
			}
			std::reverse(this->state_list.begin(), this->state_list.end());
			std::reverse(this->traj_list.begin(), this->traj_list.end());
			// rows of this->traj
			N -= this->traj_list.size() - 1;
			// this->traj
			this->traj.resize(N, 9);
			int anchor = 0;
			for (auto it = this->traj_list.cbegin(); it != this->traj_list.cend(); ++it)
			{
				this->traj.block(anchor, 0, (*it)->rows(), 9) = ArrayXXd(**it);
				anchor += (*it)->rows() - 1;
			}
		}
		return this->res;
	}

	std::shared_ptr<ArrayXXd> Simulation::interpolate(double time) const
	{
		std::shared_ptr<ArrayXXd> state = nullptr;
		if (!(time < 0.) && !(time > this->goal->time))
		{
			state = std::make_shared<ArrayXXd>(1, 9);
			int i = 1;
			while (time > this->state_list[i]->time)
				i++;
			i--; // traj_list
			// interpolate
			int j = 1;
			while (time > (*(this->traj_list[i]))(j, 0))
				j++;
			double dt = (*(this->traj_list[i]))(j, 0) - (*(this->traj_list[i]))(j - 1, 0);
			double dt1 = time - (*(this->traj_list[i]))(j - 1, 0);
			double dt2 = (*(this->traj_list[i]))(j, 0) - time;
			*state = ((*(this->traj_list[i])).row(j-1)*dt2 + (*(this->traj_list[i])).row(j)*dt1) / dt;
		}
		return state;
	}

}
