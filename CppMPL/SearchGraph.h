/*
On Road State Lattice Builder
*/

#ifndef SEARCH_GRAPH_H
#define SEARCH_GRAPH_H

// #define EIGEN_USE_MKL_ALL
// #define WITH_SQLITE3

#include <Eigen/Dense>
#include <TrajectoryGeneration.h>

#include <cmath>
#include <tuple>
#include <unordered_map>
#include <queue>
#include <memory>

#include <boost/functional/hash.hpp>


namespace std
{
	
	// hash function for tuple
	template<typename... T>
	struct hash<tuple<T...>>
	{
		size_t operator()(tuple<T...> const& arg) const noexcept
		{
			return boost::hash_value(arg);
		}
	};

}

namespace SearchGraph
{
	using namespace Eigen;
	using namespace environment;
	using namespace trajectory;
	
	//
	const std::vector<double> accerations{ -4.,-2.,0.,2. };
	const std::vector<double> v_offsets{ -1.,-0.5,0.,0.5,1. };
	const std::vector<double> times{ 1.,2.,4. };
	// const std::vector<double> weights{ 5., 10., -0.1, 10., 0.1, 0.1, 50., 5, 40., -4. };

	struct State;
	struct StatePtrCompare;


	/*
	struct MyHash {
	size_t operator()(const Sstruct& x) const { return std::hash<int>()(x.Iv); }
    };

	std::unordered_map<Sstruct,ValueType,MyHash>
	*/

	/*
namespace std {
  template <>
  struct hash<Key>
  {
    std::size_t operator()(const Key& k) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return ((hash<string>()(k.first)
               ^ (hash<string>()(k.second) << 1)) >> 1)
               ^ (hash<int>()(k.third) << 1);
    }
  };

}
	*/




	using StatePtr = std::shared_ptr<State>;
	using State_Index = std::tuple<int, int, int>;
	// State_Dict: {(r_i,r_j,v_i): State}
	using State_Dict = std::unordered_map<State_Index, StatePtr>;

	using Time_State_Index = std::tuple<int, int, int, int>;
	// Time_State_Dict: {(t_i,r_i,r_j,v_i): State}
	using Time_State_Dict = std::unordered_map<Time_State_Index, StatePtr>;

	using Traj_Index = std::tuple<StatePtr, StatePtr>;
	// Traj_Dict: {(State, State): Traj}
	using Traj_Dict = std::unordered_map<Traj_Index, ArrayXXd>;

	//
	using VecSuccs = std::vector<StatePtr>;
	//
	using PQ = std::priority_queue<StatePtr, std::vector<StatePtr>, StatePtrCompare>;

	enum class MODE {OnRoadStatic, OnRoadDynamic, OffRoadStatic, OffRoadDynamic};


	// heuristic map
	struct HeuristicMap {
		bool isdynamic;
		double start_time;
		double end_time;
		double resolution;
		unsigned int cols;
		unsigned int rows;
		unsigned int num;
		std::vector<ArrayXXd> data; // vector<ArrayXXd> data
		double delta_t;

		// constructor
		HeuristicMap(const std::vector<ArrayXXd>& maps, CostMap& cost_map) :data(maps), isdynamic(cost_map.isdynamic), start_time(cost_map.start_time), end_time(cost_map.end_time), resolution(cost_map.resolution), cols(cost_map.cols), rows(cost_map.rows), delta_t(cost_map.delta_t), num(cost_map.num) {}

		// query heuristic
		double query(StatePtr current, const Vehicle& veh, StatePtr goal) const;
	};


	struct State
	{
		//members
		// MODE mode;
		// information about trajectory
		double time;
		int t_i;
		double length;
		double x; // bondary condition
		double y; // bondary condition
		double theta; // bondary condition
		double k; // bondary condition
		double dk;
		double v; // bondary condition
		int v_i;
		double a; // bondary condition, only for initial state
		// information about road
		double r_s;
		double r_l;
		int r_i;
		int r_j;
		// total cost of trajectory that it has passed
		double cost;
		double heuristic;
		// priority = cost + heuristic, used for priority queue
		double priority;
		//
		StatePtr parent; // default - nullptr
		bool reach; //default - false
		bool extend;

		// constructors
		// default constructor
		// State();
		// goal state is used to calculate heuristic
		// use the state at the end of trajectory to construct new state

		// uesd to construct initial state, on-road only
		State(double r_s, double r_l, Road* road, double velocity, double acc=0., double cost = inf, double dk = 0.);

		State(int r_i, int r_j, Road* road, double velocity, double acc, double cost = inf, double dk = 0.) :
			State(r_i*road->grid_length, r_j*road->grid_width, road, velocity, acc, cost, dk) {}

		State(double x, double y, double theta, double k, Road* road, double velocity, double acc, double cost = inf, double dk = 0.);

		// end state of traj
		State(const ArrayXXd& traj, Road* road) :State(traj(traj.rows() - 1, 2), traj(traj.rows() - 1, 3), traj(traj.rows() - 1, 4), traj(traj.rows() - 1, 5), road, traj(traj.rows() - 1, 7), traj(traj.rows() - 1, 8)) {}

		// member functions

		// extend the current state
		// control set
		void successors(VecSuccs& outs, State_Dict& state_dict, Road* road, StatePtr goal, const std::vector<double>& as = accerations, const std::vector<double>& vs = v_offsets, const std::vector<double>& ts = times, const double* p_lims = kinematic_limits);

		// void successors(VecSuccs& outs, Time_State_Dict& state_dict, Road* road, StatePtr goal, const std::vector<double>& as = accerations, const std::vector<double>& vs = v_offsets, const std::vector<double>& ts = times, const double* p_lims = kinematic_limits) {}


		// update cost, time, length and parent...
		static bool update(StatePtr current, StatePtr parent, double cost, const ArrayXXd& traj, Traj_Dict& traj_dict, const HeuristicMap& hm, StatePtr goal, const Vehicle& veh, double delta_t=0.1);

		static double distance(const StatePtr& s1, const StatePtr& s2) { return std::abs(s1->x - s2->x) + std::abs(s1->y - s2->y) + std::abs(s1->theta - s2->theta) + std::abs(s1->k - s2->k); }

		static double time_distance(const StatePtr& s1, const StatePtr& s2) { return std::abs(s1->time - s2->time) + std::abs(s1->x - s2->x) + std::abs(s1->y - s2->y) + std::abs(s1->theta - s2->theta) + std::abs(s1->k - s2->k); }

		// static members

		// process after trajectory evaluation
		static void post_process(StatePtr current, StatePtr successor, Eval_Res& res, ArrayXXd& traj, 
			PQ& pq, State_Dict& state_dict, Traj_Dict& traj_dict, StatePtr goal, 
			const Vehicle& veh, Road* road, const CostMap& cost_map, const HeuristicMap& hm, 
			sqlite3* db, const double* weights = cost_weights);

		// static void post_process(StatePtr current, StatePtr successor, Eval_Res& res, ArrayXXd& traj,
		//	PQ& pq, Time_State_Dict& state_dict, Traj_Dict& traj_dict, StatePtr goal,
		//	const Vehicle& veh, Road* road, const CostMap& cost_map, const HeuristicMap& hm,
		//	sqlite3* db, const double* weights = cost_weights);

		// non-member functions

		// compare
		
		friend bool operator< (const State& s1, const State& s2)
		{
			return s1.priority < s2.priority;
		}
		
	};


	struct StatePtrCompare
	{
		bool operator() (const StatePtr& lhs, const StatePtr&rhs) const
		{
			return rhs->priority < lhs->priority;
		}
	};


#ifdef WITH_SQLITE3
	// compute trajectory connect s1 and s2. if s2 is not reachable, s2 will be modified identical to the end state of trajectory.
	// if trajectory just not exists, return false
	bool connect(StatePtr s1, StatePtr s2, ArrayXXd& traj, sqlite3* db);
#endif

	/*
	Initialization:
	StatePtr start = std::make_shared<State>(r_s1, r_l1, road, v1, a1, 0.,-);
	StatePtr goal = std::make_shared<State>(r_s2, r_l2, road, v2, -, -,-);
	State_Dict state_dict;
	state_dict[ std::make_tuple(start.s_i, start.s_j, start.v_i) ] = start;
	state_dict[ std::make_tuple(goal.s_i, goal.s_j, goal.v_i) ] = goal;
	PQ pq;
	pq.push(start);
	Traj_Dict traj_dict;
	*/
	bool Astar(PQ& pq, StatePtr goal, State_Dict& state_dict, Traj_Dict& traj_dict, Road* road, const Vehicle& veh, CostMap& cost_map, HeuristicMap& hm, sqlite3* db);

	// bool Astar(StatePtr start, StatePtr goal, Time_State_Dict& state_dict, Traj_Dict& traj_dict, Road* road, const Vehicle& veh, CostMap& cost_map, HeuristicMap& hm, sqlite3* db);

}

#endif
