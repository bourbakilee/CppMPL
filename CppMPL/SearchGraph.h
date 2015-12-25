#ifndef SEARCH_GRAPH_H
#define SEARCH_GRAPH_H

// define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
#include <environment.h>

namespace SearchGraph
{
	using namespace Eigen;
	using namespace environment;

	struct State
	{
		//members
		// information about trajectory
		double time;
		double length;
		double x; // bondary condition
		double y; // bondary condition
		double theta; // bondary condition
		double k; // bondary condition
		double dk;
		double v; // bondary condition
		double a; // bondary condition, only for initial state
		// information about road
		double r_s;
		double r_l;
		int r_i;
		int r_j;
		// total cost of trajectory that it has passed
		double cost;
		// priority = cost + heuristic, used for priority queue
		double priority;

		// construct functions
		State() {}
		State(int i, int j, Road* road=nullptr, State* goal=nullptr);
		State(double r_s, double r_l, Road* road = nullptr, State* goal = nullptr);
		// use the state at the end of trajectory to construct new state
		State(const ArrayXXd& traj, Road* road = nullptr, State* goal = nullptr);

		// member functions
		// heuristic function
		double heuristic(const State& goal);

		// non-member functions
		// compare
		friend bool operator< (const State& s1, const State& s2)
		{
			return s1.priority < s2.priority;
		}
		// heuristic function for heuristic seaarch algorithms
		friend double heuristic(const State& s1, const State& s2);
	};

	// compute trajectory connect s1 and s2. if s2 is not reachable, s2 will be modified identical to the end state of trajectory.
	// if trajectory just not exists, return false
	bool connect(const State& s1, State& s2, ArrayXXd& traj);

}

#endif
