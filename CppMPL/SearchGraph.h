#ifndef SEARCH_GRAPH_H
#define SEARCH_GRAPH_H

// define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
#include <TrajectoryGeneration.h>

namespace SearchGraph
{
	using namespace Eigen;
	using namespace environment;
	using namespace trajectory;

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

		// constructors
		// default constructor
		State();
		// goal state is used to calculate heuristic
		// use the state at the end of trajectory to construct new state

		// uesd to construct initial state, on-road only
		State(int i, int j, double velocity, double acc, Road* road, State* goal);
		// uesd to construct initial state, on-road only
		State(double r_s, double r_l, double velocity, double acc, Road* road, State* goal);
		// used to construct initial state, on-road and off-road
		State(double x, double y, double theta, double k, double v, double acc, State*goal, Road* road = nullptr);

		// used to construct next state, include the goal state, on-road
		State(State*prev, ArrayXXd& traj, int ref_rows, double traj_cost, Road* road, State* goal);

		// member functions

		// heuristic function
		double heuristic(State* goal);

		// non-member functions

		// compare
		friend bool operator< (const State& s1, const State& s2)
		{
			return s1.priority < s2.priority;
		}
		// heuristic function for heuristic seaarch algorithms
		friend double heuristic(State* s1, State* s2);
		// extend the current state
		// control set
		void extend();
	};

	// compute trajectory connect s1 and s2. if s2 is not reachable, s2 will be modified identical to the end state of trajectory.
	// if trajectory just not exists, return false
	bool connect(State* s1, State* s2, ArrayXXd& traj);


	/*----------------search graph---------------
	#include <boost/graph/adjacency_list.hpp> using namespace boost;
 typedef property<edge_weight_t, int> EdgeWeightProperty;
 typedef boost::adjacency_list<listS, vecS, directedS, no_property,
 EdgeWeightProperty > mygraph;
 int main() 
{ 
    mygraph g;
     add_edge (0, 1, 8, g);
     add_edge (0, 3, 18, g);
     add_edge (1, 2, 20, g);
     add_edge (2, 3, 2, g);
     add_edge (3, 1, 1, g);
     add_edge (1, 3, 7, g);
 }
	--------------------------------------*/
}

#endif
