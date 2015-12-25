#include "SearchGraph.h"

SearchGraph::State::State(int i, int j, environment::Road * road, State* goal)
{
}

SearchGraph::State::State(double r_s, double r_l, environment::Road * road, State* goal)
{
}

SearchGraph::State::State(const ArrayXXd & traj, Road * road, State* goal)
{
}

double SearchGraph::State::heuristic(const State & goal)
{
	return 0.0;
}

double SearchGraph::heuristic(const State & s1, const State & s2)
{
	return 0.0;
}

bool SearchGraph::connect(const State & s1, State & s2, ArrayXXd & traj)
{
}
