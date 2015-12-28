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
		this->heuristic = inf;
	}



	State::State(State*prev, ArrayXXd & traj, Road * road, State * goal)
	{
	}

	double State::heuristic(State * goal)
	{
		return 0.0;
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
