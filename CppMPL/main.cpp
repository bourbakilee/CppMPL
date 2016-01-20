#include "SearchGraph.h"
#include <iostream>

using namespace SearchGraph;

int main()
{
	// open database
	sqlite3* db = nullptr;
	int rc = sqlite3_open("InitialGuessTable.db", &db);
	if (rc) {
		std::cout << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
		sqlite3_close(db);
		return -1;
	}
	else {
		std::cout << "Opened database successfully" << std::endl;
	}

	// build environment

	//     vehicle
	Vehicle veh(2.94, 0.988, 1.28, 2.029);

	//     road
	//     read road center line from file
	char* center_line_file = "road_center_line.txt";
	ArrayXXd center_line;
	read_road_center_line(center_line_file, center_line);
	std::cout << "road center line: \n";
	std::cout << "Rows of center line array: "<<center_line.rows() << std::endl;
	Road road(center_line);
	std::cout << "Road constructed: " << road.length << std::endl;

	//     cost map
	char* cost_map_file = "cost_map.txt";
	ArrayXXd cost_map(500, 500);
	read_map(cost_map_file, cost_map);
	std::cout << "cost map: \n";
	std::cout << cost_map(499,499) << std::endl;
	CostMap cm(cost_map);

	//     heuristic map
	char* heuristic_map_file = "heuristic_map.txt";
	ArrayXXd heuristic_map;
	read_map(heuristic_map_file, heuristic_map);
	std::cout << "heuristic map: \n";
	std::cout << heuristic_map(499,499) << std::endl;
	std::vector<ArrayXXd> vec_hm{ heuristic_map };
	std::cout << vec_hm[0](499, 499) << std::endl;
	HeuristicMap hm(vec_hm, cm);
	std::cout << "heuristic map constructed duccessfully!\n";
	std::cout << hm.data[0](499, 499) << std::endl;

	// start, goal state
	StatePtr start = std::make_shared<State>(5., 1., &road, 8., 0., 0.);
	std::cout << "successfully construct start state:\n";
	std::cout << start->x << start->y << start->theta << start->k << start->v << std::endl;
	StatePtr goal = std::make_shared<State>(80., -1., &road, 9.);
	std::cout << "successfully construct goal state:\n";
	std::cout << goal->x << goal->y << goal->theta << goal->k << goal->v << std::endl;
	State_Dict state_dict;
	state_dict[std::make_tuple(start->r_i, start->r_j, start->v_i)] = start;
	state_dict[std::make_tuple(goal->r_i, goal->r_j, goal->v_i)] = goal;
	std::cout << "successfully construct state dict\n";
	std::cout << state_dict.size() << std::endl;
	PQ pq;
	pq.push(start);
	std::cout << "successfully construct Priority Queue\n";
	std::cout << pq.size() << std::endl;
	Traj_Dict traj_dict;
	std::cout << "successfully construct trajectory dict\n";
	std::cout << traj_dict.size() << std::endl;

	ArrayXXd traj;
	// search the graph
	// bool res = Astar(pq, goal, state_dict, traj_dict, &road, veh, cm, hm, db);
	bool res = connect(start, goal, traj, db);
	std::cout <<"compute trajectory:"<< res << std::endl;
	std::cout << "orignianl size of traj array: rows: " << traj.rows() << "  , cols: " << traj.cols() << std::endl;
	// std::cout << "Planning Successfully? :" << res << std::endl;

	if (res)
	{
		Eval_Res eval_res = eval_traj(traj, veh, cm, &road,true);
		std::cout << "size of traj array after eval: rows: " << traj.rows() << "  , cols: " << traj.cols() << std::endl;
		std::cout << "Cost: " << eval_res.first << std::endl;
		std::cout << "Truncated:  " << eval_res.second << std::endl;
	}


	sqlite3_close(db);

	return 0;
}