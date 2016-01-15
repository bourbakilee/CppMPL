// read txt files into ArrayXXd
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;

void read_array(const char* filename, ArrayXXd& array)
{
	if (array.cols() != 500 || array.rows() != 500)
		array.resize(500, 500);
	//double data[500 * 500];
	ifstream input(filename);
	for (int i; i < 500 * 500; i++)
	{
		//input >> data[i];
		// if input number is biger than 1, it should be modified to inf
		input >> array(i/500, i%500);
	}
	//cout << array << endl;
}

void read_array(const char* filename, vector<vector<double>>& array)
{
	ifstream input(filename);
	int i = 0;
	double tmp;
	vector<double> line(5,0.);
	while (!(input >> tmp).fail()) {
		i++;
		line[i - 1] = tmp;
		if (i % 5 == 0)
		{
			array.push_back(line);
		}
	}
}

void read_array(const char* filename, vector<double>& array)
{
	ifstream input(filename);
	double tmp;
	while (!(input >> tmp).fail()) {
		
			array.push_back(tmp);
	}
	ArrayXXd road_center_line(array.data());
}


int main()
{
	char* file = "cost_grayscale_map.txt";
	ArrayXXd arr(500, 500);
	read_array(file, arr);
	cout << arr << endl;
	return 0;
}

/*
#include <iostream>
#include <fstream>

using namespace std;

int main() {

double data[size of your data];

std::ifstream input("file.txt");

for (int i = 0; i < size of your data; i++) {
input >> data[i];
std::cout<< data[i]<<std::endl;
}

}
*/