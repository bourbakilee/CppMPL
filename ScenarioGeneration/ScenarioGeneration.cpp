// read txt files into ArrayXXd
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

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