#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>

#include "quick_hull.hpp"
#include "parallel_funcs.hpp"
#include "verify_conv_hull.hpp"
#include "test_config.hpp"

int my_rand_max = 100000000;

void write_input(std::vector<Point *> &input_points){
	std::ofstream input_points_file;
    input_points_file.open("../visualize/input.txt");
	for(int i = 0; i < input_points.size(); i++){
		input_points_file << "(" << input_points.at(i)->x << ", "<< input_points.at(i)->y << ")" << std::endl;
	}
    input_points_file.close();
}

void write_output(std::list<Point *> &output_points){
	// output_points.sort(compare_angles);
	std::list<Point *> output_copy(output_points);
	find_min_and_sort(output_copy);
	std::ofstream output_points_file;
    output_points_file.open("../visualize/output.txt");
	for (std::list<Point *>::iterator it = output_copy.begin(); it != output_copy.end(); ++it){
		output_points_file << "(" << (*it)->x << ", "<< (*it)->y << ")" << std::endl;
	}
	output_points_file << "(" << output_copy.front()->x << ", "<< output_copy.front()->y << ")" << std::endl;
    output_points_file.close();
}

// manually assign input points and verify using visualizer
void manual_test(void){
	std::cout << "QuickHull manual test 1:" << std::endl;
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	Point * tmp = new Point;
	tmp->x = 0; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 5; tmp->y = 5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 3; tmp->y = 4; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 6; input_points.push_back(tmp);

	write_input(input_points);
	quick_hull(input_points, output_points);
	write_output(output_points);

	verify_convex_hull(input_points, output_points);
	std::cout << "\tConvex Hull has: " << output_points.size() << " points" << std::endl;
	for(int i = 0; i < input_points.size(); i++){
		delete(input_points.at(i));
	}
}

// Used to verify correctness of convex hull output points at scale
void unit_test(){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	srand(time(NULL));
	int element_count = 10000000;
	std::cout << "QuickHull unit test:" << std::endl;
	std::cout << "\tnumber of inputs: " << element_count << std::endl;
	std::cout << "\tusing random values between 0 and " << my_rand_max << std::endl;
	for(int i = 0; i < element_count; i++){
		Point * tmp = new Point;
		int value = (rand() % my_rand_max) - (my_rand_max/2); // less than int max to reduce the chance of floating point rounding errors
		tmp->x = value;
		value = (rand() % my_rand_max) - (my_rand_max/2);
		tmp->y = value;
		input_points.push_back(tmp);
	}

	write_input(input_points);
	auto t1 = std::chrono::high_resolution_clock::now();
	quick_hull(input_points, output_points);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	std::chrono::duration<double> ms_double = t2 - t1;
	std::cout << "\tConvex Hull has: " << output_points.size() << " points" << std::endl;
	std::cout << "\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;
	write_output(output_points);
	verify_convex_hull(input_points, output_points);
	for(int i = 0; i < input_points.size(); i++){
		delete(input_points.at(i));
	}
}

// Convex Hull should be a circular polygon
void circular_polygon(){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	srand(time(NULL));
	int element_count = 1000000;
	std::cout << "Circular Polygon:" << std::endl;
	std::cout << "\tnumber of inputs: " << element_count << std::endl;
	std::cout << "\tusing random values between 0 and " << my_rand_max << std::endl;
	for(int i = 0; i < element_count; i++){
		Point * tmp = new Point;
		int value = (rand() % (my_rand_max/2));
		tmp->x = value;
		// Technique to prevent overflow: https://www.johndcook.com/blog/2010/06/02/whats-so-hard-about-finding-a-hypotenuse/
		double max = my_rand_max/2 + 1;
		double min = value;
		double r = min / max;
		tmp->y = max * sqrt(1 - r*r);
		tmp->y -= rand() % (tmp->y + 1);
		if((value % 2) == 0){
			tmp->y *= -1;
		}
		if((rand() % 2) == 0){
			tmp->x *= -1;
		}

		input_points.push_back(tmp);
	}

	write_input(input_points);
	auto t1 = std::chrono::high_resolution_clock::now();
	quick_hull(input_points, output_points);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	std::chrono::duration<double> ms_double = t2 - t1;
	std::cout << "\tConvex Hull has: " << output_points.size() << " points" << std::endl;
	std::cout << "\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;
	write_output(output_points);
	verify_convex_hull(input_points, output_points);
	for(int i = 0; i < input_points.size(); i++){
		delete(input_points.at(i));
	}
}

// Displays the execution time to compute convex hull with exponential sweep on input sizes
void benchmark_square_exp_test(void){
	std::cout << "Benchmark Square Exponential Sweep Test:" << std::endl;
	int repeat_count = 1;
	std::cout << "\trepeat_count: " << repeat_count << std::endl;
	for(int element_count = 100; element_count <= 100000000; element_count*=10){
		std::cout << "\tNumber of inputs: " << element_count << std::endl;
		for(int i = 0; i < repeat_count; i++){
			std::vector<Point *> input_points;
			std::list<Point *> output_points;

			srand(time(NULL));
			for(int i = 0; i < element_count; i++){
				Point * tmp = new Point;
				int value = (rand() % my_rand_max) - (my_rand_max/2); // less than int max to reduce the chance of floating point rounding errors
				tmp->x = value;
				value = (rand() % my_rand_max) - (my_rand_max/2);
				tmp->y = value;
				input_points.push_back(tmp);
			}

			auto t1 = std::chrono::high_resolution_clock::now();
			quick_hull(input_points, output_points);
			auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::chrono::duration<double> ms_double = t2 - t1;
			std::cout << "\t\tConvex Hull has: " << output_points.size() << " points" << std::endl;
			std::cout << "\t\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;

			verify_convex_hull(input_points, output_points);
			for(int i = 0; i < input_points.size(); i++){
				delete(input_points.at(i));
			}
		}
	}
}

// Displays the execution time to compute convex hull with exponential sweep on input sizes
void benchmark_circular_exp_test(void){
	std::cout << "Benchmark Circular Exponential Sweep Test:" << std::endl;
	int repeat_count = 1;
	std::cout << "\trepeat_count: " << repeat_count << std::endl;
	for(int element_count = 100; element_count <= 100000000; element_count*=10){
		std::cout << "\tNumber of inputs: " << element_count << std::endl;
		for(int i = 0; i < repeat_count; i++){
			std::vector<Point *> input_points;
			std::list<Point *> output_points;

			srand(time(NULL));
			for(int i = 0; i < element_count; i++){
				Point * tmp = new Point;
				int value = (rand() % (my_rand_max/2));
				tmp->x = value;
				// Technique to prevent overflow: https://www.johndcook.com/blog/2010/06/02/whats-so-hard-about-finding-a-hypotenuse/
				double max = my_rand_max/2 + 1;
				double min = value;
				double r = min / max;
				tmp->y = max * sqrt(1 - r*r);
				tmp->y -= rand() % (tmp->y + 1);
				if((value % 2) == 0){
					tmp->y *= -1;
				}
				if((rand() % 2) == 0){
					tmp->x *= -1;
				}
				input_points.push_back(tmp);
			}

			auto t1 = std::chrono::high_resolution_clock::now();
			quick_hull(input_points, output_points);
			auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::chrono::duration<double> ms_double = t2 - t1;
			std::cout << "\t\tConvex Hull has: " << output_points.size() << " points" << std::endl;
			std::cout << "\t\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;

			verify_convex_hull(input_points, output_points);
			for(int i = 0; i < input_points.size(); i++){
				delete(input_points.at(i));
			}
		}
	}
}

// Displays the execution time to compute convex hull with linear sweep on input sizes
void benchmark_square_linear_test(void){
	std::cout << "Benchmark Square Linear Sweep Test:" << std::endl;
	int repeat_count = 1;
	std::cout << "\trepeat_count: " << repeat_count << std::endl;
	for(int element_count = 10000000; element_count <= 100000000; element_count+=10000000){
		std::cout << "\tNumber of inputs: " << element_count << std::endl;
		for(int i = 0; i < repeat_count; i++){
			std::vector<Point *> input_points;
			std::list<Point *> output_points;

			srand(time(NULL));
			for(int i = 0; i < element_count; i++){
				Point * tmp = new Point;
				int value = (rand() % my_rand_max) - (my_rand_max/2); // less than int max to reduce the chance of floating point rounding errors
				tmp->x = value;
				value = (rand() % my_rand_max) - (my_rand_max/2);
				tmp->y = value;
				input_points.push_back(tmp);
			}

			auto t1 = std::chrono::high_resolution_clock::now();
			quick_hull(input_points, output_points);
			auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::chrono::duration<double> ms_double = t2 - t1;
			std::cout << "\t\tConvex Hull has: " << output_points.size() << " points" << std::endl;
			std::cout << "\t\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;

			verify_convex_hull(input_points, output_points);
			for(int i = 0; i < input_points.size(); i++){
				delete(input_points.at(i));
			}
		}
	}
}

// Displays the execution time to compute convex hull with linear sweep on input sizes
void benchmark_circular_linear_test(void){
	std::cout << "Benchmark Circular Linear Sweep Test:" << std::endl;
	int repeat_count = 1;
	std::cout << "\trepeat_count: " << repeat_count << std::endl;
	for(int element_count = 10000000; element_count <= 100000000; element_count+=10000000){
		std::cout << "\tNumber of inputs: " << element_count << std::endl;
		for(int i = 0; i < repeat_count; i++){
			std::vector<Point *> input_points;
			std::list<Point *> output_points;

			srand(time(NULL));
			for(int i = 0; i < element_count; i++){
				Point * tmp = new Point;
				int value = (rand() % (my_rand_max/2));
				tmp->x = value;
				// Technique to prevent overflow: https://www.johndcook.com/blog/2010/06/02/whats-so-hard-about-finding-a-hypotenuse/
				double max = my_rand_max/2 + 1;
				double min = value;
				double r = min / max;
				tmp->y = max * sqrt(1 - r*r);
				tmp->y -= rand() % (tmp->y + 1);
				if((value % 2) == 0){
					tmp->y *= -1;
				}
				if((rand() % 2) == 0){
					tmp->x *= -1;
				}
				input_points.push_back(tmp);
			}

			auto t1 = std::chrono::high_resolution_clock::now();
			quick_hull(input_points, output_points);
			auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::chrono::duration<double> ms_double = t2 - t1;
			std::cout << "\t\tConvex Hull has: " << output_points.size() << " points" << std::endl;
			std::cout << "\t\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;

			verify_convex_hull(input_points, output_points);
			for(int i = 0; i < input_points.size(); i++){
				delete(input_points.at(i));
			}
		}
	}
}

void poor_mans_profile(){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	srand(time(NULL));
	int element_count = 20000000;
	std::cout << "Poor Man's Profiling:" << std::endl;
	std::cout << "\tnumber of inputs: " << element_count << std::endl;
	std::cout << "\tusing random values between 0 and " << my_rand_max << std::endl;
	for(int i = 0; i < element_count; i++){
		Point * tmp = new Point;
		int value = (rand() % (my_rand_max/2));
		tmp->x = value;
		// Technique to prevent overflow: https://www.johndcook.com/blog/2010/06/02/whats-so-hard-about-finding-a-hypotenuse/
		double max = my_rand_max/2 + 1;
		double min = value;
		double r = min / max;
		tmp->y = max * sqrt(1 - r*r);
		tmp->y -= rand() % (tmp->y + 1);
		if((value % 2) == 0){
			tmp->y *= -1;
		}
		if((rand() % 2) == 0){
			tmp->x *= -1;
		}

		input_points.push_back(tmp);
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	quick_hull(input_points, output_points);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	std::chrono::duration<double> ms_double = t2 - t1;
	std::cout << "\tConvex Hull has: " << output_points.size() << " points" << std::endl;
	std::cout << "\tExecution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;
	for(int i = 0; i < input_points.size(); i++){
		delete(input_points.at(i));
	}
}

int main() {
	// to visualize & verify algorithm correctness:
	// manual_test();
	unit_test();
	// circular_polygon();

	// to benchmark the algorithm:
	// benchmark_square_exp_test();
	// benchmark_circular_exp_test();
	// benchmark_square_linear_test();
	// benchmark_circular_linear_test();

	// to profile the execution via gdb stack captures
	// poor_mans_profile();
	return(0);
}
