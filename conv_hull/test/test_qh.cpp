#include <cmath>
#include <chrono>

#include "quick_hull.hpp"
#include "parallel_funcs.hpp"
#include "verify_conv_hull.hpp"
#include "test_config.hpp"

void test_case_1(void){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	#if 0

	Point * tmp = new Point;
	tmp->x = 0; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 4; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -5; tmp->y = 6; input_points.push_back(tmp);

	tmp = new Point; tmp->x = 5; tmp->y = 5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 0; tmp->y = 10; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -1; tmp->y = 10; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 1; tmp->y = 10; input_points.push_back(tmp);

	// inside polygon
	tmp = new Point; tmp->x = 2; tmp->y = 3; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 2; tmp->y = 7; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -2; tmp->y = 3; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -2; tmp->y = 7; input_points.push_back(tmp);

	// test ray tracing
	tmp = new Point; tmp->x = 0; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 5; tmp->y = 5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 4; tmp->y = 7; input_points.push_back(tmp);

	tmp = new Point; tmp->x = 17766; tmp->y = 1191; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17766; tmp->y = 7834; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17766; tmp->y = 12999; input_points.push_back(tmp);

	tmp = new Point; tmp->x = -4596; tmp->y = -15000; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -3422; tmp->y = 17766; input_points.push_back(tmp);
	#endif

	Point * tmp = new Point;
	tmp->x = -4596; tmp->y = -15000; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 12673; tmp->y = -15000; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17037; tmp->y = -14997; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17633; tmp->y = -14981; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17693; tmp->y = -14969; input_points.push_back(tmp);

	tmp = new Point; tmp->x = 17709; tmp->y = -14854; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17761; tmp->y = -13763; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17764; tmp->y = -12724; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17766; tmp->y = 7834; input_points.push_back(tmp);

	tmp = new Point; tmp->x = 17766; tmp->y = 1191; input_points.push_back(tmp);

	// inside polygon

	tmp = new Point; tmp->x = 17766; tmp->y = 12999; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17764; tmp->y = 15911; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17757; tmp->y = 17228; input_points.push_back(tmp);

	// test ray tracing
	tmp = new Point; tmp->x = 17730; tmp->y = 17676; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 17559; tmp->y = 17747; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 15414; tmp->y = 17766; input_points.push_back(tmp);

	tmp = new Point; tmp->x = -3422; tmp->y = 17766; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -12900; tmp->y = 17765; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14703; tmp->y = 17763; input_points.push_back(tmp);

	tmp = new Point; tmp->x = -14916; tmp->y = 17756; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14984; tmp->y = 17642; input_points.push_back(tmp);

	tmp = new Point; tmp->x = -14996; tmp->y = 17160; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14999; tmp->y = 15393; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14999; tmp->y = -1822; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14997; tmp->y = -14400; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14982; tmp->y = -14798; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14960; tmp->y = -14838; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14710; tmp->y = -14922; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -14488; tmp->y = -14981; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -13890; tmp->y = -14996; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -9761; tmp->y = -14999; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -4596; tmp->y = -15000; input_points.push_back(tmp);

	quick_hull_new(input_points, output_points);

	std::cout << "Convex Hull has: " << output_points.size() << " points" << std::endl;
	std::list<Point *>::iterator it;
	for (it = output_points.begin(); it != output_points.end(); ++it){
		std::cout << "Convex Hull has point: (" << (*it)->x << "," << (*it)->y << ")" << std::endl;
	}

	// tmp = new Point; tmp->x = -3; tmp->y = 1; input_points.push_back(tmp);
	verify_convex_hull(input_points, output_points);
}

void test_case_2(void){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	srand(time(NULL));
	int element_count = 10000000;
	for(int i = 0; i < element_count; i++){
		Point * tmp = new Point;
		int value = rand() % 32767 - 15000;
		tmp->x = value;
		value = rand() % 32767 - 15000;
		tmp->y = value;
		input_points.push_back(tmp);
	}

	auto t1 = std::chrono::high_resolution_clock::now();
	quick_hull_new(input_points, output_points);
	auto t2 = std::chrono::high_resolution_clock::now();
	auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	std::chrono::duration<double> ms_double = t2 - t1;

	verify_convex_hull(input_points, output_points);
	std::cout << "Convex Hull has: " << output_points.size() << " points" << std::endl;
	std::cout << "Execution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;;
}

void test_case_3(void){


	#ifdef QUICKHULL_PARALLEL
	std::cout << "Parallel Test Case:" << std::endl;
	#else
	std::cout << "Sequential Test Case:" << std::endl;
	#endif

	int repeat_count = 1;
	for(int element_count = 1000; element_count <= 100000000; element_count*=10){
		for(int i = 0; i < repeat_count; i++){
			std::vector<Point *> input_points;
			std::list<Point *> output_points;

			srand(time(NULL));
			for(int i = 0; i < element_count; i++){
				Point * tmp = new Point;
				int value = rand() % 32767 - 15000;
				tmp->x = value;
				value = rand() % 32767 - 15000;
				tmp->y = value;
				input_points.push_back(tmp);
			}

			auto t1 = std::chrono::high_resolution_clock::now();
			quick_hull_new(input_points, output_points);
			auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::chrono::duration<double> ms_double = t2 - t1;

			verify_convex_hull(input_points, output_points);
			std::cout << "Convex Hull has: " << output_points.size() << " points" << std::endl;
			std::cout << "Execution time: " << ms_int.count() << " (ms), " <<  ms_double.count() << " (s)" << std::endl;;
			std::cout << std::endl;
		}
	}

}

int main() {
	//test_case_1();
	test_case_3();
	return(0);
}
