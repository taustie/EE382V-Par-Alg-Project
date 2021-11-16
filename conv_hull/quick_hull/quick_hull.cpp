#include <cmath>
#include <chrono>
#include <omp.h>

#include "quick_hull.hpp"
#include "parallel_funcs.hpp"
#include "verify_conv_hull.hpp"
#include "test_config.hpp"

// https://bobobobo.wordpress.com/2008/01/07/solving-linear-equations-ax-by-c-0/
Line get_line(Point p1, Point p2){
	Line result;
	result.a = p1.y - p2.y;
	result.b = p2.x - p1.x;
	result.c = (p1.x * p2.y) - (p2.x * p1.y);
	return result;
}

Vector form_zero_vector(Point p1, Point p2){
	Vector result;
	result.x = p2.x - p1.x;
	result.y = p2.y - p1.y;
	return result;
}

double get_distance_point(Point p1, Point p2){
	Vector v1 = form_zero_vector(p1, p2);
	double num = v1.x * v1.x + v1.y * v1.y;
	double result = sqrt(num);
	return result;
}

double get_distance(Line l1, Point p1){
	double num = (l1.a * p1.x) + (l1.b * p1.y) + l1.c;
	num = std::abs(num);
	double den = sqrt((l1.a * l1.a) + (l1.b * l1.b));
	double result = num / den;
	return result;
}

// true: right hand rule +z
// false: right hand rule -z
// Calculates (just z component of) v1 x v2
bool cross_prod_orientation(Vector v1, Vector v2){
	double result = v1.x * v2.y - v1.y * v2.x;
	if(result > 0){ // assumes colinear points are not distinct on convex hull
		return true;
	}
	else{
		return false;
	}
}

void quick_hull(std::vector<Point *> &input_points, std::list<Point *> &convex_hull){
	int threads = 2; // must be init to > 1
	omp_set_nested(1);
    if (omp_get_nested() != 1) {
        std::cout << "OMP_NESTED is FALSE!" << std::endl;
    }

	// Identify min and max
	Point* min_point = input_points.at(0);
	Point* max_point = input_points.at(0);

	get_max_min_points_parallel(input_points, &max_point, &min_point);

	if(threads > 1){
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				threads = omp_get_num_threads();
				#pragma omp task shared(convex_hull)
				{
					sub_hull_par(threads/2, input_points.data(), input_points.size(), min_point, max_point, convex_hull);
				}
				#pragma omp task shared(convex_hull)
				{
					sub_hull_par(threads/2, input_points.data(), input_points.size(), max_point, min_point, convex_hull);
				}
				std::cout << "\tthread count: " << threads << std::endl;
			} // don't have a double barrier - which is expensive
		}
	}
	else{
		std::cout << "\tthread count: 1" << std::endl;
		sub_hull_seq(input_points.data(), input_points.size(), min_point, max_point, convex_hull);
		sub_hull_seq(input_points.data(), input_points.size(), max_point, min_point, convex_hull);
	}
}

void get_points_on_left_sequential(std::vector<Point *> &left_points, Point ** input_points, int size, Point* p1, Point* p2){
	Vector ref_vector = form_zero_vector(*p1, *p2);
	for(long long int i = 0; i < size; i++){
		Vector test_vector = form_zero_vector(*p2, *input_points[i]);
		if (cross_prod_orientation(ref_vector, test_vector)) {
			left_points.push_back(input_points[i]);
		}
	}
}

void get_max_dist_sequential(std::vector<Point *> &left_points, Point** max, Point *p1, Point *p2){
	Line line = get_line(*p1, *p2);
	double max_dist = 0;
	double max_dist_p1 = 0;
	Point* max_point;
	for(long long int i = 0; i < left_points.size(); i++){
		double distance = get_distance(line, *left_points.at(i));
		if(max_dist < distance){
			max_dist = distance;
			max_dist_p1 = get_distance_point(*left_points.at(i), *p1);
			max_point = left_points.at(i);
		}
		else if(max_dist == distance){
			double value = get_distance_point(*left_points.at(i), *p1);
			if(value < max_dist_p1){ // either > or < works (just want to break ties with a middle colinear point)
				max_point = left_points.at(i);
				max_dist_p1 = value;
			}
		}
	}

	*max = max_point;
}

void sub_hull_seq(Point ** input_points, int size, Point* p1, Point* p2, std::list<Point *> &convex_hull){
	std::vector<Point *> left_points;
	get_points_on_left_sequential(left_points, input_points, size, p1, p2);

	if(left_points.size() < 2){ // 0 or 1 points on left, then do combine step
		if(left_points.size() == 1) {
		    #pragma omp critical
			convex_hull.push_back(left_points.at(0)); // Atomic update convex_hull
		}
        #pragma omp critical
		convex_hull.push_back(p1);
	}
	else{ // do divide step
		Point* max_point;
		get_max_dist_sequential(left_points, &max_point, p1, p2);

		sub_hull_seq(left_points.data(), left_points.size(), p1, max_point, convex_hull);
		sub_hull_seq(left_points.data(), left_points.size(), max_point, p2, convex_hull);
	}
}

void sub_hull_par(int threads, Point ** input_points, int size, Point* p1, Point* p2, std::list<Point *> &convex_hull){
	std::vector<Point *> left_points;
	get_points_on_left_parallel(left_points, input_points, size, p1, p2, threads);
	// get_points_on_left_filter_parallel(left_points, input_points, size, p1, p2);

	if(left_points.size() < 2){ // 0 or 1 points on left, then do combine step
		if(left_points.size() == 1) {
		    #pragma omp critical
			convex_hull.push_back(left_points.at(0)); // Atomic update convex_hull
		}
        #pragma omp critical
		convex_hull.push_back(p1);
	}
	else{ // do divide step
		Point* max_point;
		get_max_dist_parallel(left_points, &max_point, p1, p2, threads);

		if(threads > 1){
			#pragma omp task shared(convex_hull) //if(left_points.size() > 10)
			{
				sub_hull_par(threads/2, left_points.data(), left_points.size(), p1, max_point, convex_hull);
			}
			#pragma omp task shared(convex_hull) //if(left_points.size() > 10)
			{
				sub_hull_par(threads/2, left_points.data(), left_points.size(), max_point, p2, convex_hull);
			}
			#pragma taskwait
		}
		else{
			sub_hull_seq(left_points.data(), left_points.size(), p1, max_point, convex_hull);
			sub_hull_seq(left_points.data(), left_points.size(), max_point, p2, convex_hull);
		}
	}
}
