#include <cmath>
#include <chrono>

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
		// std::cout << "result =  " << result  << std::endl;
		return true;
	}
	else{
		return false;
	}
}

void quick_hull_new(std::vector<Point *> &input_points, std::list<Point *> &convex_hull){
	// Duplicate input points for each separate call to sub_hull -> O(n)
	std::vector<Point *> duplicate_points_1(input_points);
	std::vector<Point *> duplicate_points_2(input_points);

	// Identify min and max -> O(n)
	Point* min_point = input_points.at(0);
	Point* max_point = input_points.at(0);

	#ifdef QUICKHULL_PARALLEL
	get_max_min_points_parallel(input_points, &max_point, &min_point);
	#else
	for(std::vector<Point *>::iterator it = input_points.begin(); it != input_points.end(); ++it){
		// because colinear points are not added to the convex hull must choose min and max more carefully
		if(min_point->y > (*it)->y){ // min y point such that x is also minimal
			min_point = *it;
		}
		else if(min_point->y == (*it)->y){
			if(min_point->x > (*it)->x){
				min_point = *it;
			}
		}

		if(max_point->y < (*it)->y){ // max y point such that x is also maximal
			max_point = *it;
		}
		else if(max_point->y == (*it)->y){
			if(max_point->x < (*it)->x){
				max_point = *it;
			}
		}
	}
	#endif

	std::cout << "min point: (" << min_point->x << "," << min_point->y << ")" << std::endl;
	std::cout << "max point: (" << max_point->x << "," << max_point->y << ")" << std::endl;

	#ifdef QUICKHULL_PARALLEL
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			#endif
			sub_hull_new(duplicate_points_1, min_point, max_point, convex_hull);
			#ifdef QUICKHULL_PARALLEL
		}
		#pragma omp section
		{
			#endif
			sub_hull_new(duplicate_points_2, max_point, min_point, convex_hull);
			#ifdef QUICKHULL_PARALLEL
		}
	}
	#endif

}

void sub_hull_new(std::vector<Point *> &input_points, Point* p1, Point* p2, std::list<Point *> &convex_hull){
	// Modify input points to only be those to the left of the line
	std::vector<Point *> left_points; // must copy (in parallel each subhull call will clear the other's list)
	Vector ref_vector = form_zero_vector(*p1, *p2);

	// Sequential version is faster...Why is this???
	// #ifdef QUICKHULL_PARALLEL
	// #pragma omp parallel for schedule(static)
	// #endif
	for(long long int i = 0; i < input_points.size(); i++){
		Vector test_vector = form_zero_vector(*p2, *input_points.at(i));
        if (cross_prod_orientation(ref_vector, test_vector)) {
        	// #ifdef QUICKHULL_PARALLEL
        	 //#pragma omp critical
        	// #endif
			{
				left_points.push_back(input_points.at(i));
			}
        }
	}

	// See great resource at: https://www.py4u.net/discuss/63446

	// Option 1: use merge reduction (same as option 2)
	// #pragma omp declare reduction (merge : std::vector<Point *> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
	// #pragma omp parallel for reduction(merge: left_points)
	// for(long long int i=0; i < input_points.size(); i++) {
	// 	Vector test_vector = form_zero_vector(*p2, *input_points.at(i));
    //     if (cross_prod_orientation(ref_vector, test_vector)) {
	// 		left_points.push_back(input_points.at(i));
	// 	}
	//
	// }

	// Option 2: use private copy and merge at end
	// #pragma omp parallel
	// {
	//     std::vector<Point *> vec_private;
	// 	Vector test_vector;
	//     #pragma omp for nowait //fill vec_private in parallel
	//     for(int i=0; i < input_points.size(); i++) {
	// 		test_vector = form_zero_vector(*p2, *input_points.at(i));
	// 		if (cross_prod_orientation(ref_vector, test_vector)) {
	// 			vec_private.push_back(input_points.at(i));
	// 		}
	//     }
	//     #pragma omp critical
	//     left_points.insert(left_points.end(), vec_private.begin(), vec_private.end());
	// }



	//std::cout << "left_points.size(): " << left_points.size() << std::endl;

	// Update the convex hull
	if(left_points.size() < 2){

		// must update convex_hull atomically
		if(left_points.size() == 1) {
			#ifdef QUICKHULL_PARALLEL
		    #pragma omp critical
		    #endif
			{
				//std::cout << "Add to hull point: (" << left_points.at(0)->x << "," << left_points.at(0)->y << ")" << std::endl;
				convex_hull.push_back(left_points.at(0));
			}

		}
		#ifdef QUICKHULL_PARALLEL
        #pragma omp critical
        #endif
		{
			// std::cout << "Add to hull point: (" << p1->x << "," << p1->y << ")" << std::endl;
			convex_hull.push_back(p1);
		}
	}
	else{


		#ifdef QUICKHULL_PARALLEL
			Point* max_point;
			Line line = get_line(*p1, *p2);

			Dist_Info maxDistResult;
			maxDistResult.index = 0;
			maxDistResult.line_dist = 0.0;
			maxDistResult.p1_dist = 0.0;
			long long int index;

			#pragma omp declare reduction \
	        (maxDist:Dist_Info:omp_out=Dist_Max_Compare(omp_out, omp_in)) \
	        initializer(omp_priv = neutral_distance())

		    #pragma omp parallel for reduction(maxDist:maxDistResult)
		    for(index = 0; index < left_points.size(); index++){
				double distance = get_distance(line, *left_points.at(index));
				if(maxDistResult.line_dist < distance){
					maxDistResult.index = index;
					maxDistResult.line_dist = distance;
					maxDistResult.p1_dist = get_distance_point(*left_points.at(index), *p1);
				}
				else if(maxDistResult.line_dist == distance){
					double distance_p1 = get_distance_point(*left_points.at(index), *p1);
					if(distance_p1 < maxDistResult.p1_dist){
						maxDistResult.index = index;
						maxDistResult.p1_dist = distance_p1;
					}
				}
			}
			max_point = left_points.at(maxDistResult.index);

		#else
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
					std::cout << "found identical distance with point: x = " << left_points.at(i)->x << ", y = " << left_points.at(i)->y << std::endl;
					std::cout << "line is a: " << line.a << ", b: " << line.b << ", c: " << line.c << std::endl;
				}
			}
		#endif


		#ifdef QUICKHULL_PARALLEL
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				#endif
				sub_hull_new(left_points, p1, max_point, convex_hull);
				#ifdef QUICKHULL_PARALLEL
			}
			#pragma omp section
			{
				#endif
				sub_hull_new(left_points, max_point, p2, convex_hull);
				#ifdef QUICKHULL_PARALLEL
			}
		}
		#endif


		// #pragma omp parallel default(none) shared(left_points, p1, p2, max_point, convex_hull)
		// {
		// 	#pragma omp single
		// 	{
		// 		#pragma omp task //shared(left_points, p1, max_point, convex_hull)
		// 		sub_hull_new(left_points, p1, max_point, convex_hull);
		// 		#pragma omp task //shared(left_points, p2, max_point, convex_hull)
		// 		sub_hull_new(left_points, max_point, p2, convex_hull);
		// 	}
		// }

	}
}
