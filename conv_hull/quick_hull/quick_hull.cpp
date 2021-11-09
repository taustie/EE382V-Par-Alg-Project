#include <cmath>
#include <chrono>

#include "quick_hull.hpp"
#include "parallel_funcs.hpp"

#define QUICKHULL_PARALLEL

double dot_product(Point p1, Point p2){
	double result = (p1.x * p1.x) + (p2.y + p2.y);
	return result;
}

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

void copy_point(Point *dst, Point *src){
	dst->x = src->x;
	dst->y = src->y;
}

void quick_hull_new(std::vector<Point *> &input_points, std::list<Point *> &convex_hull){
	// Duplicate input points for each separate call to sub_hull -> O(n)
	std::vector<Point *> duplicate_points_1(input_points);
	std::vector<Point *> duplicate_points_2(input_points);

	// Identify min and max -> O(n)
	Point* min_point = input_points.at(0);
	Point* max_point = input_points.at(0);

	#ifdef QUICKHULL_PARALLEL
	std::cout << "In parallel" << std::endl;
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



	std::cout << "left_points.size(): " << left_points.size() << std::endl;

	// Update the convex hull
	if(left_points.size() < 2){

		// must update convex_hull atomically
		if(left_points.size() == 1) {
			#ifdef QUICKHULL_PARALLEL
		    #pragma omp critical
		    #endif
			{
				std::cout << "Add to hull point: (" << left_points.at(0)->x << "," << left_points.at(0)->y << ")" << std::endl;
				convex_hull.push_back(left_points.at(0));
			}

		}
		#ifdef QUICKHULL_PARALLEL
        #pragma omp critical
        #endif
		{
			std::cout << "Add to hull point: (" << p1->x << "," << p1->y << ")" << std::endl;
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

Point *G_MIN_POINT;

// comparison, not case sensitive.
bool compare_angles (Point * first, Point * second){
	if(first == G_MIN_POINT){
		return true;
	}
	if(second == G_MIN_POINT){
		return false;
	}
	Vector first_vector = form_zero_vector(*G_MIN_POINT, *first);
	double first_angle;
	if(first_vector.x == 0){
		first_angle = M_PI/2;
	}
	else{
		first_angle = (float) first_vector.y / (float) first_vector.x;
		// std::cout << "first_angle = " << first_angle << std::endl;
		first_angle = atan(first_angle);
		if(first_angle < 0){
			first_angle = M_PI + first_angle;
		}
	}

	Vector second_vector = form_zero_vector(*G_MIN_POINT, *second);
	double second_angle;
	if(second_vector.x == 0){
		second_angle = M_PI/2;
	}
	else{
		second_angle = (float) second_vector.y / (float) second_vector.x;
		// std::cout << "second_angle = " << second_angle << std::endl;
		second_angle = atan(second_angle);
		if(second_angle < 0){
			second_angle = M_PI + second_angle;
		}
	}

	return (first_angle < second_angle);
}

bool verify_ccw_angles(std::vector<Point *> &input_points, std::list<Point *> &output_copy){
	std::list<Point *>::iterator first_it = output_copy.begin();
	std::list<Point *>::iterator second_it = output_copy.begin();
	std::advance(second_it, 1);
	std::list<Point *>::iterator third_it = output_copy.begin();
	std::advance(third_it, 2);
	while(third_it != output_copy.end()){
		Vector first_vector = form_zero_vector(**first_it, **second_it);
		Vector second_vector = form_zero_vector(**second_it, **third_it);
		if (!cross_prod_orientation(first_vector, second_vector)) {
			std::cout << "Point first_it = (" << (*first_it)->x << "," << (*first_it)->y << ")" << std::endl;
			std::cout << "Point second_it = (" << (*second_it)->x << "," << (*second_it)->y << ")" << std::endl;
			std::cout << "first_vector = (" << first_vector.x << "," << first_vector.y << ")" << std::endl;
			std::cout << "second_vector = (" << second_vector.x << "," << second_vector.y << ")" << std::endl;
			std::cout << "Fail 1: points on Convex Hull are not strictly counter-clock wise" << std::endl;
			return false;
        }
		++first_it;
		++second_it;
		++third_it;
	}

	// list is not circular so check the last two cases
	std::list<Point *>::reverse_iterator rev_it = output_copy.rbegin();
	Point p1 = **rev_it;
	std::advance(rev_it, 1);
	Point p0 = **rev_it;
	first_it = output_copy.begin();
	Point p2 = **first_it;
	std::advance(first_it, 1);
	Point p3 = **first_it;

	Vector first_vector = form_zero_vector(p0, p1);
	Vector second_vector = form_zero_vector(p1, p2);
	if (!cross_prod_orientation(first_vector, second_vector)) {
		std::cout << "Fail 2: points on Convex Hull are not strictly counter-clock wise" << std::endl;
		return false;
	}

	first_vector = form_zero_vector(p1, p2);
	second_vector = form_zero_vector(p2, p3);
	if (!cross_prod_orientation(first_vector, second_vector)) {
		std::cout << "Fail 3: points on Convex Hull are not strictly counter-clock wise" << std::endl;
		return false;
	}

	return true;
}

bool verify_convex_hull(std::vector<Point *> &input_points, std::list<Point *> &output_points){
	// Find the minimum point, then sort points by angle from min point
	std::list<Point *> output_copy(output_points);
	G_MIN_POINT = output_copy.front();
	for (std::list<Point *>::iterator output_it = output_copy.begin(); output_it != output_copy.end(); ++output_it){
		if(G_MIN_POINT->y > (*output_it)->y){
			G_MIN_POINT = *output_it;
		}
		else if(G_MIN_POINT->y == (*output_it)->y){
			if(G_MIN_POINT->x > (*output_it)->x){
				G_MIN_POINT = *output_it;
			}
		}
	}

	output_copy.sort(compare_angles);
	for (std::list<Point *>::iterator output_it = output_copy.begin(); output_it != output_copy.end(); ++output_it){
		std::cout << "Point = (" << (*output_it)->x << "," << (*output_it)->y << ")" << std::endl;
	}

	// verify all edges on convevx hull are counter clockwise
	verify_ccw_angles(input_points, output_copy);

	// Verify ray tracing
	for (std::vector<Point *>::iterator input_it = input_points.begin(); input_it != input_points.end(); ++input_it){
		int on_boundary_count = 0;
		int intersection_count = 0;
		bool skip_this_point = false;
		for (std::list<Point *>::iterator output_it = output_copy.begin(); output_it != output_copy.end(); ++output_it){
			if(((*output_it)->x == (*input_it)->x) && ((*output_it)->y == (*input_it)->y)){
				skip_this_point = true;
			}
			else if((*output_it)->y == (*input_it)->y){
				if((*input_it)->x < (*output_it)->x){
					intersection_count++;
				}
			}
		}

		if(skip_this_point){
			continue;
		}

		std::list<Point *>::iterator first_it = output_copy.begin();
		std::list<Point *>::iterator second_it = output_copy.begin();
		std::advance(second_it, 1);
		while(second_it != output_copy.end()){
			double max = (*first_it)->y;
			double min = (*first_it)->y;
			if(max < (*second_it)->y){
				max = (*second_it)->y;
			}
			if(min > (*second_it)->y){
				min = (*second_it)->y;
			}

			if (( (*input_it)->y < max ) && ( (*input_it)->y > min )){
				Line line_segment = get_line(**first_it, **second_it);
				double x_value_on_line = ((-1 * line_segment.b * (*input_it)->y) - line_segment.c) / line_segment.a;
				if((*input_it)->x < x_value_on_line){
					intersection_count++;
				}
				else if((*input_it)->x == x_value_on_line){
					on_boundary_count++;
				}
			}

			first_it++;
			second_it++;
		}

		std::list<Point *>::reverse_iterator rev_it = output_copy.rbegin();
		second_it = output_copy.begin();
		{
			double max = (*rev_it)->y;
			double min = (*rev_it)->y;
			if(max < (*second_it)->y){
				max = (*second_it)->y;
			}
			if(min > (*second_it)->y){
				min = (*second_it)->y;
			}

			if(((*input_it)->y < max) && ((*input_it)->y > min)){
				Line line_segment = get_line(**rev_it, **second_it);
				double x_value_on_line = ((-1 * line_segment.b * (*input_it)->y) - line_segment.c) / line_segment.a;
				if((*input_it)->x < x_value_on_line){
					intersection_count++;
				}
				else if((*input_it)->x == x_value_on_line){
					on_boundary_count++;
				}
			}
		}

		if( (intersection_count == 1 && on_boundary_count == 0) ||
			(intersection_count == 1 && on_boundary_count == 1) ||
			(intersection_count == 0 && on_boundary_count == 1) ){
				// passes
				continue;
			}
		else{
			std::cout << "Fail: found input point outside the convex polygon!" << std::endl;
			std::cout << "intersection_count = " << intersection_count << std::endl;
			std::cout << "on_boundary_count = " << on_boundary_count << std::endl;
			std::cout << "Point input_it = (" << (*input_it)->x << "," << (*input_it)->y << ")" << std::endl;
			return false;
		}
	}

	return true;
}

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

int main() {
	//test_case_1();
	test_case_2();
	return(0);
}
