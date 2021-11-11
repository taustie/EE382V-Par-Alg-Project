#include <cmath>
#include "quick_hull.hpp"
#include "verify_conv_hull.hpp"

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

void find_min_and_sort(std::list<Point *> &output_points){
	// Find the minimum point, then sort points by angle from min point
	G_MIN_POINT = output_points.front();
	for (std::list<Point *>::iterator output_it = output_points.begin(); output_it != output_points.end(); ++output_it){
		if(G_MIN_POINT->y > (*output_it)->y){
			G_MIN_POINT = *output_it;
		}
		else if(G_MIN_POINT->y == (*output_it)->y){
			if(G_MIN_POINT->x > (*output_it)->x){
				G_MIN_POINT = *output_it;
			}
		}
	}

	output_points.sort(compare_angles);
}

bool verify_convex_hull(std::vector<Point *> &input_points, std::list<Point *> &output_points){
	// Find the minimum point, then sort points by angle from min point
	std::list<Point *> output_copy(output_points);
	find_min_and_sort(output_copy);
	
	// for (std::list<Point *>::iterator output_it = output_copy.begin(); output_it != output_copy.end(); ++output_it){
	// 	std::cout << "Point = (" << (*output_it)->x << "," << (*output_it)->y << ")" << std::endl;
	// }

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
