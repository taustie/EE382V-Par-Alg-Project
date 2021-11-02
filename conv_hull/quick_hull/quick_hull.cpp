#include <cmath>
#include <iostream>
#include <vector>
#include <list>

#include "quick_hull.hpp"

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

double get_distance(Line l1, Point p1){
	double num = (l1.a * p1.x) + (l1.b * p1.y) + l1.c;
	num = std::abs(num);
	double den = sqrt((l1.a * l1.a) + (l1.b * l1.b));
	double result = num / den;
	return result;
}

Vector form_zero_vector(Point p1, Point p2){
	Vector result;
	result.x = p2.x - p1.x;
	result.y = p2.y - p1.y;
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
	for(std::vector<Point *>::iterator it = input_points.begin(); it != input_points.end(); ++it){
		if(min_point->y > (*it)->y){
			min_point = *it;
		}
		if(max_point->y < (*it)->y){
			max_point = *it;
		}
	}

	std::cout << "min point: (" << min_point->x << "," << min_point->y << ")" << std::endl;
	std::cout << "max point: (" << max_point->x << "," << max_point->y << ")" << std::endl;

	sub_hull_new(duplicate_points_1, min_point, max_point, convex_hull);
	sub_hull_new(duplicate_points_2, max_point, min_point, convex_hull);
}

void sub_hull_new(std::vector<Point *> &input_points, Point* p1, Point* p2, std::list<Point *> &convex_hull){
	// Modify input points to only be those to the left of the line
	std::vector<Point *> left_points; // must copy (in parallel each subhull call will clear the other's list)
	Vector ref_vector = form_zero_vector(*p1, *p2);
	for(std::vector<Point *>::iterator it = input_points.begin(); it != input_points.end(); ++it){
		Vector test_vector = form_zero_vector(*p2, **it);
        if (cross_prod_orientation(ref_vector, test_vector)) {
			left_points.push_back(*it);
        }
	}

	std::cout << "left_points.size(): " << left_points.size() << std::endl;

	// Update the convex hull
	if(left_points.size() < 2){
		std::cout << "adding p1: (" << p1->x << "," << p1->y << ")" << std::endl;
		// must update convex_hull atomically
		if(left_points.size() == 1) {convex_hull.push_back(left_points.at(0));}
		convex_hull.push_back(p1);
	}
	else{
		Line line = get_line(*p1, *p2);
		double max_dist = 0;
		Point* max_point;
		for(std::vector<Point *>::iterator it = left_points.begin(); it != left_points.end(); ++it){
			double distance = get_distance(line, **it);
			if(max_dist < distance){
				max_dist = distance;
				max_point = *it;
			}
		}

		sub_hull_new(left_points, p1, max_point, convex_hull);
		sub_hull_new(left_points, max_point, p2, convex_hull);
	}
}

void test_case_1(void){
	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	Point * tmp = new Point;
	tmp->x = 0; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -1; tmp->y = 1; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 1; tmp->y = 1; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 0; tmp->y = 2; input_points.push_back(tmp);

	// inside polygon
	tmp = new Point; tmp->x = 0.1; tmp->y = 0.5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 0.1; tmp->y = 1.5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -0.1; tmp->y = 0.5; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -0.1; tmp->y = 1.5; input_points.push_back(tmp);

	quick_hull_new(input_points, output_points);

	std::cout << "Convex Hull has: " << output_points.size() << " points" << std::endl;
	std::list<Point *>::iterator it;
	for (it = output_points.begin(); it != output_points.end(); ++it){
		std::cout << "Found point: (" << (*it)->x << "," << (*it)->y << ")" << std::endl;
	}
}

int main() {
	test_case_1();
	printf("hi\n");
	return(0);
}
