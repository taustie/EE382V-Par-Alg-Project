#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

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
	num = fabs(num);
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


void quick_hull(Point * points_array, int num_points){
	Point first_points_array[num_points];
	Point second_points_array[num_points];

	Point min_point = points_array[0];
	Point max_point = points_array[0];
	for(int i = 0; i < num_points; i++){ // to do: in parallel!!!
		if(min_point.y > points_array[i].y){
			copy_point(&min_point, &points_array[i]);
		}
		if(max_point.y < points_array[i].y){
			copy_point(&max_point, &points_array[i]);
		}
		copy_point(&first_points_array[i], &points_array[i]);
		copy_point(&second_points_array[i], &points_array[i]);
	}

	Set_of_Points left_set;
	Set_of_Points right_set;
	left_set.points_array = first_points_array;
	left_set.num_points = num_points;
	right_set.points_array = second_points_array;
	right_set.num_points = num_points;

	// can reduce space complexity by freeing the array in first call to subhull
	printf("min_point: (%lf, %lf)\n", min_point.x, min_point.y);
	printf("max_point: (%lf, %lf)\n", max_point.x, max_point.y);
	printf("leaving quick_hull\n");
	Set_of_Points* left_hull = sub_hull(left_set, min_point, max_point);
	Set_of_Points* right_hull = sub_hull(right_set, max_point, min_point);
	int total_points = left_hull->num_points + right_hull->num_points;
	Point convex_hull[total_points];
	for(int i = 0; i < total_points; i++){ // to do: in parallel!!!
		int a = i;
		if(i >= left_hull->num_points){
			a = i - left_hull->num_points;
			convex_hull[i] = right_hull->points_array[a];
		}
		else{
			convex_hull[i] = left_hull->points_array[a];
		}
	}

	// view result
	printf("total_points = %d\n", total_points);
	for(int i = 0; i < total_points; i++){
		printf("point %d: (%lf, %lf)\n", i, convex_hull[i].x, convex_hull[i].y);
	}
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

Set_of_Points * sub_hull(Set_of_Points all_points, Point p1, Point p2){
	int num_points = all_points.num_points;
	int left_count = 0;


	Point point_for_found_one;
	Vector line_vector = form_zero_vector(p1, p2);
	for(int i = 0; i < num_points; i++){
		Vector point_vector = form_zero_vector(p2, all_points.points_array[i]);
		if(cross_prod_orientation(line_vector, point_vector)){
			left_count++;
			printf("found left point: (%lf, %lf)\n", all_points.points_array[i].x, all_points.points_array[i].y);
			copy_point(&point_for_found_one, &all_points.points_array[i]);
		}
	}

	printf("p1: (%lf, %lf)\n", p1.x, p1.y);
	printf("p2: (%lf, %lf)\n", p2.x, p2.y);
	printf("left_count: %d\n", left_count);

	Point *left_points;
	Set_of_Points *result;
	Set_of_Points subset;
	if(left_count == 0){
		left_points = (Point *) malloc(sizeof(Point) * (left_count + 1));
		copy_point(&left_points[0], &p1);
		result = (Set_of_Points *) malloc(sizeof(Set_of_Points));
		result->num_points = left_count + 1;
		result->points_array = left_points;
		printf("returning result point: (%lf, %lf)\n", result->points_array[0].x, result->points_array[0].y);
		return result;
	}
	else if(left_count == 1){
		left_points = (Point *) malloc(sizeof(Point) * (left_count + 1));
		copy_point(&left_points[0], &p1);
		copy_point(&left_points[1], &point_for_found_one);
		printf("other point: (%lf, %lf)\n", left_points[1].x, left_points[1].y);
		result = (Set_of_Points *) malloc(sizeof(Set_of_Points));
		result->num_points = left_count + 1;
		result->points_array = left_points;
		return result;
	}
	else{
		left_points = (Point *) malloc(sizeof(Point) * (left_count));
		subset.points_array = left_points;
		subset.num_points = left_count;
	}

	int count = 0;
	for(int i = 0; i < num_points; i++){
		Vector point_vector = form_zero_vector(p2, all_points.points_array[i]);
		if(cross_prod_orientation(line_vector, point_vector)){
			copy_point(&left_points[count], &all_points.points_array[i]);
			count++;
		}
	}

	Line line = get_line(p1, p2);
	double max_dist = 0;
	Point max_point;
	for(int i = 0; i < left_count; i++){
		double distance = get_distance(line, left_points[i]);
		if(max_dist < distance){
			max_dist = distance;
			copy_point(&max_point, &left_points[i]);
		}
	}

	Set_of_Points* left_hull = sub_hull(subset, p1, max_point);
	Set_of_Points* right_hull = sub_hull(subset, max_point, p2);
	int total_points = left_hull->num_points + right_hull->num_points;


	Point * merged_points = (Point *) malloc(sizeof(Point) * total_points);
	result = (Set_of_Points *) malloc(sizeof(Set_of_Points));
	result->num_points = total_points;
	result->points_array = merged_points;

	for(int i = 0; i < total_points; i++){ // to do: in parallel!!!
		int a = i;
		if(i >= left_hull->num_points){
			a = i - left_hull->num_points;
			merged_points[i] = right_hull->points_array[a];
		}
		else{
			merged_points[i] = left_hull->points_array[a];
		}
	}

	return result;
}

void test_case_1(void){

	std::vector<Point *> input_points;
	std::list<Point *> output_points;

	Point * tmp = new Point; tmp->x = 0; tmp->y = 0; input_points.push_back(tmp);
	tmp = new Point; tmp->x = -1; tmp->y = 1; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 1; tmp->y = 1; input_points.push_back(tmp);
	tmp = new Point; tmp->x = 0; tmp->y = 2; input_points.push_back(tmp);

	quick_hull_new(input_points, output_points);

	std::cout << "Convex Hull has: " << output_points.size() << " points" << std::endl;
	std::list<Point *>::iterator it;
	for (it = output_points.begin(); it != output_points.end(); ++it){
		std::cout << "Found point: (" << (*it)->x << "," << (*it)->y << ")" << std::endl;
	}

	// Point input[9];
	// input[0].x = 0; input[0].y = 0;
	// input[1].x = -1; input[1].y = 1;
	// input[2].x = 1; input[2].y = 1;
	// input[3].x = 0; input[3].y = 2;
	//
	// input[4].x = -0.25; input[4].y = 0.25;
	// input[5].x = 0.25; input[5].y = 0.25;
	// input[6].x = -0.25; input[6].y = 0.75;
	// input[7].x = 0.25; input[7].y = 0.75;
	//
	// input[8].x = 0.25; input[8].y = 2;
	// quick_hull(input, 9);
}

int main() {
	// call a function in another file
	test_case_1();

	// quick_hull(input, 9);

	printf("hi\n");
	return(0);
}
