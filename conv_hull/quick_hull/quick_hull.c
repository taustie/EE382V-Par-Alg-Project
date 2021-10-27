#include <math.h>
#include <stdio.h>

typedef struct point{
	double x;
	double y;
} Point;

typedef Point Vector;

typedef struct line{
	// assume line of the form ax+by+c = 0
	double a;
	double b;
	double c;
} Line; // https://brilliant.org/wiki/dot-product-distance-between-point-and-a-line/

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
}

// true: right hand rule +z
// false: right hand rule -z
// Calculates (just z component of) v1 x v2
bool cross_prod_orientation(Vector v1, Vector v2){
	double result = v1.x * v2.y - v1.y * v2.x;
	if(result >= 0){ // assumes colinear points are distinct on convext hull
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

void quick_hull(Point * points_array, int num_points){
	Points first_points_array[num_points];
	Points second_points_array[num_points];

	Point min_point = points_array[0];
	Point max_point = points_array[0];
	for(int i = 0; i < num_points; i++){ // to do: in parallel!!!
		if(min_point.y > points_array[i].y){
			copy_point(&min_point, &points_array[i]);
		}
		if(max_point.y < points_array[i].y){
			copy_point(&max_point, &points_array[i]);
		}
		copy_point(&first_points_array[i], &points_array[i])
		copy_point(&second_points_array[i], &points_array[i])
	}

	// can reduce space complexity by freeing the array in first call to subhull
	sub_hull(first_points_array, num_points, min_point, max_point);
	sub_hull(second_points_array, num_points, max_point, min_point);
}

void sub_hull(Point * points_array, int num_points, Point p1, Point p2){

}

int main() {
	// call a function in another file
	printf("hi\n");
	return(0);
}
