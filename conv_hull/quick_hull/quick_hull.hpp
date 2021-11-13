#ifndef QUICK_HULL
#define QUICK_HULL

#include <iostream>
#include <vector>
#include <list>

typedef struct point{
	long long int x;
	long long int y;
} Point;

typedef Point Vector;

typedef struct line{
	// assume line of the form ax+by+c = 0
	double a;
	double b;
	double c;
} Line; // https://brilliant.org/wiki/dot-product-distance-between-point-and-a-line/

Line get_line(Point p1, Point p2);
Vector form_zero_vector(Point p1, Point p2);
double get_distance_point(Point p1, Point p2);
double get_distance(Line l1, Point p1);
bool cross_prod_orientation(Vector v1, Vector v2);
void quick_hull(std::vector<Point *> &input_points, std::list<Point *> &convex_hull);
void sub_hull(Point ** input_points, int size, Point* p1, Point* p2, std::list<Point *> &convex_hull);

#endif
