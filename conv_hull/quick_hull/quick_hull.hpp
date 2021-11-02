#ifndef QUICK_HULL
#define QUICK_HULL

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

typedef struct set_of_points{
	Point * points_array;
	int num_points;
} Set_of_Points;

Set_of_Points * sub_hull(Set_of_Points all_points, Point p1, Point p2);
void sub_hull_new(std::vector<Point *> &input_points, Point* p1, Point* p2, std::list<Point *> &convex_hull);

#endif
