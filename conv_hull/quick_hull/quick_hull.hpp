#ifndef QUICK_HULL
#define QUICK_HULL

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

void quick_hull_new(std::vector<Point *> &input_points, std::list<Point *> &convex_hull);
void sub_hull_new(std::vector<Point *> &input_points, Point* p1, Point* p2, std::list<Point *> &convex_hull);

#endif
