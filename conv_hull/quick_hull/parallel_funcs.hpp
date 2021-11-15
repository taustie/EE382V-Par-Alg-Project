#ifndef PARALLEL_FUNCS
#define PARALLEL_FUNCS

Point * neutral_point();
Point * myMax(Point * a, Point * b);
Point * myMin(Point * a, Point * b);
void get_max_min_points_parallel(std::vector<Point *> &input_points, Point** max, Point** min);

typedef struct distance_info{
	long long int index;
	double line_dist;
	double p1_dist;
} Dist_Info;

Dist_Info neutral_distance();
Dist_Info Dist_Max_Compare(Dist_Info a, Dist_Info b);
void get_max_dist_parallel(std::vector<Point *> &left_points, Point** max, Point *p1, Point *p2, int threads);

void get_points_on_left_parallel(std::vector<Point *> &left_points, Point ** input_points, int size, Point* p1, Point* p2, int threads);
// void get_points_on_left_filter_parallel(std::vector<Point *> &left_points, Point ** input_points, int size, Point* p1, Point* p2);

#endif
