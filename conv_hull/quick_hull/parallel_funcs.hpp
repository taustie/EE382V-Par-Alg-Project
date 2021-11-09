#ifndef PARALLEL_FUNCS
#define PARALLEL_FUNCS

typedef struct distance_info{
	long long int index;
	double line_dist;
	double p1_dist;
} Dist_Info;

Point * get_empty();
Point * myMax(Point * a, Point * b);
void get_max_min_points_parallel(std::vector<Point *> &input_points, Point** max, Point** min);

Dist_Info Dist_Max_Compare(Dist_Info a, Dist_Info b);
Dist_Info neutral_distance();
void get_max_dist_parallel(std::vector<Point *> &left_points, Point** max, Point *p1, Point *p2);

#endif
