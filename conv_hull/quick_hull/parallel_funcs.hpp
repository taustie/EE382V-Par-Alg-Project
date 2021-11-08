#ifndef PARALLEL_FUNCS
#define PARALLEL_FUNCS

Point * get_empty();
Point * myMax(Point * a, Point * b);
void get_max_min_points_parallel(std::vector<Point *> &input_points, Point** max, Point** min);

#endif
