#ifndef VERIFY_CONV_HULL
#define VERIFY_CONV_HULL


void find_min_and_sort(std::list<Point *> &output_points);
bool compare_angles (Point * first, Point * second);
bool verify_ccw_angles(std::vector<Point *> &input_points, std::list<Point *> &output_copy);
bool verify_convex_hull(std::vector<Point *> &input_points, std::list<Point *> &output_points);

#endif
