#ifndef UTILS_H
#define UTILS_H

int quad(pair<int, int> p);
int orientation(pair<int, int> a, pair<int, int> b, pair<int, int> c);
int calc_left_m(vector<pair<int, int>> a);
int calc_right_m(vector<pair<int, int>> a);
double dis(pair<int, int> p1, pair<int, int> p2);
void fill_vector(vector<pair<int, int>> src, vector<pair<int, int>> des);

#endif
