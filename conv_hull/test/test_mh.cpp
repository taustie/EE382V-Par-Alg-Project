#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
//#include <bits/stdc++.h>
#include <chrono>
#include <random>
#include <algorithm>
using namespace std::chrono;
using namespace std;

#include "quick_hull.hpp"
#include "verify_conv_hull.hpp"
#include "test_config.hpp"

vector<pair<int, int>> convex_hull_par(vector<pair<int, int>> points);

// Driver code
int main( int argc, char* argv[] )
{
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100,100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 200; i++)
    {
		auto random_integer = uni(rng);
		int x = (int)random_integer;
		random_integer = uni(rng);
		int y = (int)random_integer;
		points.push_back(make_pair(x, y));
	}
    //parallel sort
    sort(points.begin(), points.end());
    auto start = high_resolution_clock::now();
    vector<pair<int, int>> final_hull = convex_hull_par(points);
    auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << duration.count() << endl;
    return 0;
}
