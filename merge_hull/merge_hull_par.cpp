#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include<bits/stdc++.h>
#include <chrono>
#include <random>
using namespace std::chrono;
using namespace std;

pair<int, int> mid;

int find_quadrant(int x, int y)
{
    if (x >= 0 && y >= 0)
    {
        return 1;
    }
    else if (x <= 0 && y >= 0)
    {
        return 2;
    }
    else if (x <= 0 && y <= 0)
    {
        return 3;
    }
    else 
    {
        return 4;
    }
}

int orientation(pair<int, int> a, pair<int, int> b, pair<int, int> c)
{
    int res = (b.second-a.second)*(c.first-b.first) - (c.second-b.second)*(b.first-a.first);
    
    if (res == 0)
    {
        return 0;
    }
    if (res > 0)
    {
        return 1;
    }
    return -1;
}

bool compare_points(pair<int, int> point1, pair<int, int> point2)
{
    int q1 = find_quadrant(point1.first - mid.first, point1.second - mid.second);
    int q2 = find_quadrant(point2.first - mid.first, point2.second - mid.second);
    if (q1 != q2) 
    {
        return (q1 < q2);
    }
        
    return (point1.second * point2.first < point1.second * point2.first);
}

vector<pair<int, int>> merge_hull(vector<pair<int, int>> left, vector<pair<int, int>> right)
{

    int l_size = left.size();
    int r_size = right.size();

    int right_m = 0;
    int left_m = 0;
    //parallel reduce

    int left_max = left[right_m].first;
    int right_max = right[left_m].first;

    #pragma omp parallel for reduction(max:left_max)
    for (int i = 1; i < l_size; i++)
    {
        if (left[i].first > left[right_m].first) 
        {
            left_max = left[i].first;
            right_m = i;
        }
    }

    //parallel reduce
    #pragma omp parallel for reduction(max:right_max)
    for (int i = 1; i < r_size; i++)
    {
        if (right[i].first < right[left_m].first) 
        {
            right_max = right[i].first;
            left_m = i;
        }
    }

    int left_c = right_m;
    int right_c = left_m;
    bool done = 0;
    while (!done)
    {
        done = 1;
        while (orientation(right[right_c], left[left_c], left[(left_c+1)%l_size]) >=0)
        {
            left_c=(left_c+1)%l_size;
        }
            
        while (orientation(left[left_c], right[right_c], right[(r_size+right_c-1)%r_size]) <=0)
        {
            right_c = (r_size+right_c-1)%r_size;
            done = 0;
        }
    }

    int uppera = left_c, upperb = right_c;
    left_c = right_m;
    right_c =left_m;
    done = 0;
    int g = 0;
    while (!done)
    {
        done = 1;
        while (orientation(left[left_c], right[right_c], right[(right_c+1)%r_size])>=0)
        {
            right_c=(right_c+1)%r_size;
        }
        while (orientation(right[right_c], left[left_c], left[(l_size+left_c-1)%l_size])<=0)
        {
            left_c=(l_size+left_c-1)%l_size;
            done=0;
        }
    }
    int lowera = left_c, lowerb = right_c;

    vector<pair<int, int>> combined_hull;
    int ind = uppera;
    combined_hull.push_back(left[uppera]);
    while (ind != lowera)
    {
        ind = (ind+1)%l_size;
        combined_hull.push_back(left[ind]);
    }
    ind = lowerb;
    combined_hull.push_back(right[lowerb]);
    while (ind != upperb)
    {
        ind = (ind+1)%r_size;
        combined_hull.push_back(right[ind]);
    }

    
    return combined_hull;
}

vector<pair<int, int>> brute_force(vector<pair<int, int>> points, int start, int end)
{
    vector<pair<int, int>> selected_points;

    #pragma omp parallel
    {
        vector<pair<int, int>> temp;
        #pragma omp for nowait schedule(static)
        for(int i = start; i <= end; i++) { 
            temp.push_back(points[i]);
        }
        #pragma omp for schedule(static) ordered
        for(int i=0; i<omp_get_num_threads(); i++) {
            #pragma omp ordered
            selected_points.insert(selected_points.end(), temp.begin(), temp.end());
        }
    }


    vector<pair<int, int>> sorted_hull;

    for (int i=0; i < selected_points.size(); i++)
    {
        for (int j=i+1; j < selected_points.size(); j++)
        {
            int px1 = selected_points[i].first;
            int py1 = selected_points[i].second;
            int px2 = selected_points[j].first;
            int py2 = selected_points[j].second;
            int a = py1 - py2;
            int b = px2 - px1;
            int c = (px1 * py2) - (py1 * px2);

            vector<int> pos_points;
            vector<int> neg_points;
            for (int k = 0; k < selected_points.size(); k++)
            {
                int side = (a*selected_points[k].first) + (b*selected_points[k].second) + c;
                if (side <= 0)
                {
                    neg_points.push_back(1);
                }
                    
                if (side >= 0)
                {
                    pos_points.push_back(1);
                }
            }
            if (pos_points.size() == selected_points.size() || neg_points.size() == selected_points.size())
            {
                sorted_hull.push_back(selected_points[i]);
                sorted_hull.push_back(selected_points[j]);
            }
        }
    }
        
    mid = {0, 0};
    for (int i = 0; i < sorted_hull.size(); i++)
    {
        mid.first += sorted_hull[i].first;
        mid.second += sorted_hull[i].second;
        sorted_hull[i].first *= sorted_hull.size();
        sorted_hull[i].second *= sorted_hull.size();
    }

    // parallel sort
    sort(sorted_hull.begin(), sorted_hull.end(), compare_points);

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < sorted_hull.size(); i++)
        {
            sorted_hull[i] = make_pair(sorted_hull[i].first/sorted_hull.size(), sorted_hull[i].second/sorted_hull.size());
        }
    }

    
    return sorted_hull;
}

// vector<pair<int, int>> find_convex_hull_seq(vector<pair<int, int>> points, int start, int end)
// {
//     if (start <= end && end - start + 1 <= 5)
//     {
//         return brute_force(points, start, end);
//     }

//     int mid_ind = (start + end) / 2;

//     vector<pair<int, int>>right_hull = find_convex_hull_seq(points, mid_ind, end);
//     vector<pair<int, int>>left_hull = find_convex_hull_seq(points, start, mid_ind-1);   

//     vector<pair<int, int>> ans =  merge_hull(left_hull, right_hull);
//     return ans;
// }

vector<pair<int, int>> find_convex_hull_par(vector<pair<int, int>> points, int start, int end, int threads)
{
    // if(threads == 1)
    // {
    //     return find_convex_hull_seq(points, start, end);
    // }

    if (end - start + 1 <= 5)
    {
        return brute_force(points, start, end);
    }

    int mid_ind = (start + end) /2 ;

    vector<pair<int, int>>right_hull;
    vector<pair<int, int>>left_hull;
    #pragma omp parallel sections
	{
		#pragma omp section
            right_hull = find_convex_hull_par(points, mid_ind, end, threads/2);

        #pragma omp section
            left_hull = find_convex_hull_par(points, start, mid_ind-1, threads - (threads/2));    
    }

    vector<pair<int, int>> ans =  merge_hull(left_hull, right_hull);
    return ans;
}

// Driver code
int main( int argc, char* argv[] )
{
    int threads = omp_get_max_threads();
	omp_set_num_threads(threads);
    omp_set_nested(1);
    
    
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100,100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 10; i++)
    {
		auto random_integer = uni(rng);
		int x = (int)random_integer;
		random_integer = uni(rng);
		int y = (int)random_integer;
		points.push_back(make_pair(x, y));
	}
    
    sort(points.begin(), points.end());
    auto start = high_resolution_clock::now();
    vector<pair<int, int>> final_hull = find_convex_hull_par(points, 0, points.size() - 1, threads);
    auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << duration.count() << endl;
    return 0;
}
