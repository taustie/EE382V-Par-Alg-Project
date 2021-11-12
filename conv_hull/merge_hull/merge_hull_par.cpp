#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <chrono>
#include <random>
#include <algorithm>
using namespace std::chrono;
using namespace std;

pair<int, int> mid;

int find_quad(pair<int, int> p)
{
	if (p.first >= 0 && p.second >= 0)
    {
        return 1;
    }
	if (p.first <= 0 && p.second >= 0)
    {
        return 2;
    }

	if (p.first <= 0 && p.second <= 0)
    {
        return 3;
    }
	return 4;
}

// Checks whether the line is crossing the polygon
int orientation(pair<int, int> a, pair<int, int> b, pair<int, int> c)
{
	int res = (b.second-a.second)*(c.first-b.first) - (c.second-b.second)*(b.first-a.first);

	if (res == 0)
    {
        return 0;
    }
	else if (res > 0)
    {
        return 1;
    }
	return -1;
}

// compare function for sorting
bool compare_points(pair<int, int> p1, pair<int, int> q1)
{
    int new_p1_x = p1.first - mid.first, new_p1_y = p1.second - mid.second;
    int new_p2_x = q1.first - mid.first, new_p2_y = q1.second - mid.second;
	pair<int, int> p = make_pair(new_p1_x, new_p1_y);
	pair<int, int> q = make_pair(new_p2_x, new_p2_y);

	int one = find_quad(p);
	int two = find_quad(q);

	if (one != two)
    {
        return (one < two);
    } else {
        return (q.second*p.first > p.second*q.first);
    }

}

vector<pair<int, int>> merge_hulls(vector<pair<int, int>> a, vector<pair<int, int>> b)
{
	// n1 -> number of points in polygon a
	// n2 -> number of points in polygon b
	int n1 = a.size(), n2 = b.size();

	int right_m = 0, left_m = 0;
	for (int i=1; i<n1; i++)
    {
        if (a[i].first > a[right_m].first)
        {
            right_m = i;
        }
    }
		

	// ib -> leftmost point of b
	for (int i=1; i<n2; i++)
    {
		if (b[i].first < b[left_m].first)
        {
			left_m=i;
        }
    }
	// finding the upper tangent
	
    int uppera, upperb, lowera, lowerb;
	
    
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            int curr_l = right_m, curr_r = left_m;
            bool crosses = 1;
            while (crosses)
            {
                crosses = 0;
                while (orientation(b[curr_r], a[curr_l], a[(curr_l+1)%n1]) >=0)
                {
                    curr_l = (curr_l + 1) % n1;
                }
                    

                while (orientation(a[curr_l], b[curr_r], b[(n2+curr_r-1)%n2]) <=0)
                {
                    curr_r = (n2+curr_r-1)%n2;
                    crosses = 1;
                }
            }
            uppera = curr_l;
            upperb = curr_r;
        }

        #pragma omp section
        {
            int curr_l = right_m, curr_r = left_m;
            bool crosses1 = 1;
            while (crosses1)
            {
                crosses1 = 0;
                while (orientation(a[curr_l], b[curr_r], b[(curr_r+1)%n2])>=0)
                {
                    curr_r=(curr_r+1)%n2;
                }
                while (orientation(b[curr_r], a[curr_l], a[(n1+curr_l-1)%n1])<=0)
                {
                    curr_l=(n1+curr_l-1)%n1;
                    crosses1 = 1;
                }
            }
            lowera = curr_l;
            lowerb = curr_r;
        }
    }

	vector<pair<int, int>> ret;

	int ind = uppera;
	ret.push_back(a[uppera]);
	while (ind != lowera)
	{
		ind = (ind+1)%n1;
		ret.push_back(a[ind]);
	}

	ind = lowerb;
	ret.push_back(b[lowerb]);
	while (ind != upperb)
	{
		ind = (ind+1)%n2;
		ret.push_back(b[ind]);
	}
	return ret;

}

vector<pair<int, int>> brute_force(vector<pair<int, int>> points)
{
	vector<pair<int, int>> convex_hull;

	for (int p1=0; p1<points.size(); p1++)
	{
		for (int p2=p1+1; p2<points.size(); p2++)
		{
			int x1 = points[p1].first, x2 = points[p2].first;
			int y1 = points[p1].second, y2 = points[p2].second;

			int a = y1-y2;
			int b = x2-x1;
			int c = x1*y2-y1*x2;
			int pos = 0, neg = 0;
			for (int k = 0; k < points.size(); k++)
			{
                if (a*points[k].first+b*points[k].second+c >= 0)
                {
                    pos++;
                }
				if (a*points[k].first+b*points[k].second+c <= 0)
                {
                    neg++;
                }
			}
			if (pos == points.size() || neg == points.size())
			{
				convex_hull.push_back(points[p1]);
				convex_hull.push_back(points[p2]);
			}
		}
	}

	// Geeks for geeks
	pair<int, int> mid = {0, 0};
	int len = convex_hull.size();
	for (int i=0; i < len; i++)
	{
		mid.first += convex_hull[i].first;
		mid.second += convex_hull[i].second;
		convex_hull[i].first *= len;
		convex_hull[i].second *= len;
	}
	sort(convex_hull.begin(), convex_hull.end(), compare_points);

    #pragma omp parallel
    {
        #pragma omp for
        for (int i=0; i < len; i++)
        {
            int new_x = convex_hull[i].first/len;
            int new_y = convex_hull[i].second/len;
            convex_hull[i] = make_pair(new_x, new_y);
        }
    }

	return convex_hull;
}

vector<pair<int, int>> convex_hull_par(vector<pair<int, int>> points)
{
	if (points.size() <= 5)
    {
        return brute_force(points);
    }

    vector<pair<int, int>> left, right;
    int mid = points.size()/2;
    #pragma omp parallel
    {
        vector<pair<int, int>> left_temp;
        #pragma omp for nowait schedule(static)
        for (int i = 0; i < mid; i++) {
            left_temp.push_back(points[i]);
        }
        #pragma omp for schedule(static) ordered
        for(int i=0; i < omp_get_num_threads(); i++) {
            #pragma omp ordered
            left.insert(left.end(), left_temp.begin(), left_temp.end());
        }
    }

    #pragma omp parallel
    {
        vector<pair<int, int>> right_temp;
        #pragma omp for nowait schedule(static)
        for (int i = mid; i < points.size(); i++) {
            right_temp.push_back(points[i]);
        }
        #pragma omp for schedule(static) ordered
        for(int i=0; i < omp_get_num_threads(); i++) {
            #pragma omp ordered
            right.insert(right.end(), right_temp.begin(), right_temp.end());
        }
    }

    vector<pair<int, int>>left_hull;
    vector<pair<int, int>>right_hull;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            left_hull = convex_hull_par(left);
        }

        #pragma omp section
        {
            right_hull = convex_hull_par(right);
        }

    }

	return merge_hulls(left_hull, right_hull);
}
