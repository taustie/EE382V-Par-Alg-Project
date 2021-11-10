#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include<bits/stdc++.h>
#include <chrono>
#include <random>
using namespace std::chrono;
using namespace std;

// pair<int, int> mid;

// int find_quadrant(int x, int y)
// {
//     if (x >= 0 && y >= 0)
//     {
//         return 1;
//     }
//     else if (x <= 0 && y >= 0)
//     {
//         return 2;
//     }
//     else if (x <= 0 && y <= 0)
//     {
//         return 3;
//     }
//     else 
//     {
//         return 4;
//     }
// }

// int orientation(pair<int, int> a, pair<int, int> b, pair<int, int> c)
// {
//     int res = (b.second-a.second)*(c.first-b.first) - (c.second-b.second)*(b.first-a.first);
    
//     if (res == 0)
//     {
//         return 0;
//     }
//     if (res > 0)
//     {
//         return 1;
//     }
//     return -1;
// }

// bool compare_points(pair<int, int> point1, pair<int, int> point2)
// {
//     int q1 = find_quadrant(point1.first - mid.first, point1.second - mid.second);
//     int q2 = find_quadrant(point2.first - mid.first, point2.second - mid.second);
//     if (q1 != q2) 
//     {
//         return (q1 < q2);
//     }
        
//     return (point1.second * point2.first < point1.second * point2.first);
// }

// vector<pair<int, int>> merge_hull(vector<pair<int, int>> left, vector<pair<int, int>> right)
// {

//     int l_size = left.size();
//     int r_size = right.size();

//     int right_m = 0;
//     int left_m = 0;

//     int left_max = left[right_m].first;
//     int right_max = right[left_m].first;

//     for (int i = 1; i < l_size; i++)
//     {
//         if (left[i].first > left[right_m].first) 
//         {
//             left_max = left[i].first;
//             right_m = i;
//         }
//     }

//     for (int i = 1; i < r_size; i++)
//     {
//         if (right[i].first < right[left_m].first) 
//         {
//             right_max = right[i].first;
//             left_m = i;
//         }
//     }

//     int left_c = right_m;
//     int right_c = left_m;
//     bool done = 0;
//     while (!done)
//     {
//         done = 1;
//         while (orientation(right[right_c], left[left_c], left[(left_c+1)%l_size]) >=0)
//         {
//             left_c=(left_c+1)%l_size;
//         }
            
//         while (orientation(left[left_c], right[right_c], right[(r_size+right_c-1)%r_size]) <=0)
//         {
//             right_c = (r_size+right_c-1)%r_size;
//             done = 0;
//         }
//     }

//     int uppera = left_c, upperb = right_c;
//     left_c = right_m;
//     right_c =left_m;
//     done = 0;
//     int g = 0;
//     while (!done)
//     {
//         done = 1;
//         while (orientation(left[left_c], right[right_c], right[(right_c+1)%r_size])>=0)
//         {
//             right_c=(right_c+1)%r_size;
//         }
//         while (orientation(right[right_c], left[left_c], left[(l_size+left_c-1)%l_size])<=0)
//         {
//             left_c=(l_size+left_c-1)%l_size;
//             done=0;
//         }
//     }
//     int lowera = left_c, lowerb = right_c;

//     vector<pair<int, int>> combined_hull;
//     int ind = uppera;
//     combined_hull.push_back(left[uppera]);
//     while (ind != lowera)
//     {
//         ind = (ind+1)%l_size;
//         combined_hull.push_back(left[ind]);
//     }
//     ind = lowerb;
//     combined_hull.push_back(right[lowerb]);
//     while (ind != upperb)
//     {
//         ind = (ind+1)%r_size;
//         combined_hull.push_back(right[ind]);
//     }

    
//     return combined_hull;
// }

// vector<pair<int, int>> brute_force(vector<pair<int, int>> points, int start, int end)
// {
//     vector<pair<int, int>> selected_points;
//     for(int i = start; i <= end; i++) { 
//         selected_points.push_back(points[i]);
//     }


//     vector<pair<int, int>> sorted_hull;

//     for (int i=0; i < selected_points.size(); i++)
//     {
//         for (int j=i+1; j < selected_points.size(); j++)
//         {
//             int px1 = selected_points[i].first;
//             int py1 = selected_points[i].second;
//             int px2 = selected_points[j].first;
//             int py2 = selected_points[j].second;
//             int a = py1 - py2;
//             int b = px2 - px1;
//             int c = (px1 * py2) - (py1 * px2);

//             vector<int> pos_points;
//             vector<int> neg_points;
//             for (int k = 0; k < selected_points.size(); k++)
//             {
//                 int side = (a*selected_points[k].first) + (b*selected_points[k].second) + c;
//                 if (side <= 0)
//                 {
//                     neg_points.push_back(1);
//                 }
                    
//                 if (side >= 0)
//                 {
//                     pos_points.push_back(1);
//                 }
//             }
//             if (pos_points.size() == selected_points.size() || neg_points.size() == selected_points.size())
//             {
//                 sorted_hull.push_back(selected_points[i]);
//                 sorted_hull.push_back(selected_points[j]);
//             }
//         }
//     }
        
//     mid = {0, 0};
//     for (int i = 0; i < sorted_hull.size(); i++)
//     {
//         mid.first += sorted_hull[i].first;
//         mid.second += sorted_hull[i].second;
//         sorted_hull[i].first *= sorted_hull.size();
//         sorted_hull[i].second *= sorted_hull.size();
//     }

//     // parallel sort
//     sort(sorted_hull.begin(), sorted_hull.end(), compare_points);
//     for (int i = 0; i < sorted_hull.size(); i++)
//     {
//         sorted_hull[i] = make_pair(sorted_hull[i].first/sorted_hull.size(), sorted_hull[i].second/sorted_hull.size());
//     }

    
//     return sorted_hull;
// }

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











// A divide and conquer program to find convex
// hull of a given set of points.
#include<bits/stdc++.h>
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









// Driver code
int main( int argc, char* argv[] )
{ 
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(-100,100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 210; i++)
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
