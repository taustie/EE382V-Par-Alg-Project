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

// stores the centre of polygon (It is made
// global because it is used in compare function)
pair<int, int> mid;

// determines the quadrant of a point
// (used in compare())
int quad(pair<int, int> p)
{
	if (p.first >= 0 && p.second >= 0)
		return 1;
	if (p.first <= 0 && p.second >= 0)
		return 2;
	if (p.first <= 0 && p.second <= 0)
		return 3;
	return 4;
}

// Checks whether the line is crossing the polygon
int orientation(pair<int, int> a, pair<int, int> b,
				pair<int, int> c)
{
	int res = (b.second-a.second)*(c.first-b.first) -
			(c.second-b.second)*(b.first-a.first);

	if (res == 0)
		return 0;
	if (res > 0)
		return 1;
	return -1;
}

// compare function for sorting
bool compare(pair<int, int> p1, pair<int, int> q1)
{
	pair<int, int> p = make_pair(p1.first - mid.first,
								p1.second - mid.second);
	pair<int, int> q = make_pair(q1.first - mid.first,
								q1.second - mid.second);

	int one = quad(p);
	int two = quad(q);

	if (one != two)
		return (one < two);
	return (p.second*q.first < q.second*p.first);
}

// Finds upper tangent of two polygons 'a' and 'b'
// represented as two vectors.
vector<pair<int, int>> merger(vector<pair<int, int> > a,
							vector<pair<int, int> > b)
{
	// n1 -> number of points in polygon a
	// n2 -> number of points in polygon b
	int n1 = a.size(), n2 = b.size();

	int ia = 0, ib = 0;
	for (int i=1; i<n1; i++)
		if (a[i].first > a[ia].first)
			ia = i;

	// ib -> leftmost point of b
	for (int i=1; i<n2; i++)
		if (b[i].first < b[ib].first)
			ib=i;

	// finding the upper tangent
	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)
	{
		done = 1;
		while (orientation(b[indb], a[inda], a[(inda+1)%n1]) >0)
		{
			inda = (inda + 1) % n1;
			cout << "while 1" << endl;
			// cout << inda << endl;
			cout << b[inda].first;
			cout << b[inda].second << endl;
		}
			
		while (orientation(a[inda], b[indb], b[(n2+indb-1)%n2]) <0)
		{
			indb = (n2+indb-1)%n2;
			done = 0;
			cout << "while 2" << endl;
			cout << b[indb].first;
			cout << b[indb].second << endl;
		}
	}

	int uppera = inda, upperb = indb;
	inda = ia, indb=ib;
	done = 0;
	int g = 0;




	while (!done)//finding the lower tangent
	{
		done = 1;
		while (orientation(a[inda], b[indb], b[(indb+1)%n2])>0)
		{
			indb=(indb+1)%n2;
			cout << "while 3" << endl;
			cout << b[indb].first;
			cout << b[indb].second << endl;
		}
			

		while (orientation(b[indb], a[inda], a[(n1+inda-1)%n1])<0)
		{
			inda=(n1+inda-1)%n1;
			done=0;
			cout << "while 4" << endl;
			cout << b[inda].first;
			cout << b[inda].second << endl;
		}
	}

	int lowera = inda, lowerb = indb;
	vector<pair<int, int>> ret;

	//ret contains the convex hull after merging the two convex hulls
	//with the points sorted in anti-clockwise order
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

// Brute force algorithm to find convex hull for a set
// of less than 6 points
vector<pair<int, int>> bruteHull(vector<pair<int, int>> a)
{
	// Take any pair of points from the set and check
	// whether it is the edge of the convex hull or not.
	// if all the remaining points are on the same side
	// of the line then the line is the edge of convex
	// hull otherwise not
	set<pair<int, int> >s;

	for (int i=0; i<a.size(); i++)
	{
		for (int j=i+1; j<a.size(); j++)
		{
			int x1 = a[i].first, x2 = a[j].first;
			int y1 = a[i].second, y2 = a[j].second;

			int a1 = y1-y2;
			int b1 = x2-x1;
			int c1 = x1*y2-y1*x2;
			int pos = 0, neg = 0;
			for (int k=0; k<a.size(); k++)
			{
				if (a1*a[k].first+b1*a[k].second+c1 <= 0)
					neg++;
				if (a1*a[k].first+b1*a[k].second+c1 >= 0)
					pos++;
			}
			if (pos == a.size() || neg == a.size())
			{
				s.insert(a[i]);
				s.insert(a[j]);
			}
		}
	}

	vector<pair<int, int>>ret;
	for (auto e:s)
		ret.push_back(e);

	// Sorting the points in the anti-clockwise order
	mid = {0, 0};
	int n = ret.size();
	for (int i=0; i<n; i++)
	{
		mid.first += ret[i].first;
		mid.second += ret[i].second;
		ret[i].first *= n;
		ret[i].second *= n;
	}
	sort(ret.begin(), ret.end(), compare);
	for (int i=0; i<n; i++)
		ret[i] = make_pair(ret[i].first/n, ret[i].second/n);

	return ret;
}

// Returns the convex hull for the given set of points
vector<pair<int, int>> divide(vector<pair<int, int>> a)
{
	// If the number of points is less than 6 then the
	// function uses the brute algorithm to find the
	// convex hull
	if (a.size() <= 5)
		return bruteHull(a);

	// left contains the left half points
	// right contains the right half points
	vector<pair<int, int>>left, right;
	for (int i=0; i<a.size()/2; i++)
		left.push_back(a[i]);
	for (int i=a.size()/2; i<a.size(); i++)
		right.push_back(a[i]);

	// convex hull for the left and right sets
	vector<pair<int, int>>left_hull = divide(left);
	vector<pair<int, int>>right_hull = divide(right);

	// merging the convex hulls
	return merger(left_hull, right_hull);
}









// Driver code
int main( int argc, char* argv[] )
{ 
    std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> uni(1,100);
    vector<pair<int, int>> points;
	for(int i = 0; i < 10000; i++)
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
    vector<pair<int, int>> final_hull = divide(points);
    auto end = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(end - start);
	cout << duration.count() << endl;
    return 0;
}
